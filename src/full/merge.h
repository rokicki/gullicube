// Merge many independent orbit solutions into one move sequence, choosing the
// interleaving greedily to maximise cancellation at the seams.
//
// Each algorithm instance carries its moves (brobdicube mv: face, depth, twist)
// and, for each orbit component it touches, a bitmap of the positions it moves.
// Within a component, instance j must stay after an earlier instance i when
// their bitmaps overlap (otherwise they commute).  Instances of different
// components act on disjoint pieces and may be interleaved freely.
//
// Moves are kept as physical layers: axis (U/D=0, L/R=1, F/B=2), position
// 1..N along the axis measured from the U/L/F side, and a quarter-turn amount
// in the U/L/F turning direction.  Same-axis layers commute, so the sequence
// is a stack of axis groups; cost = number of nonzero layer turns.
#pragma once
#include "mv.h"
#include <algorithm>
#include <cstdint>
#include <vector>
#include <atomic>
#include <chrono>
#include <unordered_map>
#include <queue>
#include <tuple>
#include <thread>
#include <mutex>

struct PLayer { int8_t axis; int16_t pos; int8_t amt; };

struct MergeInst {
  std::vector<PLayer> moves;
  std::vector<std::pair<int, uint64_t>> deps;  // (component, moved-position bitmap)
};

class Merger {
 public:
  explicit Merger(int N) : N(N) {}

  PLayer toLayer(const mv &m) const {
    static const int AX[6] = {0, 1, 2, 1, 2, 0}, POS[6] = {1, 1, 1, 0, 0, 0};  // ULFRBD
    PLayer p;
    p.axis = AX[m.face];
    bool pos = POS[m.face];
    p.pos = pos ? m.dep : N + 1 - m.dep;
    p.amt = pos ? m.twist : (4 - m.twist) & 3;
    return p;
  }
  mv toMove(const PLayer &p) const {
    static const short FACE[3] = {0, 1, 2};  // U, L, F
    return {p.pos, FACE[p.axis], (short)p.amt};
  }

  int add(const std::vector<mv> &moves, std::vector<std::pair<int, uint64_t>> deps) {
    MergeInst in;
    for (auto &m : moves) in.moves.push_back(toLayer(m));
    in.deps = std::move(deps);
    insts.push_back(std::move(in));
    return insts.size() - 1;
  }

  // Greedy merge.  Returns the merged moves; rawMoves / mergedMoves are filled.
  // Time budget for the greedy merge (seconds, 0 = none): the number of
  // candidates evaluated per step adapts (halving when the projected time
  // exceeds the budget, doubling when well under, within [4, 256]).
  double budget = 0;
  std::vector<mv> merge(long &rawMoves, long &mergedMoves, std::vector<int> *order = nullptr) {
    if (order) order->clear();
    size_t n = insts.size();
    // dependency DAG within components (instances are added in solution order)
    std::vector<std::vector<int>> succ(n);
    std::vector<int> npred(n, 0);
    {
      std::vector<std::vector<std::pair<int, uint64_t>>> comp;  // per component: (inst, mask)
      for (size_t i = 0; i < n; i++)
        for (auto &d : insts[i].deps) {
          if ((int)comp.size() <= d.first) comp.resize(d.first + 1);
          for (auto &prev : comp[d.first])
            if (prev.second & d.second) {
              succ[prev.first].push_back(i);
              npred[i]++;
            }
          comp[d.first].push_back({(int)i, d.second});
        }
    }
    // Available instances are indexed by the (axis, layer, twist) of each move
    // in their first same-axis group: only an instance sharing a layer with
    // the tail's last group can save anything (opposite twists cancel, saving
    // 2; others merge, saving 1).  Each step looks in the cancelling buckets
    // first, then the merging ones, evaluating the exact saving for a few
    // members (a full cancellation can cascade); if nothing saves, the oldest
    // available instance is taken.  Near O(1) per step instead of a scan of
    // every available instance.
    rawMoves = 0;
    for (auto &in : insts) rawMoves += in.moves.size();
    int mergeLim = getenv("MERGELIM") ? atoi(getenv("MERGELIM")) : 256;  // candidates evaluated per step (16: 1.5% worse at 256^3)
    const auto tStart = std::chrono::steady_clock::now();
    auto bkey = [](int axis, int pos, int amt) { return ((uint64_t)axis << 32) | ((uint64_t)pos << 2) | (uint64_t)amt; };
    std::unordered_map<uint64_t, std::vector<int>> bucket;
    std::vector<char> placed(n, 0);
    std::vector<int> fifo;
    size_t fifoHead = 0;
    auto makeAvail = [&](int i) {
      fifo.push_back(i);
      const auto &mv = insts[i].moves;
      for (size_t k = 0; k < mv.size() && mv[k].axis == mv[0].axis; k++) bucket[bkey(mv[k].axis, mv[k].pos, mv[k].amt)].push_back(i);
    };
    for (size_t i = 0; i < n; i++)
      if (!npred[i]) makeAvail(i);
    tail.clear();
    cost = 0;
    long steps = 0;
    size_t nPlaced = 0;
    while (nPlaced < n) {
      if ((++steps & 1023) == 0) {
        if (stop && stop->load(std::memory_order_relaxed)) break;
        if (budget > 0) {  // pace: projected total time against the budget
          double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count();
          double proj = el * n / std::max<size_t>(1, nPlaced);
          if (proj > budget && mergeLim > 4) mergeLim /= 2;
          else if (proj < 0.5 * budget && mergeLim < 256) mergeLim *= 2;
        }
      }
      int best = -1, bestSave = 0, evaluated = 0;
      if (!tail.empty()) {
        const Group &g = tail.back();
        for (int pass = 0; pass < 2 && evaluated < mergeLim; pass++)  // pass 0: cancelling twists, 1: merging
          for (auto &l : g.layers) {
            for (int amt = 1; amt <= 3; amt++) {
              bool cancels = ((l.second + amt) & 3) == 0;
              if (cancels != (pass == 0)) continue;
              auto it = bucket.find(bkey(g.axis, l.first, amt));
              if (it == bucket.end()) continue;
              auto &v = it->second;
              for (size_t k = v.size(); k-- > 0 && evaluated < mergeLim;) {
                int i = v[k];
                if (placed[i]) { v[k] = v.back(); v.pop_back(); continue; }  // lazy removal
                evaluated++;
                int save = (int)insts[i].moves.size() - delta(insts[i].moves);
                if (save > bestSave) { bestSave = save; best = i; }
              }
            }
          }
      }
      if (best < 0) {  // nothing saves: the oldest available instance
        while (fifoHead < fifo.size() && placed[fifo[fifoHead]]) fifoHead++;
        if (fifoHead == fifo.size()) break;  // (cannot happen for a DAG)
        best = fifo[fifoHead];
      }
      placed[best] = 1;
      nPlaced++;
      if (order) order->push_back(best);
      for (auto &p : insts[best].moves) push(tail, cost, p);
      for (int s : succ[best])
        if (--npred[s] == 0) makeAvail(s);
    }
    mergedMoves = cost;
    std::vector<mv> out;
    for (auto &g : tail)
      for (auto &l : g.layers) out.push_back(toMove({g.axis, l.first, l.second}));
    return out;
  }

  // Merge beam: level by level (one instance placed per level), keeping the
  // `width` cheapest partial interleavings.  A state is the set of placed
  // instances, the per-instance count of unplaced predecessors, the last few
  // axis groups of its merged sequence (only those can interact with what
  // comes next) and its cost; equal (placed set, tail) states are merged.
  // Returns the best complete order's merged sequence.
  std::vector<mv> mergeBeam(int width, long &rawMoves, long &mergedMoves, std::vector<int> *order = nullptr) {
    const size_t n = insts.size();
    std::vector<std::vector<int>> succ(n);
    std::vector<uint16_t> npred0(n, 0);
    {
      std::vector<std::vector<std::pair<int, uint64_t>>> comp;
      for (size_t i = 0; i < n; i++)
        for (auto &d : insts[i].deps) {
          if ((int)comp.size() <= d.first) comp.resize(d.first + 1);
          for (auto &prev : comp[d.first])
            if (prev.second & d.second) { succ[prev.first].push_back(i); npred0[i]++; }
          comp[d.first].push_back({(int)i, d.second});
        }
    }
    rawMoves = 0;
    size_t maxLen = 1;
    for (auto &in : insts) { rawMoves += in.moves.size(); maxLen = std::max(maxLen, in.moves.size()); }
    const size_t K = maxLen + 1;  // tail groups that can interact with the next instance
    struct St { std::vector<uint64_t> placed; std::vector<uint16_t> npred; std::vector<Group> tail; long cost; int node; };
    struct Hist { int parent, inst; };
    std::vector<Hist> hist;
    const size_t nw = (n + 63) / 64;
    std::vector<St> cur(1);
    cur[0].placed.assign(nw, 0);
    cur[0].npred = npred0;
    cur[0].cost = 0;
    cur[0].node = -1;
    auto sigOf = [&](const St &s) {
      uint64_t h = 0x9E3779B97F4A7C15ULL;
      for (auto w : s.placed) { h ^= w; h *= 0xBF58476D1CE4E5B9ULL; h ^= h >> 31; }
      if (!s.tail.empty()) {
        const Group &g = s.tail.back();
        h ^= (uint64_t)(g.axis + 1) * 0x94D049BB133111EBULL;
        for (auto &l : g.layers) { h ^= (uint64_t)(l.first * 4 + l.second) + 0x632BE59BD9B4E019ULL; h *= 0xBF58476D1CE4E5B9ULL; }
      }
      return h;
    };
    for (size_t level = 0; level < n; level++) {
      if (stop && stop->load(std::memory_order_relaxed)) break;
      struct Cand { long cost; int s; int inst; };
      std::vector<Cand> cands;
      for (int si = 0; si < (int)cur.size(); si++) {
        const St &s = cur[si];
        for (size_t i = 0; i < n; i++) {
          if (s.placed[i / 64] >> (i % 64) & 1) continue;
          if (s.npred[i]) continue;
          std::vector<Group> t = s.tail;
          long c = 0;
          for (auto &p : insts[i].moves) push(t, c, p);
          cands.push_back({s.cost + c, si, (int)i});
        }
      }
      if (cands.empty()) break;
      std::sort(cands.begin(), cands.end(), [](const Cand &x, const Cand &y) { return x.cost < y.cost; });
      std::vector<St> next;
      std::unordered_map<uint64_t, int> seen;
      for (auto &c : cands) {
        if ((int)next.size() >= width) break;
        const St &s = cur[c.s];
        St ns;
        ns.placed = s.placed;
        ns.placed[c.inst / 64] |= 1ULL << (c.inst % 64);
        ns.tail = s.tail;
        long dc = 0;
        for (auto &p : insts[c.inst].moves) push(ns.tail, dc, p);
        if (ns.tail.size() > K) ns.tail.erase(ns.tail.begin(), ns.tail.end() - K);
        ns.cost = s.cost + dc;
        uint64_t sig = sigOf(ns);
        if (seen.count(sig)) continue;  // same placed set and tail, not cheaper (sorted)
        seen[sig] = 1;
        ns.npred = s.npred;
        for (int q : succ[c.inst]) ns.npred[q]--;
        hist.push_back({s.node, c.inst});
        ns.node = hist.size() - 1;
        next.push_back(std::move(ns));
      }
      cur.swap(next);
    }
    // best complete state
    int best = -1;
    for (int si = 0; si < (int)cur.size(); si++)
      if (best < 0 || cur[si].cost < cur[best].cost) best = si;
    std::vector<int> ord;
    for (int x = best >= 0 ? cur[best].node : -1; x >= 0; x = hist[x].parent) ord.push_back(hist[x].inst);
    std::reverse(ord.begin(), ord.end());
    if (ord.size() != n) { long r, m2; return merge(r, m2, order); }  // interrupted: fall back to greedy
    tail.clear();
    cost = 0;
    for (int i : ord)
      for (auto &p : insts[i].moves) push(tail, cost, p);
    mergedMoves = cost;
    if (order) *order = ord;
    std::vector<mv> out;
    for (auto &g : tail)
      for (auto &l : g.layers) out.push_back(toMove({g.axis, l.first, l.second}));
    return out;
  }


  // Pairwise merge (greedy edge, like building a Huffman tree): every
  // instance starts as a piece of its own; repeatedly the best seam anywhere
  // -- the end of one piece followed by the start of another -- joins the two
  // pieces, as long as the order constraints still allow some order.  Then
  // the pieces are laid out in an order the constraints allow.
  //  - Seam savings are exact from instance to instance (the last instance of
  //    one piece, the first of the other): same-axis groups meet at the seam,
  //    and when both cancel completely the next groups inward meet, and so
  //    on.  The final sequence is built and costed exactly.
  //  - Each piece end keeps its best seam (from the (axis, layer, twist)
  //    buckets of first groups, at most pairLim candidates scanned) in a heap;
  //    a stale best is recomputed when it reaches the top.
  //  - Pieces keep a topological order under the constraints (Pearce-Kelly).
  //    A seam P -> Q ("immediately before") is refused if a constraint path
  //    leads from Q to P, or from P to Q through another piece; both searches
  //    stay between P and Q in the order.  An accepted seam reorders only that
  //    region.  A search past pairVisitCap pieces refuses the seam (safe).
  int pairLim = 256;
  int pairVisitCap = 64;
  int pairWindow = 0;
  // restarts: with pairSeed != 0, equal seams are taken in a seeded random
  // order, and pairNoise (eighths of a move) of seeded noise is added to each
  // seam's priority; the joins themselves stay exact
  uint64_t pairSeed = 0;
  int pairNoise = 0;     // experiment: refuse seams between pieces further apart in the order (0: off)  // measured at 128^3: 64 costs 0.6% against 4096, at half the time
  // with `budget` (seconds, shared with merge()): every 1024 heap steps the
  // finish time is projected from the joins so far; over budget, the scan
  // limit and search cap halve (down to 8); over twice the budget once at
  // the minimum, joining stops and the pieces so far are laid out.
  long pairStatJoins = 0, pairStatRejected = 0, pairStatCapped = 0, pairStatVisits = 0, pairStatEdges = 0, pairStatScored = 0, pairStatEdgesOk = 0, pairStatRefDirect = 0, pairStatRefLong = 0, pairStatRefBetween = 0;
  std::vector<mv> mergePairwise(long &rawMoves, long &mergedMoves, std::vector<int> *order = nullptr) {
    if (order) order->clear();
    const int n = insts.size();
    rawMoves = 0;
    for (auto &in : insts) rawMoves += in.moves.size();
    // order constraints (as merge() builds them), both directions, per instance
    std::vector<std::vector<int>> outs(n), ins(n);
    {
      std::vector<std::vector<std::pair<int, uint64_t>>> comp;
      for (int i = 0; i < n; i++)
        for (auto &d : insts[i].deps) {
          if ((int)comp.size() <= d.first) comp.resize(d.first + 1);
          for (auto &prev : comp[d.first])
            if (prev.second & d.second) { outs[prev.first].push_back(i); ins[i].push_back(prev.first); }
          comp[d.first].push_back({i, d.second});
        }
    }
    // same-axis groups of each instance: moves [goff[g], goff[g+1]) for g in [gbeg[i], gbeg[i+1])
    std::vector<int> gbeg(n + 1), goff;
    for (int i = 0; i < n; i++) {
      gbeg[i] = goff.size();
      const auto &m = insts[i].moves;
      for (size_t k = 0; k < m.size(); k++)
        if (k == 0 || m[k].axis != m[k - 1].axis) goff.push_back(k);
    }
    gbeg[n] = goff.size();
    goff.push_back(0);  // sentinel (unused)
    auto gEnd = [&](int i, int g) { return g + 1 < gbeg[i + 1] ? goff[g + 1] : (int)insts[i].moves.size(); };
    // exact saving of j directly after i
    auto seamSave = [&](int i, int j) {
      const auto &a = insts[i].moves, &b = insts[j].moves;
      int save = 0;
      for (int gi = gbeg[i + 1] - 1, gj = gbeg[j]; gi >= gbeg[i] && gj < gbeg[j + 1]; gi--, gj++) {
        const int a0 = goff[gi], a1 = gEnd(i, gi), b0 = goff[gj], b1 = gEnd(j, gj);
        if (a[a0].axis != b[b0].axis) break;
        int cancelled = 0;
        for (int h = b0; h < b1; h++)
          for (int t = a0; t < a1; t++)
            if (a[t].pos == b[h].pos) {
              if (((a[t].amt + b[h].amt) & 3) == 0) { save += 2; cancelled++; }
              else save += 1;
              break;
            }
        if (cancelled != a1 - a0 || cancelled != b1 - b0) break;  // a group survives: the cascade stops
      }
      return save;
    };
    // pieces: union-find; at the root, head and tail instance and its place in the order
    std::vector<int> parent(n), phead(n), ptail(n), succ(n, -1), pred(n, -1), ord(n);
    std::vector<size_t> outsKept(n, 0), insKept(n, 0);
    for (int i = 0; i < n; i++) parent[i] = phead[i] = ptail[i] = ord[i] = i;  // constraints go forward in index
    auto find = [&](int x) {
      while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
      return x;
    };
    // buckets of piece heads by (axis, layer, twist) of their first group
    auto bkey = [&](int axis, int pos, int amt) { return ((size_t)axis * (N + 2) + pos) * 4 + amt; };
    std::vector<std::vector<int>> bucket(3 * (size_t)(N + 2) * 4);
    for (int j = 0; j < n; j++)
      if (gbeg[j] < gbeg[j + 1])
        for (int h = goff[gbeg[j]]; h < gEnd(j, gbeg[j]); h++) {
          const auto &p = insts[j].moves[h];
          bucket[bkey(p.axis, p.pos, p.amt)].push_back(j);
        }
    std::vector<std::vector<int>> rejected(n);  // per tail instance: heads refused by the constraints
    auto mix = [&](uint64_t x) {
      x += pairSeed * 0x9E3779B97F4A7C15ULL;
      x ^= x >> 31; x *= 0xBF58476D1CE4E5B9ULL; x ^= x >> 29; x *= 0x94D049BB133111EBULL; x ^= x >> 32;
      return x;
    };
    auto tieJ = [&](int j) -> uint64_t { return j < 0 ? ~0ULL : pairSeed ? mix(j) : (uint64_t)j; };
    // heap priority of a seam: 8 x saving, plus seeded noise
    auto prio = [&](int s, int i, int j) {
      return 8 * s + (pairSeed && pairNoise ? (int)(mix(((uint64_t)i << 32) ^ (uint32_t)j ^ 0x5555) % (pairNoise + 1)) : 0);
    };
    auto bestFrom = [&](int i) {
      int best = 0, bj = -1, evaluated = 0;
      if (gbeg[i] == gbeg[i + 1]) return std::make_pair(0, -1);
      const auto &a = insts[i].moves;
      const int ri = find(i), g = gbeg[i + 1] - 1;
      for (int pass = 0; pass < 2 && evaluated < pairLim; pass++)  // pass 0: cancelling twists first
        for (int t = goff[g]; t < gEnd(i, g) && evaluated < pairLim; t++)
          for (int amt = 1; amt <= 3 && evaluated < pairLim; amt++) {
            if ((((a[t].amt + amt) & 3) == 0) != (pass == 0)) continue;
            auto &v = bucket[bkey(a[t].axis, a[t].pos, amt)];
            for (size_t k = v.size(); k-- > 0 && evaluated < pairLim;) {
              const int j = v[k];
              if (pred[j] >= 0) { v[k] = v.back(); v.pop_back(); continue; }  // no longer a head
              if (find(j) == ri) continue;
              if (std::find(rejected[i].begin(), rejected[i].end(), j) != rejected[i].end()) continue;
              evaluated++;
              pairStatScored++;
              const int s = seamSave(i, j);
              if (s > best || (s == best && s > 0 && tieJ(j) < tieJ(bj))) { best = s; bj = j; }
            }
          }
      return std::make_pair(best, bj);
    };
    // amortized compaction of a root's edge list: roots, no self, no repeats
    std::vector<int> stamp(n, 0);
    int gen = 0;
    auto compact = [&](std::vector<int> &L, size_t &kept, int self) {
      if (L.size() <= 2 * kept + 16) return;
      gen++;
      size_t w = 0;
      for (int x : L) {
        int q = find(x);
        if (q == self || stamp[q] == gen) continue;
        stamp[q] = gen;
        L[w++] = q;
      }
      L.resize(w);
      kept = w;
    };
    // seam P -> Q (roots): 0 refused, 1 accepted (and the order updated), -1 cap
    std::vector<int> Fset, Bset, stack, slots;
    auto tryJoin = [&](int P, int Q) {
      const int lo = std::min(ord[P], ord[Q]), hi = std::max(ord[P], ord[Q]);
      if (pairWindow > 0 && hi - lo > pairWindow) return 0;
      const int L = ord[P] < ord[Q] ? P : Q, Hn = L == P ? Q : P;
      int visits = 0;
      long edgesHere = 0;
      // forward from L inside (lo, hi): reaching Hn is a contradiction, except
      // the direct edge P -> Q
      gen++;
      Fset.clear();
      stack.assign(1, L);
      stamp[L] = gen;
      while (!stack.empty()) {
        int r = stack.back();
        stack.pop_back();
        if (++visits > pairVisitCap) return -1;
        pairStatEdges += outs[r].size();
        edgesHere += outs[r].size();
        for (int x : outs[r]) {
          int q = find(x);
          if (q == r || stamp[q] == gen) continue;
          if (q == Hn) {
            if (r == P && L == P) continue;  // direct P -> Q: fine
            if (L == P) pairStatRefBetween++;       // P ~> X ~> Q
            else if (r == Q) pairStatRefDirect++;   // Q -> P directly
            else pairStatRefLong++;                 // Q ~> X ~> P
            return 0;
          }
          if (ord[q] > hi) continue;
          stamp[q] = gen;
          Fset.push_back(q);
          stack.push_back(q);
        }
      }
      // backward from Hn inside (lo, hi)
      Bset.clear();
      stack.assign(1, Hn);
      stamp[Hn] = gen;
      while (!stack.empty()) {
        int r = stack.back();
        stack.pop_back();
        if (++visits > pairVisitCap) return -1;
        pairStatEdges += ins[r].size();
        edgesHere += ins[r].size();
        for (int x : ins[r]) {
          int q = find(x);
          if (q == r || q == L || stamp[q] == gen) continue;  // stamp: also skips F (disjoint by now)
          if (ord[q] < lo) continue;
          stamp[q] = gen;
          Bset.push_back(q);
          stack.push_back(q);
        }
      }
      pairStatVisits += visits;
      pairStatEdgesOk += edgesHere;
      // new order in the region: B, the joined piece, F (each keeping its order)
      slots.clear();
      for (int q : Bset) slots.push_back(ord[q]);
      for (int q : Fset) slots.push_back(ord[q]);
      slots.push_back(ord[P]);
      slots.push_back(ord[Q]);
      std::sort(slots.begin(), slots.end());
      auto byOrd = [&](int x, int y) { return ord[x] < ord[y]; };
      std::sort(Bset.begin(), Bset.end(), byOrd);
      std::sort(Fset.begin(), Fset.end(), byOrd);
      size_t k = 0;
      for (int q : Bset) ord[q] = slots[k++];
      const int joinedOrd = slots[k++];
      k++;  // a spare slot
      for (int q : Fset) ord[q] = slots[k++];
      // contract
      int big = outs[P].size() + ins[P].size() >= outs[Q].size() + ins[Q].size() ? P : Q, small = big == P ? Q : P;
      parent[small] = big;
      outs[big].insert(outs[big].end(), outs[small].begin(), outs[small].end());
      ins[big].insert(ins[big].end(), ins[small].begin(), ins[small].end());
      std::vector<int>().swap(outs[small]);
      std::vector<int>().swap(ins[small]);
      compact(outs[big], outsKept[big], big);
      compact(ins[big], insKept[big], big);
      ord[big] = joinedOrd;
      phead[big] = phead[P];
      ptail[big] = ptail[Q];
      return 1;
    };
    using E = std::tuple<int, int64_t, int, int>;  // (priority, tie, i, j): highest priority, then tie
    auto tieI = [&](int i) -> int64_t { return pairSeed ? (int64_t)(mix(~(uint64_t)i) >> 1) : -(int64_t)i; };
    std::priority_queue<E> heap;
    for (int i = 0; i < n; i++) {
      auto [s, j] = bestFrom(i);
      if (j >= 0) heap.push({prio(s, i, j), tieI(i), i, j});
    }
    long steps = 0;
    const auto tStart = std::chrono::steady_clock::now();
    const long joinsMax = n;  // at most n - 1 joins; projection by joins done
    while (!heap.empty()) {
      if ((++steps & 1023) == 0) {
        if (stop && stop->load(std::memory_order_relaxed)) break;
        if (budget > 0 && pairStatJoins > 0) {
          double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count();
          double proj = el * joinsMax / pairStatJoins;
          if (proj > budget && (pairLim > 8 || pairVisitCap > 8)) { pairLim = std::max(8, pairLim / 2); pairVisitCap = std::max(8, pairVisitCap / 2); }
          else if (pairLim <= 8 && pairVisitCap <= 8 && el > 2 * budget) break;
        }
      }
      auto [pr, tie, i, j] = heap.top();
      heap.pop();
      if (succ[i] >= 0) continue;  // no longer a piece end
      const int P = find(i), Q = find(j);
      if (pred[j] >= 0 || P == Q) {  // stale target: look again
        auto [s2, j2] = bestFrom(i);
        if (j2 >= 0) heap.push({prio(s2, i, j2), tieI(i), i, j2});
        continue;
      }
      const int ok = tryJoin(P, Q);
      if (ok != 1) {
        pairStatRejected++;
        if (ok < 0) pairStatCapped++;
        rejected[i].push_back(j);
        auto [s2, j2] = bestFrom(i);
        if (j2 >= 0) heap.push({prio(s2, i, j2), tieI(i), i, j2});
        continue;
      }
      pairStatJoins++;
      succ[i] = j;
      pred[j] = i;
      const int ti = ptail[find(i)];  // the new piece end looks for its best seam
      auto [s2, j2] = bestFrom(ti);
      if (j2 >= 0) heap.push({prio(s2, ti, j2), tieI(ti), ti, j2});
    }
    // lay out the pieces in their topological order
    std::vector<int> roots;
    for (int i = 0; i < n; i++) if (find(i) == i) roots.push_back(i);
    std::sort(roots.begin(), roots.end(), [&](int x, int y) { return ord[x] < ord[y]; });
    std::vector<int> seq;
    for (int r : roots)
      for (int x = phead[r]; x >= 0; x = succ[x]) seq.push_back(x);
    if ((int)seq.size() != n) { long r, m; return merge(r, m, order); }  // cannot happen; be safe
    tail.clear();
    cost = 0;
    for (int i : seq)
      for (auto &p : insts[i].moves) push(tail, cost, p);
    mergedMoves = cost;
    if (order) *order = seq;
    std::vector<mv> out;
    for (auto &g : tail)
      for (auto &l : g.layers) out.push_back(toMove({g.axis, l.first, l.second}));
    return out;
  }


  // ---------------------------------------------------------------- best merge
  // Constraints as merge() builds them: pr[i] must precede i, su[i] follow it.
  void constraints(std::vector<std::vector<int>> &pr, std::vector<std::vector<int>> &su) const {
    const int n = insts.size();
    pr.assign(n, {});
    su.assign(n, {});
    std::vector<std::vector<std::pair<int, uint64_t>>> comp;
    for (int i = 0; i < n; i++)
      for (auto &d : insts[i].deps) {
        if ((int)comp.size() <= d.first) comp.resize(d.first + 1);
        for (auto &p : comp[d.first])
          if (p.second & d.second) { pr[i].push_back(p.first); su[p.first].push_back(i); }
        comp[d.first].push_back({i, d.second});
      }
  }
  // exact merged cost of an order
  long orderCost(const std::vector<int> &ord) const {
    std::vector<Group> t;
    long c = 0;
    for (int i : ord)
      for (auto &p : insts[i].moves) push(t, c, p);
    return c;
  }
  // Levels: each instance on the latest level its successors allow; each
  // level merged pairwise (no constraints inside a level); levels in order.
  std::vector<int> levelsOrder() const {
    const int n = insts.size();
    std::vector<std::vector<int>> pr, su;
    constraints(pr, su);
    std::vector<int> lev(n, 0), out;
    int L = 0;
    for (int i = n - 1; i >= 0; i--) {
      for (int s : su[i]) lev[i] = std::max(lev[i], lev[s] + 1);
      L = std::max(L, lev[i]);
    }
    std::vector<std::vector<int>> byLev(L + 1);
    for (int i = 0; i < n; i++) byLev[L - lev[i]].push_back(i);
    for (auto &ids : byLev) {
      Merger m(N);
      for (int i : ids) { MergeInst J; J.moves = insts[i].moves; m.insts.push_back(std::move(J)); }
      long r, x;
      std::vector<int> o;
      m.mergePairwise(r, x, &o);
      for (int k : o) out.push_back(ids[k]);
    }
    return out;
  }
  // Local improvement of an order: move a run of 1..maxRun consecutive
  // instances to just after a partner its first instance seams with, or just
  // before a partner its last instance seams with (partners from buckets of
  // first / last groups, candLim each way), when every predecessor stays
  // before the run and every successor after it (order labels) and the
  // pairwise seam savings rise.  Passes until nothing gains (at most maxPass).
  int impRun = 3, impPass = 20, impCand = 64;
  // improvement stops at this time (mergeBest sets it from `budget`)
  std::chrono::steady_clock::time_point impDeadline = std::chrono::steady_clock::time_point::max();
  std::vector<int> improveOrder(const std::vector<int> &start) const {
    const int n = insts.size();
    if (n < 2) return start;
    std::vector<std::vector<int>> pr, su;
    constraints(pr, su);
    std::vector<std::vector<int>> goff(n);
    for (int i = 0; i < n; i++) {
      const auto &m = insts[i].moves;
      for (size_t k = 0; k < m.size(); k++)
        if (k == 0 || m[k].axis != m[k - 1].axis) goff[i].push_back(k);
      goff[i].push_back(m.size());
    }
    auto save = [&](int i, int j) {
      if (i < 0 || j < 0) return 0;
      const auto &A = insts[i].moves, &B = insts[j].moves;
      const auto &ga = goff[i], &gb = goff[j];
      int s = 0;
      for (int gi = (int)ga.size() - 2, gj = 0; gi >= 0 && gj + 1 < (int)gb.size(); gi--, gj++) {
        const int a0 = ga[gi], a1 = ga[gi + 1], b0 = gb[gj], b1 = gb[gj + 1];
        if (A[a0].axis != B[b0].axis) break;
        int cancelled = 0;
        for (int h = b0; h < b1; h++)
          for (int t = a0; t < a1; t++)
            if (A[t].pos == B[h].pos) {
              if (((A[t].amt + B[h].amt) & 3) == 0) { s += 2; cancelled++; }
              else s += 1;
              break;
            }
        if (cancelled != a1 - a0 || cancelled != b1 - b0) break;
      }
      return s;
    };
    auto bkey = [&](int axis, int pos, int amt) { return ((size_t)axis * (N + 2) + pos) * 4 + amt; };
    std::vector<std::vector<int>> headB(3 * (size_t)(N + 2) * 4), tailB(3 * (size_t)(N + 2) * 4);
    for (int i = 0; i < n; i++) {
      const auto &m = insts[i].moves;
      if (m.empty()) continue;
      for (int k = goff[i][0]; k < goff[i][1]; k++) headB[bkey(m[k].axis, m[k].pos, m[k].amt)].push_back(i);
      const int g = goff[i].size() - 2;
      for (int k = goff[i][g]; k < goff[i][g + 1]; k++) tailB[bkey(m[k].axis, m[k].pos, m[k].amt)].push_back(i);
    }
    std::vector<int> nx(n, -1), pv(n, -1);
    std::vector<long double> lab(n);
    for (int k = 0; k < n; k++) {
      lab[start[k]] = k;
      if (k) { pv[start[k]] = start[k - 1]; nx[start[k - 1]] = start[k]; }
    }
    int first = start[0];
    for (int pass = 0; pass < impPass; pass++) {
      if (stop && stop->load(std::memory_order_relaxed)) break;
      long accepted = 0;
      bool late = false;
      for (int x0 = 0; x0 < n && !late; x0++) {
        if ((x0 & 255) == 0 && std::chrono::steady_clock::now() > impDeadline) late = true;
        int xe = x0;
        for (int len = 1; len <= impRun; len++) {
          if (len > 1) { xe = nx[xe]; if (xe < 0) break; }
          const int P = pv[x0], Nn = nx[xe];
          const long double r0 = lab[x0], r1 = lab[xe];
          auto inRun = [&](int y) { return y >= 0 && lab[y] >= r0 && lab[y] <= r1; };
          long double lo = -1e30L, hi = 1e30L;
          bool hasPred = false, hasSucc = false;
          for (int y = x0;; y = nx[y]) {
            for (int p : pr[y]) if (!inRun(p)) { lo = std::max(lo, lab[p]); hasPred = true; }
            for (int t : su[y]) if (!inRun(t)) { hi = std::min(hi, lab[t]); hasSucc = true; }
            if (y == xe) break;
          }
          const int removeGain = save(P, Nn) - save(P, x0) - save(xe, Nn);
          int bestGain = 0, bu = -2, bv = -2, evaluated = 0;
          auto consider = [&](int u, int v) {
            if (inRun(u) || inRun(v) || (u == P && v == Nn)) return;
            if (u >= 0 ? lab[u] < lo : hasPred) return;
            if (v >= 0 ? lab[v] > hi : hasSucc) return;
            const int g = removeGain + save(u, x0) + save(xe, v) - save(u, v);
            if (g > bestGain) { bestGain = g; bu = u; bv = v; }
          };
          {
            const auto &m = insts[x0].moves;
            for (int k = goff[x0][0]; k < goff[x0][1]; k++)
              for (int amt = 1; amt <= 3; amt++)
                for (int u : tailB[bkey(m[k].axis, m[k].pos, amt)]) {
                  if (++evaluated > impCand) break;
                  consider(u, nx[u] == x0 ? Nn : nx[u]);
                }
            evaluated = 0;
            const auto &me = insts[xe].moves;
            const int g = goff[xe].size() - 2;
            for (int k = goff[xe][g]; k < goff[xe][g + 1]; k++)
              for (int amt = 1; amt <= 3; amt++)
                for (int v : headB[bkey(me[k].axis, me[k].pos, amt)]) {
                  if (++evaluated > impCand) break;
                  consider(pv[v] == xe ? P : pv[v], v);
                }
          }
          if (bestGain <= 0) continue;
          if (P >= 0) nx[P] = Nn; else first = Nn;
          if (Nn >= 0) pv[Nn] = P;
          if (bu >= 0) nx[bu] = x0; else first = x0;
          pv[x0] = bu;
          nx[xe] = bv;
          if (bv >= 0) pv[bv] = xe;
          const long double Lb = bu >= 0 ? lab[bu] : (bv >= 0 ? lab[bv] - (len + 1) : 0), Rb = bv >= 0 ? lab[bv] : Lb + len + 1;
          if (Rb - Lb < 1e-6L) {  // relabel everything
            long double v = 0;
            for (int y = first; y >= 0; y = nx[y]) lab[y] = v++;
          } else {
            const long double step = (Rb - Lb) / (len + 1);
            long double v = Lb + step;
            for (int y = x0;; y = nx[y]) { lab[y] = v; v += step; if (y == xe) break; }
          }
          accepted++;
          break;
        }
      }
      if (!accepted || late) break;
    }
    std::vector<int> out;
    for (int x = first; x >= 0; x = nx[x]) out.push_back(x);
    return out;
  }
  // The merge: pairwise (levels above pairMaxInst instances), improved; then
  // seeded pairwise restarts, each improved, while the time `budget` allows
  // and at most `restarts` of them; the cheapest order wins.
  int restarts = 256;
  int threads = 1;  // restarts run on this many threads (the result does not depend on it)
  int pairMaxInst = 20000;
  long statRestarts = 0;
  std::vector<mv> mergeBest(long &rawMoves, long &mergedMoves, std::vector<int> *order = nullptr) {
    const auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    rawMoves = 0;
    for (auto &in : insts) rawMoves += in.moves.size();
    const int n = insts.size();
    impDeadline = budget > 0 ? t0 + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(budget))
                             : std::chrono::steady_clock::time_point::max();
    std::vector<int> best;
    if (n <= pairMaxInst) {
      long r, m;
      pairSeed = 0;
      mergePairwise(r, m, &best);
    } else
      best = levelsOrder();
    best = improveOrder(best);
    long bestCost = orderCost(best);
    statRestarts = 0;
    if (n <= pairMaxInst && restarts > 0) {
      // seeded restarts on `threads` threads; each improved; the cheapest order
      // wins, ties to the lower seed (so the result does not depend on timing)
      const double one = elapsed();  // a restart costs about as much as the first run
      std::atomic<int> nextSeed{1};
      std::atomic<long> done{0};
      std::mutex mu;
      int bestSeed = 0;
      auto worker = [&]() {
        Merger m = *this;
        m.budget = 0;
        m.stop = stop;
        for (;;) {
          if (stop && stop->load(std::memory_order_relaxed)) break;
          if (budget > 0 && elapsed() + one > budget) break;
          const int k = nextSeed++;
          if (k > restarts) break;
          m.pairSeed = k;
          long r, x;
          std::vector<int> o;
          m.mergePairwise(r, x, &o);
          o = m.improveOrder(o);
          const long c = m.orderCost(o);
          done++;
          std::lock_guard<std::mutex> g(mu);
          if (c < bestCost || (c == bestCost && k < bestSeed)) { bestCost = c; bestSeed = k; best = std::move(o); }
        }
      };
      const int T = std::max(1, std::min(threads, restarts));
      std::vector<std::thread> th;
      for (int t = 1; t < T; t++) th.emplace_back(worker);
      worker();
      for (auto &t : th) t.join();
      statRestarts = done;
    }
    pairSeed = 0;
    tail.clear();
    cost = 0;
    for (int i : best)
      for (auto &p : insts[i].moves) push(tail, cost, p);
    mergedMoves = cost;
    if (order) *order = best;
    std::vector<mv> out;
    for (auto &g : tail)
      for (auto &l : g.layers) out.push_back(toMove({g.axis, l.first, l.second}));
    return out;
  }

  std::vector<MergeInst> insts;
  int N;
  const std::atomic<bool> *stop = nullptr;  // checked every 4096 merge steps

 public:  // (Group and push are public for devtools experiments)
  struct Group {
    int8_t axis;
    std::vector<std::pair<int16_t, int8_t>> layers;  // (pos, amount), amount != 0
  };
  static void push(std::vector<Group> &t, long &c, const PLayer &p) {
    if (!t.empty() && t.back().axis == p.axis) {
      auto &L = t.back().layers;
      for (size_t k = 0; k < L.size(); k++)
        if (L[k].first == p.pos) {
          int na = (L[k].second + p.amt) & 3;
          if (na) L[k].second = na;
          else {
            L.erase(L.begin() + k);
            c--;
            if (L.empty()) t.pop_back();
          }
          return;
        }
      L.push_back({p.pos, p.amt});
      c++;
    } else {
      t.push_back({p.axis, {{p.pos, p.amt}}});
      c++;
    }
  }
 private:
  std::vector<Group> tail;
  long cost = 0;
  // cost increase of appending mv (only the tail can interact)
  int delta(const std::vector<PLayer> &mvs) const {
    size_t keep = std::min(tail.size(), mvs.size() + 1);
    std::vector<Group> t(tail.end() - keep, tail.end());
    long c = 0;
    for (auto &p : mvs) push(t, c, p);
    return (int)c;
  }
};
