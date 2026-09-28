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

  std::vector<MergeInst> insts;
  const std::atomic<bool> *stop = nullptr;  // checked every 4096 merge steps

 private:
  struct Group {
    int8_t axis;
    std::vector<std::pair<int16_t, int8_t>> layers;  // (pos, amount), amount != 0
  };
  int N;
  std::vector<Group> tail;
  long cost = 0;

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
  // cost increase of appending mv (only the tail can interact)
  int delta(const std::vector<PLayer> &mvs) const {
    size_t keep = std::min(tail.size(), mvs.size() + 1);
    std::vector<Group> t(tail.end() - keep, tail.end());
    long c = 0;
    for (auto &p : mvs) push(t, c, p);
    return (int)c;
  }
};
