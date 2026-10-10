// Table beam for one wing orbit over a pool of wing algorithms (centres free).
//
// Pool lines: "COST s0..s47  moves" (slots 24..47 unused); after algorithm a,
// wing slot j holds what was at src[j]; '2' in the moves is the orbit's depth.
// State: v[slot] = piece (24 distinct pieces), solved when v[j] == j.
//
// Wings in place after a: count of slots p with v[p] == dest_a(p), where dest_a
// is the inverse of src_a.  With the piece ids as 5 bit-planes V0..V4 over the
// 24 slots and dest_a as planes T0..T4, that is
// popcount(~((V0^T0)|..|(V4^T4)) & 24 bits): a few XORs and one popcount.
//
// Level-by-level beam with a direct-mapped next-level table (Zobrist hash over
// the algorithm's moved slots, seeded slot mapping and tie-break), a table of
// seen hashes, and a finish table (hash of every state one algorithm solves ->
// cheapest such algorithm) probed for every node and admitted child.  Score:
// mul * (wings in place) - 4 * (moves so far).
#pragma once
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <functional>
#include <algorithm>
#include <random>
#include <string>
#include <vector>
#include <thread>
#include <chrono>
#include "mintree.h"

struct WingPool {
  std::vector<std::array<uint8_t, 24>> src;
  std::vector<uint8_t> cost;
  std::vector<std::array<uint32_t, 5>> T;  // dest planes
  std::vector<uint32_t> suppMask;          // moved slots as a bit mask
  // same-axis groups (seq.h encoding, depths 1..2 per side) at the two ends,
  // for cancellation at the seam with the previous algorithm
  std::vector<uint16_t> firstG, lastG;
  std::vector<uint8_t> ngroups;
  // moves (depth 1 or 2, face ULFRBD, twist): stored when loaded from an
  // expanded pool file, else generated on demand from a reference (a symmetry
  // class and image) by movesFn
  using Moves = std::vector<std::array<int8_t, 3>>;
  std::vector<Moves> moves;
  std::vector<uint32_t> ref;
  std::function<Moves(uint32_t)> movesFn;
  Moves movesOf(size_t a) const { return moves.empty() ? movesFn(ref[a]) : moves[a]; }

  static bool parseMoves(const char *q, Moves &mv3) {  // "2F' U ..." (depth prefix 2 = the orbit's depth)
    while (*q) {
      while (*q == ' ') q++;
      if (!*q || *q == '\n') break;
      int d = 1;
      if (*q >= '0' && *q <= '9') d = *q++ - '0';
      const char *F = "ULFRBD";
      const char *fp = *q ? strchr(F, *q) : nullptr;
      if (!fp) return false;
      q++;
      int tw = 1;
      if (*q == '2') { tw = 2; q++; } else if (*q == '\'') { tw = 3; q++; }
      mv3.push_back({(int8_t)d, (int8_t)(fp - F), (int8_t)tw});
    }
    return true;
  }
  // expanded pool file: "COST s0..s47  moves"
  bool load(const char *fn, int maxCost = 99) {
    FILE *f = fopen(fn, "r");
    if (!f) return false;
    char buf[4096];
    while (fgets(buf, sizeof buf, f)) {
      char *p = buf;
      int c = (int)strtol(p, &p, 10);
      if (c > maxCost) continue;
      std::array<uint8_t, 24> s;
      for (int k = 0; k < 48; k++) {
        int x = (int)strtol(p, &p, 10);
        if (k < 24) s[k] = x;
      }
      Moves mv3;
      if (!parseMoves(p, mv3)) { fclose(f); return false; }
      add(c, s, mv3);
      moves.push_back(mv3);
    }
    fclose(f);
    return true;
  }
  // one algorithm (moves are only used for the seam groups here)
  void add(int c, const std::array<uint8_t, 24> &s, const Moves &mv3) {
    std::array<uint32_t, 5> planes{};
    uint32_t sm = 0;
    for (int j = 0; j < 24; j++) {
      int from = s[j];  // piece from slot `from` lands on j: dest(from) = j
      for (int b = 0; b < 5; b++)
        if (j >> b & 1) planes[b] |= 1u << from;
      if (from != j) sm |= 1u << j;
    }
    std::vector<uint16_t> g;
    static const int AX[6] = {0, 1, 2, 1, 2, 0}, SD[6] = {0, 0, 0, 1, 1, 1};
    for (auto &m : mv3) {
      int ax = AX[m[1]], slot = SD[m[1]] * 3 + (m[0] - 1);
      uint16_t amt = (uint16_t)(m[2] << (2 * slot));
      if (!g.empty() && (g.back() >> 12) == ax) {
        uint16_t lo = (g.back() ^ amt) & 0x555, carry = (g.back() & amt & 0x555) << 1;
        uint16_t sum = (lo | ((g.back() ^ amt ^ carry) & 0xAAA)) & 0xFFF;
        if (!sum) g.pop_back(); else g.back() = (ax << 12) | sum;
      } else g.push_back((ax << 12) | amt);
    }
    firstG.push_back(g.empty() ? 0xFFFF : g.front());
    lastG.push_back(g.empty() ? 0xFFFF : g.back());
    ngroups.push_back(g.size());
    suppMask.push_back(sm);
    src.push_back(s);
    cost.push_back(c);
    T.push_back(planes);
  }
  size_t size() const { return cost.size(); }
};

struct WingTableBeam {
  static constexpr uint32_t FIN = 1u << 31;  // result index flag: algorithm from finPool
  const WingPool *pool = nullptr;
  uint64_t Z[24][24];
  int mul = 8;  // quarter-moves per wing in place (2 moves; 4..10 measured, 8 best)
  uint64_t seed = 1;
  const std::atomic<bool> *stop = nullptr;
  int threads = 1;  // threads per beam level (set by the caller: spare threads beyond one per orbit)
  // finish table
  std::vector<uint64_t> fkey;
  std::vector<uint32_t> falg;
  uint64_t fmask = 0;
  bool useFinish = true;
  int keepSolved = 2;  // reject algorithms moving more than (keepSolved - 1) wings already in place (0: off)

  const WingPool *finPool = nullptr;  // finish-table algorithms (may be a larger pool than the beam's)
  // per slot: bitset over the pool's algorithms that move that slot (for keepSolved)
  size_t nw = 0;
  std::vector<uint64_t> BS;  // 24 * nw words
  void build(const WingPool &p, const WingPool *fin = nullptr) {
    pool = &p;
    finPool = fin ? fin : &p;
    const WingPool &FP = *finPool;
    std::mt19937_64 r(777);
    for (auto &z : Z) for (auto &x : z) x = r();
    size_t sz = 1;
    while (sz < 2 * FP.size()) sz <<= 1;
    fkey.assign(sz, 0);
    falg.assign(sz, 0);
    fmask = sz - 1;
    nw = (p.size() + 63) / 64;
    BS.assign(24 * nw, 0);
    for (size_t a = 0; a < p.size(); a++)
      for (uint32_t m = p.suppMask[a]; m; m &= m - 1) BS[__builtin_ctz(m) * nw + a / 64] |= 1ULL << (a % 64);
    for (size_t a = 0; a < FP.size(); a++) {
      uint64_t k = 0;  // the state a solves: slot src[j] holds piece j
      for (int j = 0; j < 24; j++) k ^= Z[FP.src[a][j]][j];
      for (size_t i = k & fmask;; i = (i + 1) & fmask) {
        if (fkey[i] == k) { if (FP.cost[a] < FP.cost[falg[i]]) falg[i] = a; break; }
        if (fkey[i] == 0) { fkey[i] = k; falg[i] = a; break; }
      }
    }
  }
  long findFinish(uint64_t k) const {
    for (size_t i = k & fmask;; i = (i + 1) & fmask) {
      if (fkey[i] == k) return falg[i];
      if (fkey[i] == 0) return -1;
    }
  }
  uint64_t zob(const std::array<uint8_t, 24> &v) const {
    uint64_t h = 0;
    for (int j = 0; j < 24; j++) h ^= Z[j][v[j]];
    return h;
  }

  // Solve one orbit; returns pool indices (empty if already solved; ok=false on failure).
  const WingPool &poolOf(uint32_t r) const { return r & FIN ? *finPool : *pool; }
  std::vector<uint32_t> solve(const std::array<uint8_t, 24> &start, int width, bool &ok, uint64_t runSeed = 0) const {
    const uint64_t seed = runSeed ? runSeed : this->seed;
    ok = false;
    const WingPool &P = *pool;
    struct Node { std::array<uint8_t, 24> v; uint64_t h; int g4, parent; uint32_t alg; uint16_t tail; };
    struct Slot { int score; int parent; uint32_t alg; int g4; uint64_t h; uint16_t tail; };
    // seam: moves saved when algorithm (first group f) follows a sequence whose last
    // group is t; nt = the new last group if the algorithm is a single group
    auto gnz = [](uint16_t amt) { return __builtin_popcount((amt | (amt >> 1)) & 0x555); };
    auto seam = [&](uint16_t t, uint16_t f, uint16_t &merged) {
      merged = 0xFFFF;
      if (t == 0xFFFF || f == 0xFFFF || (t >> 12) != (f >> 12)) return 0;
      uint16_t a = t & 0xFFF, b = f & 0xFFF;
      uint16_t lo = (a ^ b) & 0x555, carry = (a & b & 0x555) << 1;
      uint16_t sum = (lo | ((a ^ b ^ carry) & 0xAAA)) & 0xFFF;
      merged = sum ? (uint16_t)((t & 0xF000) | sum) : 0xFFFF;
      return gnz(sum) - gnz(a) - gnz(b);
    };
    const int TS = std::max(1, width);
    const uint64_t seedMix = seed * 0x9E3779B97F4A7C15ULL;
    auto slotOf = [&](uint64_t h) {
      uint64_t x = h ^ seedMix;
      x ^= x >> 31; x *= 0xBF58476D1CE4E5B9ULL; x ^= x >> 29;
      return (size_t)(((unsigned __int128)x * (uint64_t)TS) >> 64);
    };
    auto tieKey = [&](uint64_t h) { uint64_t x = (h + seedMix) * 0xD6E8FEB86659FD93ULL; return x ^ (x >> 32); };
    // every node's (parent, algorithm) for the solution walk, 8 bytes each (a deque:
    // no doubling copies); full nodes only for the current level (cur, whose ids
    // are the contiguous level[0] ..) and the one being built (nxt)
    struct Step { int parent; uint32_t alg; };
    std::deque<Step> hist;
    std::vector<Node> cur, nxt;
    hist.push_back({-1, 0});
    cur.push_back({start, zob(start), 0, -1, 0, 0xFFFF});
    std::vector<int> level = {0};
    // duplicate filter sized to the beam: a fixed 16k overflowed at width 2048
    // and cancelling algorithm pairs then cycled until the depth cap
    int SEEN = 1 << 14;
    while (SEEN < 64 * TS && SEEN < (1 << 22)) SEEN <<= 1;
    std::vector<uint64_t> seen(SEEN, 0);
    seen[cur[0].h & (SEEN - 1)] = cur[0].h;
    int best4 = 1 << 30, bestNode = -1;
    long bestFin = -1;  // finishing algorithm after bestNode (or -1)
    long bestFinAlg = -1;  // algorithm between bestNode and the finish (-1: none)
    const size_t n = P.size();
    constexpr uint32_t M24 = (1u << 24) - 1;
    // Per-thread state for one level: its own candidate table (merged slot by
    // slot afterwards: each slot keeps its best candidate by score, then the
    // seeded tie-break, so the table hardly depends on the split) and its own
    // best solution.
    struct Ctx {
      std::vector<Slot> table;
      MinTree minTree;
      int filled = 0, weakest = INT32_MIN, best4 = 1 << 30, bestNode = -1;
      long bestFin = -1, bestFinAlg = -1;
      std::vector<uint64_t> ones, twos;
    };
    const int T = std::max(1, std::min(threads, 64));
    std::vector<Ctx> ctx(T);
    for (auto &c : ctx) c.table.resize(TS);
    auto runNodes = [&](Ctx &cx, size_t from, size_t to) {
      for (auto &sl : cx.table) sl.score = INT32_MIN;
      cx.minTree.reset(TS);
      cx.filled = 0;
      cx.weakest = INT32_MIN;
      cx.best4 = best4;
      cx.bestNode = -1;
      for (size_t li = from; li < to; li++) {
        if ((li & 63) == 0 && stop && stop->load(std::memory_order_relaxed)) return;  // Ctrl-C: within a level too
        const int ni = level[li];
        const Node &nd = cur[ni - level[0]];
        uint32_t V[5] = {0, 0, 0, 0, 0};
        int c0 = 0;
        for (int p = 0; p < 24; p++) {
          for (int b = 0; b < 5; b++) V[b] |= (uint32_t)(nd.v[p] >> b & 1) << p;
          c0 += nd.v[p] == p;
        }
        if (c0 == 24) {
          if (nd.g4 < cx.best4) { cx.best4 = nd.g4; cx.bestNode = ni; cx.bestFin = cx.bestFinAlg = -1; }
          continue;
        }
        if (useFinish) {
          long f = findFinish(nd.h);
          uint16_t mg;
          int fc4 = f >= 0 ? 4 * (finPool->cost[f] + seam(nd.tail, finPool->firstG[f], mg)) : 0;
          if (f >= 0 && nd.g4 + fc4 < cx.best4) { cx.best4 = nd.g4 + fc4; cx.bestNode = ni; cx.bestFin = f; cx.bestFinAlg = -1; }
        }
        // rejected: algorithms moving >= keepSolved wings already in place
        // (bitset "at least 1" / "at least 2" over the correct slots' bitsets)
        const uint64_t *rejp = nullptr;
        if (keepSolved == 1 || keepSolved == 2) {
          cx.ones.assign(nw, 0);
          cx.twos.assign(nw, 0);
          if (n % 64) cx.twos[nw - 1] = cx.ones[nw - 1] = ~0ULL << (n % 64);  // padding never survives
          for (int p = 0; p < 24; p++)
            if (nd.v[p] == p) {
              const uint64_t *b = &BS[p * nw];
              for (size_t w = 0; w < nw; w++) { cx.twos[w] |= cx.ones[w] & b[w]; cx.ones[w] |= b[w]; }
            }
          rejp = keepSolved == 1 ? cx.ones.data() : cx.twos.data();
        }
        auto score = [&](size_t a) {
          const auto &Tp = P.T[a];
          int c = __builtin_popcount(~((V[0] ^ Tp[0]) | (V[1] ^ Tp[1]) | (V[2] ^ Tp[2]) | (V[3] ^ Tp[3]) | (V[4] ^ Tp[4])) & M24);
          uint16_t merged;
          const int sv = seam(nd.tail, P.firstG[a], merged);
          if (P.cost[a] + sv <= 0) return;  // adds no moves: only undoes (and lets the beam cycle)
          const int g4 = nd.g4 + 4 * (P.cost[a] + sv);
          const uint16_t ntail = P.ngroups[a] == 1 ? merged == 0xFFFF && sv == 0 ? P.lastG[a] : merged : P.lastG[a];
          int sc = mul * c - g4;
          if (sc < cx.weakest || g4 >= cx.best4) return;
          uint64_t h = nd.h;
          const auto &src = P.src[a];
          for (uint32_t m = P.suppMask[a]; m; m &= m - 1) {
            const int j = __builtin_ctz(m);
            h ^= Z[j][nd.v[j]] ^ Z[j][nd.v[src[j]]];
          }
          if (seen[h & (SEEN - 1)] == h) return;
          if (useFinish) {
            long f = findFinish(h);
            uint16_t mg;
            int fc4 = f >= 0 ? 4 * (finPool->cost[f] + seam(ntail, finPool->firstG[f], mg)) : 0;
            if (f >= 0 && g4 + fc4 < cx.best4) { cx.best4 = g4 + fc4; cx.bestNode = ni; cx.bestFin = f; cx.bestFinAlg = a; }
          }
          const size_t si = slotOf(h);
          Slot &sl = cx.table[si];
          if (sl.score == INT32_MIN) cx.filled++;
          else if (sc < sl.score || (sc == sl.score && tieKey(h) >= tieKey(sl.h))) return;
          sl = {sc, ni, (uint32_t)a, g4, h, ntail};
          cx.minTree.set(si, sc);
          if (cx.filled == TS) cx.weakest = cx.minTree.min();
        };
        if (rejp) {
          for (size_t w = 0; w < nw; w++)
            for (uint64_t bits = ~rejp[w]; bits; bits &= bits - 1) score(w * 64 + __builtin_ctzll(bits));
        } else
          for (size_t a = 0; a < n; a++) score(a);
      }
    };
    const auto tBeam0 = std::chrono::steady_clock::now();
    for (int depth = 0; depth < 200 && !level.empty(); depth++) {
      if (getenv("WINGSTAT"))
        fprintf(stderr, "    depth %d: level %zu, %.2fs, best %d\n", depth, level.size(),
                std::chrono::duration<double>(std::chrono::steady_clock::now() - tBeam0).count(), best4);
      if (stop && stop->load(std::memory_order_relaxed)) return {};
      // split the level over the threads (one thread for small levels)
      const int U = level.size() >= (size_t)(8 * T) ? T : 1;
      if (U == 1) runNodes(ctx[0], 0, level.size());
      else {
        std::vector<std::thread> th;
        for (int u = 1; u < U; u++) th.emplace_back(runNodes, std::ref(ctx[u]), level.size() * u / U, level.size() * (u + 1) / U);
        runNodes(ctx[0], 0, level.size() / U);
        for (auto &t : th) t.join();
      }
      if (stop && stop->load(std::memory_order_relaxed)) return {};
      for (int u = 0; u < U; u++)
        if (ctx[u].bestNode >= 0 && ctx[u].best4 < best4) {
          best4 = ctx[u].best4; bestNode = ctx[u].bestNode; bestFin = ctx[u].bestFin; bestFinAlg = ctx[u].bestFinAlg;
        }
      std::vector<Slot> &table = ctx[0].table;
      for (int u = 1; u < U; u++)
        for (int i = 0; i < TS; i++) {
          const Slot &o = ctx[u].table[i];
          Slot &sl = table[i];
          if (o.score == INT32_MIN) continue;
          if (sl.score == INT32_MIN || o.score > sl.score || (o.score == sl.score && tieKey(o.h) < tieKey(sl.h))) sl = o;
        }
      const int curBase = level[0];
      level.clear();
      nxt.clear();
      for (auto &sl : ctx[0].table) {
        if (sl.score == INT32_MIN) continue;
        Node c;
        const Node &pn = cur[sl.parent - curBase];
        for (int j = 0; j < 24; j++) c.v[j] = pn.v[P.src[sl.alg][j]];
        c.h = sl.h; c.g4 = sl.g4; c.parent = sl.parent; c.alg = sl.alg; c.tail = sl.tail;
        seen[c.h & (SEEN - 1)] = c.h;
        hist.push_back({sl.parent, sl.alg});
        nxt.push_back(c);
        level.push_back(hist.size() - 1);
      }
      cur.swap(nxt);
    }
    if (getenv("WINGSTAT")) {
      long scored = 0;
      fprintf(stderr, "  wing beam width %d: %zu nodes, best %d quarter-moves\n", TS, hist.size(), best4);
      (void)scored;
    }
    if (bestNode < 0) return {};
    std::vector<uint32_t> r;
    for (int x = bestNode; hist[x].parent >= 0; x = hist[x].parent) r.push_back(hist[x].alg);
    std::reverse(r.begin(), r.end());
    if (bestFinAlg >= 0) r.push_back(bestFinAlg);
    if (bestFin >= 0) r.push_back((uint32_t)bestFin | FIN);
    ok = true;
    return r;
  }
};

