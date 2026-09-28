// Endgame by enumeration + one-algorithm finish table.
//
// Finish table: Zobrist hash of every colour state some single algorithm solves
// -> cheapest such algorithm.  Search: enumerate up to `depth` algorithms whose
// support lies in (wrong set + at most one correct position), each making net
// progress, and look every reached state up in the table.  Branch and bound
// on total cost.
#pragma once
#include "pool.h"
#include <unordered_map>

struct EG2 {
  Pool cand;  // small-support algorithms used for enumeration
  std::unordered_map<uint64_t, std::vector<uint32_t>> bySupport;
  int maxSup = 0;
  // finish table (built from its own pool, may be larger)
  Pool fin;
  uint64_t Z[48][6];
  std::vector<uint64_t> keys;
  std::vector<uint32_t> vals;
  uint64_t mask = 0;
  int finMaxMis = 0;

  uint64_t zob(const std::array<uint8_t, 48> &c) const {
    uint64_t h = 0;
    for (int d = 0; d < 48; d++) h ^= Z[d][c[d]];
    return h;
  }
  void build() {
    cand.sortByCost();
    for (size_t a = 0; a < cand.size(); a++) {
      maxSup = std::max(maxSup, __builtin_popcountll(cand.support[a]));
      bySupport[cand.support[a]].push_back(a);
    }
    std::mt19937_64 r(99);
    for (auto &z : Z) for (auto &x : z) x = r();
    size_t sz = 1;
    while (sz < 2 * fin.size()) sz <<= 1;
    keys.assign(sz, 0);
    vals.assign(sz, 0);
    mask = sz - 1;
    size_t used = 0;
    for (size_t b = 0; b < fin.size(); b++) {
      std::array<uint8_t, 48> t;
      for (int d = 0; d < 48; d++) t[fin.src[b][d]] = target(d);
      int mis = 0;
      for (int d = 0; d < 48; d++) mis += t[d] != target(d);
      finMaxMis = std::max(finMaxMis, mis);
      uint64_t k = zob(t);
      for (size_t i = k & mask;; i = (i + 1) & mask) {
        if (keys[i] == k) { if (fin.cost[b] < fin.cost[vals[i]]) vals[i] = b; break; }
        if (keys[i] == 0) { keys[i] = k; vals[i] = b; used++; break; }
      }
    }
    fprintf(stderr, "EG2: %zu enumeration algs (%zu supports), finish table %zu algs -> %zu states\n", cand.size(),
            bySupport.size(), fin.size(), used);
  }
  int find(uint64_t k) const {
    for (size_t i = k & mask;; i = (i + 1) & mask) {
      if (keys[i] == k) return vals[i];
      if (keys[i] == 0) return -1;
    }
  }
  void collect(uint64_t W, uint64_t cur, int from, int size, std::vector<uint32_t> &out, int budget) const {
    if (size >= 2) {
      auto it = bySupport.find(cur);
      if (it != bySupport.end())
        for (auto a : it->second) {
          if (cand.cost[a] > budget) break;
          out.push_back(a);
        }
    }
    if (size == maxSup) return;
    for (int p = from; p < 48; p++)
      if (W >> p & 1) collect(W, cur | 1ULL << p, p + 1, size + 1, out, budget);
  }

  struct Sol {
    int cost = 1 << 30;
    std::vector<uint32_t> path;  // indices into cand
    int finAlg = -1;             // index into fin
  };

  int topK = 1 << 30;  // children recursed into per node (all get a table lookup)
  double rate = 2.4;   // moves per wrong piece, for ranking children

  // Node s (already looked up by the caller).  Enumerate progress-making
  // candidates, look each child up in the finish table, and recurse into the
  // best topK children while depth remains.
  void rec(const State &s, uint64_t h, int g, int depth, std::vector<uint32_t> &path, Sol &best, long &work) const {
    if (depth == 0 || g + 8 >= best.cost) return;
    uint64_t W = 0;
    for (int p = 0; p < 48; p++)
      if (s.col[p] != target(p)) W |= 1ULL << p;
    int m = __builtin_popcountll(W);
    std::vector<uint32_t> cs;
    int budget = best.cost - g - 1;
    collect(W, 0, 0, 0, cs, budget);
    for (int e = 0; e < 48; e++)
      if (!(W >> e & 1)) collect(W, 1ULL << e, 0, 1, cs, budget);
    int c0 = 48 - m;
    struct Ch { double key; uint32_t a; uint64_t h; int gain; };
    std::vector<Ch> kids;
    for (auto a : cs) {
      work++;
      int ng = g + cand.cost[a];
      if (ng >= best.cost) continue;
      int gain = correctAfter(cand.S[a], s.M) - c0;
      if (gain <= 0) continue;
      uint64_t h2 = h;
      const auto &sr = cand.src[a];
      for (int d = 0; d < 48; d++)
        if (sr[d] != d) h2 ^= Z[d][s.col[d]] ^ Z[d][s.col[sr[d]]];
      if (gain == m) {
        best.cost = ng; best.path = path; best.path.push_back(a); best.finAlg = -1;
        continue;
      }
      if (m - gain <= finMaxMis && ng + 4 < best.cost) {
        int b = find(h2);
        if (b >= 0 && ng + fin.cost[b] < best.cost) {
          best.cost = ng + fin.cost[b]; best.path = path; best.path.push_back(a); best.finAlg = b;
        }
      }
      if (depth > 1) kids.push_back({ng + rate * (m - gain), a, h2, gain});
    }
    if (depth <= 1) return;
    std::sort(kids.begin(), kids.end(), [](const Ch &x, const Ch &y) { return x.key < y.key; });
    if ((int)kids.size() > topK) kids.resize(topK);
    for (auto &k : kids) {
      int ng = g + cand.cost[k.a];
      if (ng + 8 >= best.cost) continue;
      path.push_back(k.a);
      rec(s.apply(cand, k.a), k.h, ng, depth - 1, path, best, work);
      path.pop_back();
    }
  }
  Sol solve(const State &s, int depth, int upper, long *workOut = nullptr) const {
    Sol best;
    best.cost = upper;
    std::vector<uint32_t> path;
    long work = 0;
    uint64_t h = zob(s.col);
    if (s.correct() == 48) { best.cost = 0; return best; }
    if (48 - s.correct() <= finMaxMis) {
      int b = find(h);
      if (b >= 0 && fin.cost[b] < best.cost) { best.cost = fin.cost[b]; best.finAlg = b; }
    }
    rec(s, h, 0, depth, path, best, work);
    if (workOut) *workOut = work;
    return best;
  }
};
