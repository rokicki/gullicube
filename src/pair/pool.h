// Pool of pair algorithms in the layout the search wants.
//
// State: col[48], colour held at each position; target colour of p is (p%24)/4.
// An algorithm with source list s (after it, position d holds what was at s_d)
// is stored as 6 source masks  S_c = { s_d : target(d) == c }.  Then the number
// of correct positions after applying it to a state with colour masks M_c is
//     sum_c popcount(S_c & M_c)
// -- 6 AND + 6 POPCNT, no gathers.
#pragma once
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <random>
#include "seq.h"

inline int target(int p) { return (p % 24) / 4; }

struct Pool {
  std::vector<std::array<uint64_t, 6>> S;
  std::vector<std::array<uint8_t, 48>> src;
  std::vector<uint8_t> cost;
  std::vector<uint64_t> support;
  std::vector<std::string> text;
  // 1 if the algorithm permutes pieces of the diagonal orbits (i,i)/(j,j):
  // after the off-diagonal solve, the diagonal state must be updated by
  // exactly these algorithms before the diagonals are solved.  Read from a
  // trailing " D0"/" D1" column (diagflag -w); 1 (unknown) if absent.
  std::vector<uint8_t> diag;
  // merged group form of each algorithm: groups goff[a]..goff[a+1] of gdat
  std::vector<uint32_t> goff{0};
  std::vector<uint16_t> gdat;
  std::vector<uint16_t> head, second;  // first/second group key (0xFFFF if none)
  size_t size() const { return cost.size(); }
  const uint16_t *gbegin(size_t a) const { return gdat.data() + goff[a]; }
  const uint16_t *gend(size_t a) const { return gdat.data() + goff[a + 1]; }

  // append one algorithm: cost, src permutation, moves (text kept if keepText)
  void addEntry(int len, const std::array<uint8_t, 48> &s, const std::string &mv, uint8_t dflag, bool keepText = true) {
    std::array<uint64_t, 6> m{};
    uint64_t sup = 0;
    for (int d = 0; d < 48; d++) {
      m[target(d)] |= 1ULL << s[d];
      if (s[d] != d) sup |= 1ULL << d;
    }
    auto g = to_groups(parse_alg(mv.c_str()));
    head.push_back(g.size() > 0 ? g[0] : 0xFFFF);
    second.push_back(g.size() > 1 ? g[1] : 0xFFFF);
    gdat.insert(gdat.end(), g.begin(), g.end());
    goff.push_back(gdat.size());
    S.push_back(m);
    src.push_back(s);
    cost.push_back(len);
    support.push_back(sup);
    text.push_back(keepText ? mv : std::string());
    diag.push_back(dflag);
  }
  // append one algorithm given as moves (no text kept)
  void addEntry(int len, const std::array<uint8_t, 48> &s, const std::vector<Move> &mvs, uint8_t dflag) {
    std::array<uint64_t, 6> m{};
    uint64_t sup = 0;
    for (int d = 0; d < 48; d++) {
      m[target(d)] |= 1ULL << s[d];
      if (s[d] != d) sup |= 1ULL << d;
    }
    auto g = to_groups(mvs);
    head.push_back(g.size() > 0 ? g[0] : 0xFFFF);
    second.push_back(g.size() > 1 ? g[1] : 0xFFFF);
    gdat.insert(gdat.end(), g.begin(), g.end());
    goff.push_back(gdat.size());
    S.push_back(m);
    src.push_back(s);
    cost.push_back(len);
    support.push_back(sup);
    text.push_back(std::string());
    diag.push_back(dflag);
  }
  void load(const char *fn, int maxlen = 99, bool keepText = true) {
    std::ifstream in(fn);
    if (!in) { fprintf(stderr, "can't open %s\n", fn); exit(1); }
    std::string line;
    while (std::getline(in, line)) {
      const char *p = line.c_str();
      char *e;
      int len = strtol(p, &e, 10);
      if (e == p) continue;
      if (len > maxlen) continue;
      p = e;
      std::array<uint8_t, 48> s;
      for (int d = 0; d < 48; d++) {
        s[d] = strtol(p, &e, 10);
        p = e;
      }
      while (*p == ' ') p++;
      std::string mv(p);
      uint8_t dflag = 1;
      if (mv.size() >= 3 && mv[mv.size() - 3] == ' ' && mv[mv.size() - 2] == 'D' &&
          (mv.back() == '0' || mv.back() == '1')) {
        dflag = mv.back() - '0';
        mv.resize(mv.size() - 3);
      }
      addEntry(len, s, mv, dflag, keepText);
    }
  }
  // sort by cost (stable), so a scan can stop once cost alone exceeds a bound
  void sortByCost() {
    std::vector<size_t> idx(size());
    for (size_t i = 0; i < idx.size(); i++) idx[i] = i;
    std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return cost[a] < cost[b]; });
    Pool q;
    for (auto i : idx) {
      q.S.push_back(S[i]);
      q.src.push_back(src[i]);
      q.cost.push_back(cost[i]);
      q.support.push_back(support[i]);
      q.text.push_back(std::move(text[i]));
      q.diag.push_back(diag[i]);
      q.head.push_back(head[i]);
      q.second.push_back(second[i]);
      q.gdat.insert(q.gdat.end(), gbegin(i), gend(i));
      q.goff.push_back(q.gdat.size());
    }
    *this = std::move(q);
  }
};

struct State {
  std::array<uint8_t, 48> col;
  std::array<uint64_t, 6> M;
  void setMasks() {
    M = {};
    for (int p = 0; p < 48; p++) M[col[p]] |= 1ULL << p;
  }
  int correct() const {
    int c = 0;
    for (int p = 0; p < 48; p++) c += col[p] == target(p);
    return c;
  }
  State apply(const Pool &pool, size_t a) const {
    State r;
    for (int d = 0; d < 48; d++) r.col[d] = col[pool.src[a][d]];
    r.setMasks();
    return r;
  }
  uint64_t hash() const {
    uint64_t h = 0;
    for (int c = 0; c < 6; c++) h = (h ^ M[c]) * 0x9E3779B97F4A7C15ULL + c;
    return h ^ (h >> 29);
  }
};

inline int correctAfter(const std::array<uint64_t, 6> &S, const std::array<uint64_t, 6> &M) {
  return __builtin_popcountll(S[0] & M[0]) + __builtin_popcountll(S[1] & M[1]) +
         __builtin_popcountll(S[2] & M[2]) + __builtin_popcountll(S[3] & M[3]) +
         __builtin_popcountll(S[4] & M[4]) + __builtin_popcountll(S[5] & M[5]);
}

inline State randomState(std::mt19937_64 &rng) {
  State s;
  for (int o = 0; o < 2; o++) {
    std::array<uint8_t, 24> a;
    for (int i = 0; i < 24; i++) a[i] = i / 4;
    std::shuffle(a.begin(), a.end(), rng);
    for (int i = 0; i < 24; i++) s.col[o * 24 + i] = a[i];
  }
  s.setMasks();
  return s;
}
