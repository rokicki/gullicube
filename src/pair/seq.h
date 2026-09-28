// Move sequences as stacks of same-axis groups, so that concatenating
// algorithms can merge and cancel moves at the seams.
//
// Layers on one axis commute, so a maximal run of same-axis moves is just an
// amount (0..3) per layer.  With depths 1..3 from each of the two faces that is
// 6 slots x 2 bits = 12 bits; a group key is axis<<12 | amounts.  The slice-turn
// cost of a group is its number of nonzero slots.
#pragma once
#include "cube.h"
#include <cstdint>
#include <vector>
#include <string>

static const int AXIS[6] = {0, 1, 2, 1, 2, 0};  // U L F R B D
static const int SIDE[6] = {0, 0, 0, 1, 1, 1};

inline int gnz(uint16_t amt) {
  uint16_t x = (amt | (amt >> 1)) & 0x555;
  return __builtin_popcount(x);
}
inline uint16_t gadd(uint16_t a, uint16_t b) {
  uint16_t lo = (a ^ b) & 0x555;
  uint16_t carry = (a & b & 0x555) << 1;
  uint16_t hi = (a ^ b ^ carry) & 0xAAA;
  return lo | hi;
}
inline int gaxis(uint16_t key) { return key >> 12; }
inline uint16_t gamt(uint16_t key) { return key & 0xFFF; }

inline std::vector<uint16_t> to_groups(const std::vector<Move> &alg) {
  std::vector<uint16_t> g;
  for (auto &m : alg) {
    int ax = AXIS[m.face];
    int slot = SIDE[m.face] * 3 + (m.depth - 1);
    uint16_t amt = (uint16_t)(m.amt << (2 * slot));
    if (!g.empty() && gaxis(g.back()) == ax) {
      uint16_t s = gadd(gamt(g.back()), amt);
      if (s == 0) g.pop_back();
      else g.back() = (ax << 12) | s;
    } else g.push_back((ax << 12) | amt);
  }
  return g;
}

// append groups b onto stack a, merging at the seam(s); returns new cost
inline void merge_into(std::vector<uint16_t> &a, const std::vector<uint16_t> &b) {
  for (auto k : b) {
    if (!a.empty() && gaxis(a.back()) == gaxis(k)) {
      uint16_t s = gadd(gamt(a.back()), gamt(k));
      if (s == 0) a.pop_back();
      else a.back() = (gaxis(k) << 12) | s;
    } else a.push_back(k);
  }
}
inline int groups_cost(const std::vector<uint16_t> &a) {
  int c = 0;
  for (auto k : a) c += gnz(gamt(k));
  return c;
}
inline std::string groups_str(const std::vector<uint16_t> &a) {
  static const int FACEOF[3][2] = {{0, 5}, {1, 3}, {2, 4}};
  std::string s;
  for (auto k : a)
    for (int slot = 0; slot < 6; slot++) {
      int amt = (gamt(k) >> (2 * slot)) & 3;
      if (!amt) continue;
      Move m{(uint8_t)FACEOF[gaxis(k)][slot / 3], (uint8_t)(slot % 3 + 1), (uint8_t)amt};
      if (!s.empty()) s += ' ';
      s += move_str(m);
    }
  return s;
}
