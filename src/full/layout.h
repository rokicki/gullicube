// Bridge between the pair-solver layout (pair/cube.h: 48 positions, orbit A
// then B, target colour (p%24)/4) and brobdicube's xcube, where centre orbit
// (x,y) lives in cube.rows[y-1][x-1].a[4*face + quadrant] and involves slice
// depths x+1 and y+1.
//
// A pair algorithm written with slices "2"/"3" acts on xorbits (1,2) and (2,1);
// substituting 2 -> a+1, 3 -> b+1 makes it act on xorbits (a,b) and (b,a) with
// the same index layout.  phi maps pair-solver positions to brobdicube
// positions: 0..23 = xorbit(1,2), 24..47 = xorbit(2,1).
#pragma once
#include "xcube.h"
#include <array>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace ps {  // pair-solver geometry (from pair/cube.h, trimmed)
static const int NORMAL[6][3] = {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}, {1, 0, 0}, {0, 0, -1}, {0, -1, 0}};
// quarter-turn maps on the 48 pair positions for faces x depths 1..3 (n = 8)
struct PairGeom {
  int dst[6][4][48];  // dst[f][d][p] = new position of the piece at p
  PairGeom() {
    const int n = 8;
    struct St { int p[3], face, row, col; };
    std::vector<St> st;
    for (int f = 0; f < 6; f++)
      for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++) {
          St s;
          int a = -(n - 1) + 2 * r, b = -(n - 1) + 2 * c;
          int ax = NORMAL[f][0] ? 0 : (NORMAL[f][1] ? 1 : 2);
          s.p[ax] = NORMAL[f][ax] * n;
          s.p[(ax + 1) % 3] = a;
          s.p[(ax + 2) % 3] = b;
          s.face = f; s.row = r; s.col = c;
          st.push_back(s);
        }
    auto dep = [&](int v) { return std::min(v, n - 1 - v) + 1; };
    auto find = [&](const int *q) {
      for (size_t i = 0; i < st.size(); i++)
        if (st[i].p[0] == q[0] && st[i].p[1] == q[1] && st[i].p[2] == q[2]) return (int)i;
      return -1;
    };
    // quarter-turn map on all stickers
    auto turn = [&](int f, int d, int i) {
      const int *N = NORMAL[f];
      int ax = N[0] ? 0 : (N[1] ? 1 : 2), sgn = N[ax];
      int layer = sgn * (n - 1 - 2 * (d - 1));
      int q[3] = {st[i].p[0], st[i].p[1], st[i].p[2]};
      int cub = q[ax];
      if (std::abs(q[ax]) == n) cub = q[ax] - (q[ax] > 0 ? 1 : -1);
      if (cub != layer) return i;
      int o1 = (ax + 1) % 3, o2 = (ax + 2) % 3, r[3] = {q[0], q[1], q[2]};
      if (sgn > 0) { r[o1] = q[o2]; r[o2] = -q[o1]; }
      else { r[o1] = -q[o2]; r[o2] = q[o1]; }
      return find(r);
    };
    // orbits of the (2,3) pair: union-find over all turns
    int S = st.size();
    std::vector<int> uf(S);
    for (int i = 0; i < S; i++) uf[i] = i;
    std::function<int(int)> fd = [&](int x) { return uf[x] == x ? x : uf[x] = fd(uf[x]); };
    for (int f = 0; f < 6; f++)
      for (int d = 1; d <= n; d++)
        for (int i = 0; i < S; i++) uf[fd(i)] = fd(turn(f, d, i));
    int orbA = -1, orbB = -1;
    for (int i = 0; i < S; i++)
      if (st[i].face == 0 && st[i].row == 1 && st[i].col == 2) orbA = fd(i);
    for (int i = 0; i < S; i++) {
      int dr = dep(st[i].row), dc = dep(st[i].col);
      if (((dr == 2 && dc == 3) || (dr == 3 && dc == 2)) && fd(i) != orbA) orbB = fd(i);
    }
    std::vector<int> idx(S, -1), pos;
    for (int w = 0; w < 2; w++)
      for (int i = 0; i < S; i++)
        if (fd(i) == (w ? orbB : orbA)) { idx[i] = pos.size(); pos.push_back(i); }
    for (int f = 0; f < 6; f++)
      for (int d = 1; d <= 3; d++)
        for (int p = 0; p < 48; p++) dst[f][d][p] = idx[turn(f, d, pos[p])];
  }
};
}  // namespace ps

// apply a move list to an xcube, substituting slice depths via depthOf(2 or 3)
inline void applyMapped(xcube &c, const std::vector<mv> &alg, const std::function<int(int)> &depthOf) {
  for (auto m : alg) {
    if (m.dep > 1) m.dep = depthOf(m.dep);
    domove(m, c, true);
  }
}

// Read centre orbit (x,y) (1-based) in physical quadrant order, the convention
// of xcube::get(): quadrant oz of face f is a[4f + ((oz + facerots[f]) & 3)].
// (xcube::normalizefacerots() rotates the other way, so we never call it.)
inline void readOrbit(const xcube &c, int x, int y, uint8_t out[24]) {
  const xorbit &o = c.rows[y - 1][x - 1];
  for (int f = 0; f < 6; f++)
    for (int oz = 0; oz < 4; oz++) out[4 * f + oz] = o.a[4 * f + ((oz + c.facerots[f]) & 3)];
}

struct Layout {
  int phi[48];  // pair-solver position -> brobdicube position (0..23 xorbit(1,2), 24..47 xorbit(2,1))
  // brobdicube quarter-turn maps on the same 48 positions
  static inline int faceTw = 1, sliceTw = 1;  // brobdicube twist matching our clockwise quarter turn
  static void brobDst(int f, int d, int out[48]) {
    xcube c(3, 3, 'U');
    domove({d, (short)f, (short)(d == 1 ? faceTw : sliceTw)}, c, true);
    uint8_t X[2][24];
    readOrbit(c, 1, 2, X[0]);
    readOrbit(c, 2, 1, X[1]);
    for (int w = 0; w < 2; w++)
      for (int k = 0; k < 24; k++) out[w * 24 + X[w][k]] = w * 24 + k;  // piece X[k] now at k
  }
  bool build() {
    ps::PairGeom G;
    static int B[6][4][48];
    for (int f = 0; f < 6; f++)
      for (int d = 1; d <= 3; d++) brobDst(f, d, B[f][d]);
    for (int w = 0; w < 2; w++) {  // pair-solver orbit A (w=0) and B (w=1)
      bool ok = false;
      for (int cand = 0; cand < 48 && !ok; cand++) {
        int base = w * 24;
        if (cand / 4 % 6 != 0) continue;  // base position 0/24 is on U
        int m[48];
        for (auto &x : m) x = -1;
        m[base] = cand;
        std::vector<int> q = {base};
        bool bad = false;
        while (!q.empty() && !bad) {
          int p = q.back();
          q.pop_back();
          for (int f = 0; f < 6 && !bad; f++)
            for (int d = 1; d <= 3 && !bad; d++) {
              int p2 = G.dst[f][d][p], b2 = B[f][d][m[p]];
              if (m[p2] < 0) { m[p2] = b2; q.push_back(p2); }
              else if (m[p2] != b2) bad = true;
            }
        }
        for (int p = base; p < base + 24 && !bad; p++)
          if (m[p] < 0 || (m[p] % 24) / 4 != (p % 24) / 4) bad = true;
        if (!bad) {
          for (int p = base; p < base + 24; p++) phi[p] = m[p];
          ok = true;
        }
      }
      if (!ok) return false;
    }
    return true;
  }
};

// Write a centre orbit back in physical quadrant order (inverse of readOrbit).
inline void writeOrbit(xcube &c, int x, int y, const uint8_t in[24]) {
  xorbit &o = c.rows[y - 1][x - 1];
  for (int f = 0; f < 6; f++)
    for (int oz = 0; oz < 4; oz++) o.a[4 * f + ((oz + c.facerots[f]) & 3)] = in[4 * f + oz];
}

// Effect of a pair algorithm (slices '2'/'3') on the two diagonal orbits it can
// touch: D2 on (a,a) and D3 on (b,b) for pair (a,b), in src form over the 24
// physical positions (after the alg, position k holds what was at src[k]).
// Computed on an 8x8x8 with (a,b) = (1,2); the layout is the same for every
// (a,b) (checked by diaglayouttest).
inline void diagEffect(const std::vector<mv> &alg, std::array<uint8_t, 24> &D2, std::array<uint8_t, 24> &D3) {
  xcube c(3, 3, 'U');
  for (auto m : alg) domove(m, c, true);  // depths 2/3 as written
  uint8_t o[24];
  readOrbit(c, 1, 1, o);
  for (int k = 0; k < 24; k++) D2[k] = o[k];
  readOrbit(c, 2, 2, o);
  for (int k = 0; k < 24; k++) D3[k] = o[k];
}
