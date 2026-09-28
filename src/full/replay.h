// Cache-friendly replay of a move sequence onto an xcube's centre orbits.
//
// domove() applies each slice turn by walking a whole row and a whole column
// of the centre grid; for a column that is one cache miss per row, so long
// slice-heavy sequences are slow on big cubes.  Face turns are lazy in xcube
// (facerots only).  Here we record the facerots value in force at every slice
// turn, then replay the sequence row by row (in parallel): row y only sees
// its own column element of every slice, plus the full-row effect of slices
// at its own depth.  The rotation tables are brobdicube's (xcube.cpp domove).
//
// Only the centre grid and facerots are updated: the caller is responsible
// for wings, corners and middle edges (e.g. a wing solution leaves corners and
// middle edges untouched and the wings solved).  Middle-layer turns are not
// supported (not used by the wing phase).
#pragma once
#include "xcube.h"
#include <atomic>
#include <thread>
#include <vector>

struct CentreReplay {
  struct Rot { int kind; int f[4], o[4]; };  // kind 4 = quarter cycle, 2 = double swap
  // [face*3 + twist] for non-middle slice turns: horizontal and vertical parts
  static const Rot &hor(int k) {
    static const Rot H[19] = {{},
      {4, {1, 2, 3, 4}, {0, 0, 0, 0}}, {2, {1, 2, 3, 4}, {0, 0, 0, 0}}, {4, {4, 3, 2, 1}, {0, 0, 0, 0}},
      {4, {2, 0, 4, 5}, {3, 3, 1, 3}}, {2, {2, 5, 4, 0}, {3, 3, 1, 3}}, {4, {2, 5, 4, 0}, {3, 3, 1, 3}},
      {4, {0, 1, 5, 3}, {2, 1, 0, 3}}, {2, {0, 3, 5, 1}, {2, 3, 0, 1}}, {4, {0, 3, 5, 1}, {2, 3, 0, 1}},
      {4, {2, 5, 4, 0}, {1, 1, 3, 1}}, {2, {2, 5, 4, 0}, {1, 1, 3, 1}}, {4, {2, 0, 4, 5}, {1, 1, 3, 1}},
      {4, {0, 3, 5, 1}, {0, 1, 2, 3}}, {2, {0, 3, 5, 1}, {0, 1, 2, 3}}, {4, {0, 1, 5, 3}, {0, 3, 2, 1}},
      {4, {4, 3, 2, 1}, {2, 2, 2, 2}}, {2, {1, 2, 3, 4}, {2, 2, 2, 2}}, {4, {4, 1, 2, 3}, {2, 2, 2, 2}}};
    return H[k];
  }
  static const Rot &ver(int k) {
    static const Rot V[19] = {{},
      {4, {1, 2, 3, 4}, {1, 1, 1, 1}}, {2, {1, 2, 3, 4}, {1, 1, 1, 1}}, {4, {4, 3, 2, 1}, {1, 1, 1, 1}},
      {4, {2, 0, 4, 5}, {0, 0, 2, 0}}, {2, {2, 5, 4, 0}, {0, 0, 2, 0}}, {4, {2, 5, 4, 0}, {0, 0, 2, 0}},
      {4, {0, 1, 5, 3}, {3, 2, 1, 0}}, {2, {0, 3, 5, 1}, {3, 0, 1, 2}}, {4, {0, 3, 5, 1}, {3, 0, 1, 2}},
      {4, {2, 5, 4, 0}, {2, 2, 0, 2}}, {2, {2, 5, 4, 0}, {2, 2, 0, 2}}, {4, {2, 0, 4, 5}, {2, 2, 0, 2}},
      {4, {0, 3, 5, 1}, {1, 2, 3, 0}}, {2, {0, 3, 5, 1}, {1, 2, 3, 0}}, {4, {0, 1, 5, 3}, {1, 0, 3, 2}},
      {4, {4, 3, 2, 1}, {3, 3, 3, 3}}, {2, {1, 2, 3, 4}, {3, 3, 3, 3}}, {4, {4, 1, 2, 3}, {3, 3, 3, 3}}};
    return V[k];
  }
  static inline void rot(xorbit &xo, const Rot &r, const uint8_t fr[6]) {
    uchar &c1 = xo.a[4 * r.f[0] + ((r.o[0] + fr[r.f[0]]) & 3)];
    uchar &c2 = xo.a[4 * r.f[1] + ((r.o[1] + fr[r.f[1]]) & 3)];
    uchar &c3 = xo.a[4 * r.f[2] + ((r.o[2] + fr[r.f[2]]) & 3)];
    uchar &c4 = xo.a[4 * r.f[3] + ((r.o[3] + fr[r.f[3]]) & 3)];
    if (r.kind == 4) { auto t = c1; c1 = c2; c2 = c3; c3 = c4; c4 = t; }
    else { std::swap(c1, c3); std::swap(c2, c4); }
  }

  // Apply mvs to the centre grid and facerots of cube.
  static void apply(xcube &cube, const std::vector<mv> &mvs, int threads) {
    struct Slice { int d; uint8_t k; uint8_t fr[6]; };
    std::vector<Slice> sl;
    uint8_t fr[6];
    for (int f = 0; f < 6; f++) fr[f] = cube.facerots[f];
    const int lim = cube.Mx + cube.My + 3;
    for (mv m : mvs) {
      if (m.dep > lim - m.dep) {  // same normalisation as domove
        m.dep = lim - m.dep;
        m.face = oppface[m.face];
        m.twist = (4 - m.twist) & 3;
      }
      if (m.dep == 1) { fr[m.face] = (fr[m.face] + 4 - m.twist) & 3; continue; }
      int d = m.dep - 2;
      if (d >= cube.My) error("CentreReplay: middle-layer turns not supported");
      Slice s;
      s.d = d;
      s.k = m.face * 3 + m.twist;
      for (int f = 0; f < 6; f++) s.fr[f] = fr[f];
      sl.push_back(s);
    }
    const int My = cube.My, Mx = cube.Mx;
    std::atomic<int> nx{0};
    std::vector<std::thread> th;
    for (int t = 0; t < std::max(1, threads); t++)
      th.emplace_back([&] {
        for (int y; (y = nx++) < My;) {
          auto &row = cube.rows[y];
          for (const Slice &s : sl) {
            if (s.d == y)
              for (int x = 0; x < Mx; x++) rot(row[x], hor(s.k), s.fr);
            rot(row[s.d], ver(s.k), s.fr);
          }
        }
      });
    for (auto &t : th) t.join();
    for (int f = 0; f < 6; f++) cube.facerots[f] = fr[f];
  }
};
