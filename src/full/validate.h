// Validation helpers shared by mkpools and mine: symmetries, algorithm text,
// full-cube snapshots, and the "changes nothing but this centre orbit" check.
#pragma once
#include "layout.h"
#include <functional>
#include <string>
#include <vector>

// 48 cube symmetries as face maps + handedness
struct Sym { int face[6]; bool refl; };
inline std::vector<Sym> symmetries() {
  std::vector<Sym> r;
  int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  int psign[6] = {1, -1, -1, 1, 1, -1};
  for (int pi = 0; pi < 6; pi++)
    for (int s = 0; s < 8; s++) {
      Sym y;
      for (int f = 0; f < 6; f++) {
        int v[3];
        for (int k = 0; k < 3; k++) v[k] = ps::NORMAL[f][perms[pi][k]] * ((s >> k) & 1 ? -1 : 1);
        for (int g = 0; g < 6; g++)
          if (ps::NORMAL[g][0] == v[0] && ps::NORMAL[g][1] == v[1] && ps::NORMAL[g][2] == v[2]) y.face[f] = g;
      }
      y.refl = psign[pi] * (__builtin_popcount(s) & 1 ? -1 : 1) < 0;
      r.push_back(y);
    }
  return r;
}

inline std::string algText(const std::vector<mv> &a) {
  std::string s;
  for (auto m : a) {
    if (!s.empty()) s += ' ';
    if (m.dep > 1) s += std::to_string(m.dep);
    s += "ULFRBD"[m.face];
    if (m.twist == 2) s += '2';
    if (m.twist == 3) s += '\'';
  }
  return s;
}

// snapshot of everything in an xcube, for "did anything else change" checks
struct Snap {
  std::vector<uint8_t> centres, wings, rest;
  explicit Snap(const xcube &c) {
    for (int y = 1; y <= c.My; y++)
      for (int x = 1; x <= c.Mx; x++) {
        uint8_t o[24];
        readOrbit(c, x, y, o);
        centres.insert(centres.end(), o, o + 24);
      }
    for (auto &w : c.wingorbits) wings.insert(wings.end(), w.a, w.a + 24);
    rest.insert(rest.end(), c.tc.corner.a, c.tc.corner.a + 24);
    rest.insert(rest.end(), c.tc.edge.a, c.tc.edge.a + 24);
    rest.insert(rest.end(), c.tc.centercenters, c.tc.centercenters + 6);
  }
};

// apply alg with depth map; return perm on centre orbit (x,y) if only it changed
inline bool centreOnly(int Mx, int My, const std::vector<mv> &a, std::function<int(int)> dm, int x, int y,
                       std::vector<uint8_t> &perm) {
  xcube c(Mx, My, 'U');
  Snap s0(c);
  applyMapped(c, a, dm);
  Snap s1(c);
  if (s0.wings != s1.wings || s0.rest != s1.rest) return false;
  int idx = ((y - 1) * Mx + (x - 1)) * 24;
  for (size_t i = 0; i < s0.centres.size(); i++)
    if ((i < (size_t)idx || i >= (size_t)idx + 24) && s0.centres[i] != s1.centres[i]) return false;
  perm.assign(s1.centres.begin() + idx, s1.centres.begin() + idx + 24);
  return true;
}

// Wing-phase check: alg (placeholder depth 2 -> depth dm(2)) changes nothing but
// wing orbit k (0-based, depth k+2) among corners, middle edges and all wing
// orbits (centres are free in the wing phase).  perm = the permutation of orbit k.
inline bool wingOnly(int Mx, int My, const std::vector<mv> &a, std::function<int(int)> dm, int k,
                     std::vector<uint8_t> &perm) {
  xcube c(Mx, My, 'U');
  xcube c0 = c;
  applyMapped(c, a, dm);
  for (int j = 0; j < 24; j++)
    if (c.tc.corner.a[j] != c0.tc.corner.a[j] || c.tc.edge.a[j] != c0.tc.edge.a[j]) return false;
  for (size_t w = 0; w < c.wingorbits.size(); w++)
    if ((int)w != k)
      for (int j = 0; j < 24; j++)
        if (c.wingorbits[w].a[j] != c0.wingorbits[w].a[j]) return false;
  perm.assign(c.wingorbits[k].a, c.wingorbits[k].a + 24);
  return true;
}
