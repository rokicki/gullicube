// Symmetry-class pools: one representative per class of the 48 cube symmetries
// x inverse, expanded on load.
//
// Every pool here is closed under the cube symmetries (and inversion), and
// keeps the cheapest algorithm per permutation, so a class is determined by any
// member and all its images cost the same.  A class file holds "COST moves"
// lines (placeholder depths as in the expanded pool); loading regenerates the
// images.
//
// Orbit model: how single moves act on the orbit's positions (src form: after
// the move, position j holds what was at src[j]).  For each symmetry s there is
// a relabelling sigma of the positions with perm(s(A)) = sigma o perm(A) o
// sigma^-1 for every algorithm A; it is derived once from the move perms.
#pragma once
#include "validate.h"
#include "../pair/seq.h"
#include "../pair/pool.h"
#include "wingbeam.h"
#include <array>
#include <functional>
#include <string>
#include <vector>

struct OrbitModel {
  int P = 24;                                          // positions (24 single orbit, 48 pair)
  std::vector<std::vector<uint8_t>> mp;                // [face*16 + depth*4 + twist] -> src perm
  std::vector<Sym> syms;                               // 48 symmetries (validate.h)
  std::vector<std::vector<uint8_t>> sigma, sigmaInv;   // per symmetry
  int maxDepth = 2;                                    // placeholder depths used (1..maxDepth)

  const std::vector<uint8_t> &mv(int f, int d, int t) const {
    const auto &r = mp[f * 16 + d * 4 + t];
    if (r.empty()) { fprintf(stderr, "OrbitModel: move depth %d not in the model\n", d); exit(1); }
    return r;
  }
  // perm of a move list (placeholder depths)
  std::vector<uint8_t> permOf(const std::vector<Move> &a) const {
    std::vector<uint8_t> p(P), r(P);
    for (int k = 0; k < P; k++) p[k] = k;
    for (auto &m : a) {
      const auto &q = mv(m.face, m.depth, m.amt);
      for (int k = 0; k < P; k++) r[k] = p[q[k]];
      p.swap(r);
    }
    return p;
  }
  // symmetry image of a move list (img = sym * 2 + inverse)
  std::vector<Move> image(std::vector<Move> b, int img) const {
    if (img & 1) {
      std::reverse(b.begin(), b.end());
      for (auto &m : b) m.amt = (4 - m.amt) & 3;
    }
    const Sym &s = syms[img >> 1];
    for (auto &m : b) {
      m.face = s.face[m.face];
      if (s.refl) m.amt = (4 - m.amt) & 3;
    }
    return b;
  }
  // perm of image img of an algorithm with perm p (by relabelling, no replay)
  void imagePerm(const uint8_t *p, int img, uint8_t *out) const {
    const auto &sg = sigma[img >> 1], &si = sigmaInv[img >> 1];
    uint8_t q[48];
    if (img & 1) for (int k = 0; k < P; k++) q[p[k]] = k;  // inverse
    else for (int k = 0; k < P; k++) q[k] = p[k];
    // perm(s(A)) = sigma o q o sigma^-1 in src form: out[sigma[j]] = sigma[q[j]]
    for (int j = 0; j < P; j++) out[sg[j]] = sg[q[j]];
    (void)si;
  }
  // derive sigma for every symmetry from the move perms; returns false if some
  // symmetry has no consistent relabelling
  bool build(std::function<std::vector<uint8_t>(int f, int d, int t)> movePerm, int positions, int depths) {
    P = positions;
    maxDepth = depths;
    mp.assign(6 * 16, std::vector<uint8_t>());
    for (int f = 0; f < 6; f++)
      for (int d = 1; d <= depths; d++)
        for (int t = 1; t <= 3; t++) mp[f * 16 + d * 4 + t] = movePerm(f, d, t);
    syms = symmetries();
    sigma.clear();
    sigmaInv.clear();
    for (auto &s : syms) {
      // sigma with src_{s(m)}[sigma[j]] = sigma[src_m[j]] for every move m;
      // positions form one component per orbit, so seed each component in
      // turn (depth-first over the seed choices) and propagate
      std::vector<int> sg(P, -1);
      std::function<bool()> solve = [&]() -> bool {
        int j0 = -1;
        for (int j = 0; j < P; j++) if (sg[j] < 0) { j0 = j; break; }
        if (j0 < 0) return true;
        std::vector<char> used(P, 0);
        for (int j = 0; j < P; j++) if (sg[j] >= 0) used[sg[j]] = 1;
        for (int x0 = 0; x0 < P; x0++) {
          if (used[x0]) continue;
          std::vector<int> saved = sg;
          sg[j0] = x0;
          std::vector<int> stack = {j0};
          bool ok = true;
          while (!stack.empty() && ok) {
            int j = stack.back();
            stack.pop_back();
            for (int f = 0; f < 6 && ok; f++)
              for (int d = 1; d <= depths && ok; d++)
                for (int t = 1; t <= 3 && ok; t++) {
                  const auto &a = mv(f, d, t);
                  int tt = s.refl ? (4 - t) & 3 : t;
                  const auto &b = mv(s.face[f], d, tt);
                  int from = a[j], to = b[sg[j]];
                  if (sg[from] < 0) {
                    for (int k = 0; k < P; k++) if (sg[k] == to) ok = false;  // injective
                    if (!ok) break;
                    sg[from] = to;
                    stack.push_back(from);
                  } else if (sg[from] != to) ok = false;
                }
          }
          if (ok && solve()) return true;
          sg = saved;
        }
        return false;
      };
      if (!solve()) return false;
      std::vector<uint8_t> S(P), SI(P);
      for (int j = 0; j < P; j++) { S[j] = sg[j]; SI[sg[j]] = j; }
      for (int f = 0; f < 6; f++)  // final check on every move
        for (int d = 1; d <= depths; d++)
          for (int t = 1; t <= 3; t++) {
            const auto &a = mv(f, d, t);
            int tt = s.refl ? (4 - t) & 3 : t;
            const auto &b = mv(s.face[f], d, tt);
            for (int j = 0; j < P; j++) if (S[a[j]] != b[S[j]]) return false;
          }
      sigma.push_back(S);
      sigmaInv.push_back(SI);
    }
    return true;
  }
  // class key: least image perm (lexicographic) over the 96 images
  void classKey(const uint8_t *p, uint8_t *key) const {
    uint8_t cur[48];
    bool first = true;
    for (int img = 0; img < 96; img++) {
      imagePerm(p, img, cur);
      if (first || std::lexicographical_compare(cur, cur + P, key, key + P)) { std::copy(cur, cur + P, key); first = false; }
    }
  }
};

// Orbit models built from xcube (every sticker distinct).  Placeholder depth 2
// is the orbit's depth; for mids, 3 is the middle layer.
inline OrbitModel diagModel() {  // diagonal orbit (1,1) of 8^3; placeholders 2 and 3 are both depth 2
  OrbitModel m;
  m.build([](int f, int d, int t) {
    xcube c(3, 3, 'U');
    domove({(short)(d == 1 ? 1 : 2), (short)f, (short)t}, c, true);
    uint8_t o[24];
    readOrbit(c, 1, 1, o);
    return std::vector<uint8_t>(o, o + 24);
  }, 24, 3);
  return m;
}
inline OrbitModel midModel() {  // mid orbit (4,1) of 9^3: '2' -> depth 2, '3' -> middle (5)
  OrbitModel m;
  m.build([](int f, int d, int t) {
    xcube c(4, 3, 'U');
    domove({(short)(d == 1 ? 1 : d == 2 ? 2 : 5), (short)f, (short)t}, c, true);
    uint8_t o[24];
    readOrbit(c, 4, 1, o);
    return std::vector<uint8_t>(o, o + 24);
  }, 24, 3);
  return m;
}
inline OrbitModel pairModel() {  // the pair solver's 48 positions: orbits (1,2)/(2,1) of 8^3 via Layout::phi
  static Layout L;
  static bool built = L.build();
  (void)built;
  OrbitModel m;
  m.build([](int f, int d, int t) {
    xcube c(3, 3, 'U');
    domove({(short)d, (short)f, (short)t}, c, true);  // placeholder depth d = actual depth (2 = a+1, 3 = b+1)
    uint8_t X0[24], X1[24];
    readOrbit(c, 1, 2, X0);
    readOrbit(c, 2, 1, X1);
    // position p (xcube position phi[p]) now holds the piece from xcube position lab; back to pair positions
    std::vector<int> inv(48);
    for (int p = 0; p < 48; p++) inv[L.phi[p]] = p;
    std::vector<uint8_t> src(48);
    for (int p = 0; p < 48; p++) {
      int x = L.phi[p];
      int lab = x < 24 ? X0[x] : 24 + X1[x - 24];
      src[p] = inv[lab];
    }
    return src;
  }, 48, 3);
  return m;
}
inline OrbitModel wingModel() {  // wing orbit 0 (depth 2) of 6^3
  OrbitModel m;
  m.build([](int f, int d, int t) {
    xcube c(2, 2, 'U');
    domove({(short)d, (short)f, (short)t}, c, true);
    return std::vector<uint8_t>(c.wingorbits[0].a, c.wingorbits[0].a + 24);
  }, 24, 2);
  return m;
}

// A class file: representatives (cost + moves).
struct ClassPool {
  std::vector<uint8_t> cost;
  std::vector<uint8_t> dflag;  // pair pools: 1 if the algorithm moves diagonal pieces (symmetry-invariant)
  std::vector<std::vector<Move>> rep;
  bool load(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) return false;
    char buf[4096];
    while (fgets(buf, sizeof buf, f)) {
      if (buf[0] == '#') continue;
      char *p;
      int c = (int)strtol(buf, &p, 10);
      if (p == buf) continue;
      std::string t(p);
      while (!t.empty() && (t.back() == '\n' || t.back() == ' ')) t.pop_back();
      uint8_t df = 1;
      if (t.size() >= 3 && t.compare(t.size() - 3, 2, " D") == 0 && (t.back() == '0' || t.back() == '1')) {
        df = t.back() - '0';
        t.resize(t.size() - 3);
      }
      cost.push_back(c);
      dflag.push_back(df);
      rep.push_back(parse_alg(t));
    }
    fclose(f);
    return true;
  }
  size_t size() const { return cost.size(); }
  // visit every distinct image of every class: fn(perm, cost, ref = class * 96 + img)
  // (images equal to an earlier image of the same class are skipped)
  // With distinct = false, images that coincide (classes with symmetry) are
  // visited more than once -- fine for hash-table inserts, and much faster.
  template <class F> void expand(const OrbitModel &m, F fn, bool distinct = true) const {
    expandRange(m, 0, rep.size(), fn, distinct);
  }
  template <class F> void expandRange(const OrbitModel &m, size_t lo, size_t hi, F fn, bool distinct = true) const {
    uint8_t p[48], cur[48];
    std::vector<uint64_t> seen;
    for (size_t r = lo; r < hi; r++) {
      auto base = m.permOf(rep[r]);
      std::copy(base.begin(), base.end(), p);
      seen.clear();
      for (int img = 0; img < 96; img++) {
        m.imagePerm(p, img, cur);
        if (distinct) {
          uint64_t h = 1469598103934665603ULL;
          for (int k = 0; k < m.P; k++) h = (h ^ cur[k]) * 1099511628211ULL;
          if (std::find(seen.begin(), seen.end(), h) != seen.end()) continue;  // 64-bit hash of the image
          seen.push_back(h);
        }
        fn(cur, cost[r], (uint32_t)(r * 96 + img));
      }
    }
  }
  // move text of image `ref` (for output), placeholder depths
  std::string text(const OrbitModel &m, uint32_t ref) const {
    return groups_str(to_groups(m.image(rep[ref / 96], ref % 96)));
  }
};

// The image `ref` (= class * 96 + image) of a class pool as a pool entry of dst
// (single orbit: positions 24..47 identity); returns its index.
inline uint32_t appendImage(const ClassPool &cp, const OrbitModel &m, uint32_t ref, Pool &dst) {
  auto base = m.permOf(cp.rep[ref / 96]);
  uint8_t p[48];
  m.imagePerm(base.data(), ref % 96, p);
  std::array<uint8_t, 48> s;
  for (int k = 0; k < 48; k++) s[k] = k < m.P ? p[k] : k;
  dst.addEntry(cp.cost[ref / 96], s, cp.text(m, ref), cp.dflag[ref / 96]);
  return dst.cost.size() - 1;
}

// A wing pool from a symmetry-class file: every distinct image with cost <=
// maxCost; each entry keeps its scan data and a reference (class * 96 + image);
// moves are generated on demand from the class representative.
inline void loadWingPoolFromClasses(WingPool &wp, const ClassPool &cp, const OrbitModel &m, int maxCost = 99) {
  cp.expand(m, [&](const uint8_t *p, int c, uint32_t ref) {
    if (c > maxCost) return;
    std::array<uint8_t, 24> s;
    std::copy(p, p + 24, s.begin());
    WingPool::Moves mv3;
    for (auto &x : m.image(cp.rep[ref / 96], ref % 96)) mv3.push_back({(int8_t)x.depth, (int8_t)x.face, (int8_t)x.amt});
    wp.add(c, s, mv3);
    wp.ref.push_back(ref);
  });
  wp.movesFn = [&cp, &m](uint32_t ref) {
    WingPool::Moves mv3;
    for (auto &x : m.image(cp.rep[ref / 96], ref % 96)) mv3.push_back({(int8_t)x.depth, (int8_t)x.face, (int8_t)x.amt});
    return mv3;
  };
}

// A Pool with every distinct image of a class file (single orbit: positions
// P..47 identity).
inline void poolFromClasses(Pool &dst, const ClassPool &cp, const OrbitModel &m) {
  cp.expand(m, [&](const uint8_t *p, int c, uint32_t ref) {
    std::array<uint8_t, 48> s;
    for (int k = 0; k < 48; k++) s[k] = k < m.P ? p[k] : k;
    dst.addEntry(c, s, m.image(cp.rep[ref / 96], ref % 96), cp.dflag[ref / 96]);
  });
}
inline bool endsWith(const std::string &s, const char *suf) {
  size_t n = strlen(suf);
  return s.size() >= n && s.compare(s.size() - n, n, suf) == 0;
}
