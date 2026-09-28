// Wings-first solver (centres free, corners and middle edges preserved).
//
// State: brobdicube's wing orbits, one per depth k = 2..My+1 (index k-2),
// each a[slot] = piece, solved when a[i] == i.
//
// Moves: slice turns at wing depths (any face, any twist) are free; face turns
// are only used as setups on a stack and are later undone by their exact
// inverse, so the face turns of the whole sequence cancel as a word and corners
// and middle edges are untouched.  Centres get scrambled, which is fine: they
// are solved afterwards.
//
// The search works in the virtual frame V = the wing state as it would be if
// all open setups were closed now.  A slice S applied under open setup P acts on
// V as P^-1 S P; opening/closing a setup only changes P.  Goal: V solved and the
// stack empty.
#pragma once
#include "xcube.h"
#include "merge.h"
#include "wingbeam.h"
#include <array>
#include <unordered_set>
#include <atomic>

using Perm24 = std::array<uint8_t, 24>;  // src form: after the move, slot j holds what was at src[j]

struct WingGens {
  Perm24 face[6][4];   // face f, twist t (1..3): same action on every wing orbit
  Perm24 slice[6][4];  // slice at wing depth (any) from face f, twist t: acts on that depth's orbit
  static Perm24 compose(const Perm24 &first, const Perm24 &then) {  // apply first, then 'then'
    Perm24 r;
    for (int j = 0; j < 24; j++) r[j] = first[then[j]];
    return r;
  }
  static Perm24 ident() { Perm24 r; for (int j = 0; j < 24; j++) r[j] = j; return r; }
  void build() {
    for (int f = 0; f < 6; f++)
      for (int t = 1; t <= 3; t++) {
        xcube c(4, 4, 'U');  // N = 10: wing orbits at depths 2..5
        domove({1, (short)f, (short)t}, c, true);
        for (int j = 0; j < 24; j++) face[f][t][j] = c.wingorbits[0].a[j];
        xcube d(4, 4, 'U');
        domove({3, (short)f, (short)t}, d, true);  // depth 3 -> wing orbit 1
        for (int j = 0; j < 24; j++) slice[f][t][j] = d.wingorbits[1].a[j];
        // sanity: the slice touches only its own orbit, the face touches all alike
        for (int k = 0; k < 4; k++)
          for (int j = 0; j < 24; j++) {
            if (k != 1 && d.wingorbits[k].a[j] != j) error("slice touches another wing orbit");
            if (c.wingorbits[k].a[j] != c.wingorbits[0].a[j]) error("face acts differently on orbits");
          }
      }
  }
};

// Exact wing 3-cycles with centres free: for every 3-cycle of wing slots, the
// cheapest X [S,Y] X' or X [Y,S] X' (setup X of up to 2 face turns, slice S,
// face turn Y).  Used to finish the last few wings of each orbit.
struct WingFinisher {
  struct Macro { int cost = 1 << 20; std::vector<std::array<int, 3>> seq; };  // (isSlice, face, twist)
  // key: c[t2] = t1, c[t3] = t2, c[t1] = t3  (piece at t1 -> t2, t2 -> t3, t3 -> t1)
  std::vector<Macro> table = std::vector<Macro>(24 * 24 * 24);
  static int key(int t1, int t2, int t3) { return (t1 * 24 + t2) * 24 + t3; }
  void build(const WingGens &G) {
    static const int AX[6] = {0, 1, 2, 1, 2, 0};
    std::vector<std::vector<std::array<int, 2>>> setups = {{}};  // face turns
    for (int f = 0; f < 6; f++)
      for (int t = 1; t <= 3; t++) {
        setups.push_back({{f, t}});
        for (int g = 0; g < 6; g++)
          if (g != f)
            for (int u = 1; u <= 3; u++) setups.push_back({{f, t}, {g, u}});
      }
    for (auto &X : setups) {
      Perm24 P = WingGens::ident(), Pi = WingGens::ident();
      for (auto &m : X) P = WingGens::compose(P, G.face[m[0]][m[1]]);
      for (int i = (int)X.size() - 1; i >= 0; i--) Pi = WingGens::compose(Pi, G.face[X[i][0]][4 - X[i][1]]);
      for (int sf = 0; sf < 6; sf++)
        for (int st = 1; st <= 3; st++)
          for (int yf = 0; yf < 6; yf++) {
            if (AX[yf] == AX[sf]) continue;
            for (int yt = 1; yt <= 3; yt++)
              for (int type = 0; type < 2; type++) {
                const Perm24 &S = G.slice[sf][st], &Si = G.slice[sf][4 - st], &Y = G.face[yf][yt], &Yi = G.face[yf][4 - yt];
                Perm24 m = type == 0 ? WingGens::compose(WingGens::compose(WingGens::compose(S, Y), Si), Yi)
                                     : WingGens::compose(WingGens::compose(WingGens::compose(Y, S), Yi), Si);
                Perm24 c = WingGens::compose(WingGens::compose(P, m), Pi);
                int mv[24], n = 0;
                for (int j = 0; j < 24; j++) if (c[j] != j) { if (n < 24) mv[n] = j; n++; }
                if (n != 3) continue;
                // c[j] = source slot: piece at c[j] goes to j
                int a = mv[0];
                int t2 = a, t1 = c[a], t3 = -1;
                for (int q = 0; q < 3; q++) if (mv[q] != t1 && mv[q] != t2) t3 = mv[q];
                int cost = 4 + 2 * X.size();
                Macro M;
                M.cost = cost;
                for (auto &x : X) M.seq.push_back({0, x[0], x[1]});
                auto S1 = std::array<int, 3>{1, sf, st}, S2 = std::array<int, 3>{1, sf, 4 - st};
                auto Y1 = std::array<int, 3>{0, yf, yt}, Y2 = std::array<int, 3>{0, yf, 4 - yt};
                if (type == 0) { M.seq.push_back(S1); M.seq.push_back(Y1); M.seq.push_back(S2); M.seq.push_back(Y2); }
                else { M.seq.push_back(Y1); M.seq.push_back(S1); M.seq.push_back(Y2); M.seq.push_back(S2); }
                for (int i = (int)X.size() - 1; i >= 0; i--) M.seq.push_back({0, X[i][0], 4 - X[i][1]});
                // the same 3-cycle has three rotations of (t1,t2,t3)
                int ks[3] = {key(t1, t2, t3), key(t2, t3, t1), key(t3, t1, t2)};
                for (int kk : ks)
                  if (cost < table[kk].cost) table[kk] = M;
              }
          }
    }
    int have = 0;
    for (int a = 0; a < 24; a++) for (int b = 0; b < 24; b++) for (int c = 0; c < 24; c++)
      if (a != b && b != c && a != c && table[key(a, b, c)].cost < (1 << 20)) have++;
    if (have != 24 * 23 * 22) fprintf(stderr, "wing finisher covers %d of %d 3-cycles\n", have, 24 * 23 * 22);
  }
  // Table from a pool of 3-cycle algorithms (brobdicube's wing commutators,
  // wing.pool: 48-slot format, slots 0..23 = the wing orbit, slice '2').
  template <class PoolT> void buildFromPool(const PoolT &pool, std::vector<std::vector<std::array<int, 3>>> seqs) {
    for (size_t a = 0; a < pool.size(); a++) {
      const auto &c = pool.src[a];
      int mvd[3], n = 0;
      for (int j = 0; j < 24; j++) if (c[j] != j) { if (n < 3) mvd[n] = j; n++; }
      if (n != 3) continue;
      int t2 = mvd[0], t1 = c[t2], t3 = -1;
      for (int q = 0; q < 3; q++) if (mvd[q] != t1 && mvd[q] != t2) t3 = mvd[q];
      Macro M;
      M.cost = pool.cost[a];
      M.seq = seqs[a];
      int ks[3] = {key(t1, t2, t3), key(t2, t3, t1), key(t3, t1, t2)};
      for (int kk : ks)
        if (M.cost < table[kk].cost) table[kk] = M;
    }
    int have = 0;
    for (int a = 0; a < 24; a++) for (int b = 0; b < 24; b++) for (int c = 0; c < 24; c++)
      if (a != b && b != c && a != c && table[key(a, b, c)].cost < (1 << 20)) have++;
    if (have != 24 * 23 * 22) fprintf(stderr, "wing finisher covers %d of %d 3-cycles\n", have, 24 * 23 * 22);
  }
  // Solve one orbit (even permutation) by 3-cycles; append moves (slice depth dep).
  int finish(Perm24 v, int dep, std::vector<mv> &out, const WingGens *G = nullptr) const {
    int cost = 0;
    {  // odd permutation: one slice quarter turn first (3-cycles are even)
      Perm24 p = v;
      int par = 0;
      for (int i = 0; i < 24; i++) while (p[i] != i) { std::swap(p[i], p[p[i]]); par ^= 1; }
      if (par) {
        if (!G) return -1;
        Perm24 o = v;
        for (int j = 0; j < 24; j++) v[j] = o[G->slice[0][1][j]];
        out.push_back({dep, 0, 1});
        cost++;
      }
    }
    auto apply = [&](int t1, int t2, int t3) {  // piece at t1 -> t2, t2 -> t3, t3 -> t1
      const Macro &M = table[key(t1, t2, t3)];
      for (auto &s : M.seq) out.push_back({s[0] ? dep : 1, (short)s[1], (short)s[2]});
      cost += M.cost;
      uint8_t a = v[t1], b = v[t2], c = v[t3];
      v[t2] = a; v[t3] = b; v[t1] = c;
    };
    for (int i = 0; i < 24; i++)
      while (v[i] != i) {
        int t1 = 0;
        while (v[t1] != i) t1++;   // slot holding piece i
        int t2 = i, p = v[i], t3 = p;  // send the piece now at i to its own home
        if (t3 == t1) {  // 2-cycle (i t1): break it with another unsolved slot
          t3 = -1;
          for (int s = i + 1; s < 24 && t3 < 0; s++) if (v[s] != s && s != t1) t3 = s;
          if (t3 < 0) return -1;  // odd permutation left: should not happen
        }
        apply(t1, t2, t3);
      }
    return cost;
  }
};

// Beam search over (virtual wing state, setup stack).
//
// Actions, all in the current setup context P (the stack):
//   push / pop a face-turn setup                    (1 move)
//   single slice turn on one orbit                  (1 move; needed for parity)
//   commutator macro [S,X] = S X S' X' or [X,S] = X S X' S', applied at once to
//   every wing orbit where it helps: the face turns are shared, so on m orbits
//   it costs 2 + 2m moves (brobdicube's parallel trick, chosen per state).
// Score: moves so far + open setups (cost to close) + lambda * misplaced wings.
struct WingBeam {
  const WingGens *G;
  int N, My;
  int width = 16, maxStack = 2, maxWorse = 1;
  double lambda = 1.5;  // measured best on 4^3 wings (1.0 / 1.25 / 2 / 3 all worse), 2026-09-26
  const WingFinisher *fin = nullptr;
  const std::atomic<bool> *stop = nullptr;  // checked once per level
  int finishBelow = 12;  // try the exact finisher on nodes with at most this many wrong wings
  uint64_t seed = 0;     // seeded tie-break among equal scores (different runs explore different paths)

  struct Layer { int8_t axis; int16_t pos; int8_t amt; };
  static Layer layerOf(int N, int face, int dep, int twist) {
    static const int AX[6] = {0, 1, 2, 1, 2, 0}, POS[6] = {1, 1, 1, 0, 0, 0};
    Layer l;
    l.axis = AX[face];
    l.pos = POS[face] ? dep : N + 1 - dep;
    l.amt = POS[face] ? twist : (4 - twist) & 3;
    return l;
  }
  // last axis group of the emitted sequence, for exact merged move counts
  struct Tail {
    int8_t axis = -1;
    std::array<std::pair<int16_t, int8_t>, 64> L;
    int n = 0;
    int add(const Layer &m) {
      if (axis != m.axis) { axis = m.axis; n = 0; L[n++] = {m.pos, m.amt}; return 1; }
      for (int i = 0; i < n; i++)
        if (L[i].first == m.pos) {
          int a = (L[i].second + m.amt) & 3;
          if (a) { L[i].second = a; return 0; }
          L[i] = L[--n];
          if (n == 0) axis = -1;
          return -1;
        }
      if (n < 64) L[n++] = {m.pos, m.amt};
      return 1;
    }
  };
  struct Node {
    std::vector<uint8_t> V;    // My * 24 virtual wing state
    std::vector<uint8_t> stk;  // face*4 + twist
    Tail tail;
    int g, W, parent;
    std::vector<mv> moves;     // moves emitted by the action that led here
  };
  static uint64_t hashNode(const Node &n) {
    uint64_t h = 1469598103934665603ULL;
    for (auto c : n.V) h = (h ^ c) * 1099511628211ULL;
    for (auto c : n.stk) h = (h ^ (c + 101)) * 1099511628211ULL;
    return h;
  }
  Perm24 conj(const Perm24 &P, const Perm24 &op, const Perm24 &Pi) const {
    return WingGens::compose(WingGens::compose(P, op), Pi);
  }
  // base commutator macros (identity context): index (sf*4+st)*... -> perm
  struct Macro { int sf, st, xf, xt, type; Perm24 m; };
  std::vector<Macro> macros;
  void buildMacros() {
    static const int AX[6] = {0, 1, 2, 1, 2, 0};
    for (int sf = 0; sf < 6; sf++)
      for (int st = 1; st <= 3; st++)
        for (int xf = 0; xf < 6; xf++) {
          if (AX[xf] == AX[sf]) continue;
          for (int xt = 1; xt <= 3; xt++)
            for (int type = 0; type < 2; type++) {
              const Perm24 &S = G->slice[sf][st], &Si = G->slice[sf][4 - st];
              const Perm24 &X = G->face[xf][xt], &Xi = G->face[xf][4 - xt];
              Perm24 m = type == 0 ? WingGens::compose(WingGens::compose(WingGens::compose(S, X), Si), Xi)
                                   : WingGens::compose(WingGens::compose(WingGens::compose(X, S), Xi), Si);
              macros.push_back({sf, st, xf, xt, type, m});
            }
        }
  }

  std::vector<mv> solve(const std::vector<Perm24> &start0) {
    if (macros.empty()) buildMacros();
    // parity first: commutators are even, so each odd-parity orbit gets one
    // slice quarter turn up front (as brobdicube does)
    std::vector<Perm24> start = start0;
    std::vector<mv> prefix;
    for (int k = 0; k < My; k++) {
      Perm24 p = start[k];
      int par = 0;
      for (int i = 0; i < 24; i++)
        while (p[i] != i) { std::swap(p[i], p[p[i]]); par ^= 1; }
      if (par) {
        Perm24 old = start[k];
        for (int j = 0; j < 24; j++) start[k][j] = old[G->slice[0][1][j]];
        prefix.push_back({k + 2, 0, 1});
      }
    }
    std::vector<Node> nodes;
    Node root;
    root.V.resize(My * 24);
    root.W = 0;
    for (int k = 0; k < My; k++)
      for (int j = 0; j < 24; j++) {
        root.V[k * 24 + j] = start[k][j];
        root.W += start[k][j] != j;
      }
    root.g = 0;
    root.parent = -1;
    nodes.push_back(root);
    std::vector<int> frontier = {0};
    std::unordered_set<uint64_t> seen = {hashNode(root)};
    int best = 1 << 30, bestNode = -1;
    std::vector<mv> bestFinish;  // finisher moves appended after bestNode (after closing its stack)
    auto tryFinish = [&](int ni) {
      const Node &nd = nodes[ni];
      if (!fin || nd.W > finishBelow) return;
      std::vector<mv> tail;
      int cost = nd.g;
      for (int q = (int)nd.stk.size() - 1; q >= 0; q--) { tail.push_back({1, (short)(nd.stk[q] / 4), (short)(4 - nd.stk[q] % 4)}); cost++; }
      for (int k = 0; k < My; k++) {
        Perm24 v;
        bool solved = true;
        for (int j = 0; j < 24; j++) { v[j] = nd.V[k * 24 + j]; solved &= v[j] == j; }
        if (solved) continue;
        int c = fin->finish(v, k + 2, tail, G);
        if (c < 0) return;
        cost += c;
      }
      if (cost < best) { best = cost; bestNode = ni; bestFinish = tail; }
    };
    // candidate: kind 0 push(a=f,b=t), 1 pop, 2 single slice (a=f,b=t,c=orbit), 3 macro (a=sf*4+st, b=xf*4+xt, c=type, orbit mask in K)
    struct Cand { double score; int parent, kind, a, b, c, g, W; };
    for (int level = 0; level < 20000 && !frontier.empty(); level++) {
      if (stop && stop->load(std::memory_order_relaxed)) return {};
      std::vector<Cand> cs;
      for (int ni : frontier) {
        const Node &nd = nodes[ni];
        if (nd.W == 0 && nd.stk.empty()) {
          if (nd.g < best) { best = nd.g; bestNode = ni; bestFinish.clear(); }
          continue;
        }
        tryFinish(ni);
        int s = nd.stk.size();
        Perm24 P = WingGens::ident(), Pi = WingGens::ident();
        for (auto m : nd.stk) P = WingGens::compose(P, G->face[m / 4][m % 4]);
        for (int i = s - 1; i >= 0; i--) Pi = WingGens::compose(Pi, G->face[nd.stk[i] / 4][4 - nd.stk[i] % 4]);
        auto gainOf = [&](const Perm24 &c, int k) {
          const uint8_t *v = &nd.V[k * 24];
          int d = 0;
          for (int j = 0; j < 24; j++)
            if (c[j] != j) d += (v[c[j]] != j) - (v[j] != j);
          return -d;
        };
        if (s < maxStack)
          for (int f = 0; f < 6; f++) {
            if (s && nd.stk.back() / 4 == f) continue;
            for (int t = 1; t <= 3; t++) {
              Tail tl = nd.tail;
              int ng = nd.g + tl.add(layerOf(N, f, 1, t));
              cs.push_back({ng + (s + 1) + lambda * nd.W, ni, 0, f, t, 0, ng, nd.W});
            }
          }
        if (s) {
          int m = nd.stk.back();
          Tail tl = nd.tail;
          int ng = nd.g + tl.add(layerOf(N, m / 4, 1, 4 - m % 4));
          cs.push_back({ng + (s - 1) + lambda * nd.W, ni, 1, 0, 0, 0, ng, nd.W});
        }
        for (int sf = 0; sf < 6; sf++)
          for (int st = 1; st <= 3; st++) {
            Perm24 cS = conj(P, G->slice[sf][st], Pi);
            for (int k = 0; k < My; k++) {  // single slice
              int gk = gainOf(cS, k);
              if (-gk > maxWorse) continue;
              Tail tl = nd.tail;
              int ng = nd.g + tl.add(layerOf(N, sf, k + 2, st));
              cs.push_back({ng + s + lambda * (nd.W - gk), ni, 2, sf, st, k, ng, nd.W - gk});
            }
          }
        // commutator macros, conjugated into the current context; scored on their support only
        for (int mi = 0; mi < (int)macros.size(); mi++) {
          const Macro &M = macros[mi];
          Perm24 c = conj(P, M.m, Pi);
          int sup[24], ns = 0;
          for (int j = 0; j < 24; j++) if (c[j] != j) sup[ns++] = j;
          int tot = 0, nk = 0;
          for (int k = 0; k < My; k++) {
            const uint8_t *v = &nd.V[k * 24];
            int d = 0;
            for (int q = 0; q < ns; q++) { int j = sup[q]; d += (v[c[j]] != j) - (v[j] != j); }
            if (d < 0) { tot -= d; nk++; }
          }
          if (!nk) continue;
          // merged cost: slices of the helped orbits (same face/twist, commuting), face turn, inverses
          Tail tl = nd.tail;
          int ng = nd.g;
          auto emitS = [&](int tw) {
            for (int k = 0; k < My; k++) {
              const uint8_t *v = &nd.V[k * 24];
              int d = 0;
              for (int q = 0; q < ns; q++) { int j = sup[q]; d += (v[c[j]] != j) - (v[j] != j); }
              if (d < 0) ng += tl.add(layerOf(N, M.sf, k + 2, tw));
            }
          };
          if (M.type == 0) { emitS(M.st); ng += tl.add(layerOf(N, M.xf, 1, M.xt)); emitS(4 - M.st); ng += tl.add(layerOf(N, M.xf, 1, 4 - M.xt)); }
          else { ng += tl.add(layerOf(N, M.xf, 1, M.xt)); emitS(M.st); ng += tl.add(layerOf(N, M.xf, 1, 4 - M.xt)); emitS(4 - M.st); }
          cs.push_back({ng + s + lambda * (nd.W - tot), ni, 3, mi, 0, 0, ng, nd.W - tot});
        }
      }
      if (cs.empty()) break;
      if (bestNode >= 0 && level > 0) {
        double mn = 1e18;
        for (auto &c : cs) mn = std::min(mn, (double)c.g + nodes[c.parent].stk.size());
        if (mn >= best) break;
      }
      if (seed)  // scores are multiples of 0.5 or more apart; an offset below 1e-3 only breaks ties
        for (auto &c : cs) {
          uint64_t x = seed ^ ((uint64_t)c.parent << 32) ^ ((uint64_t)c.kind << 24) ^ ((uint64_t)c.a << 16) ^
                       ((uint64_t)c.b << 8) ^ (uint64_t)c.c;
          x *= 0x9E3779B97F4A7C15ULL; x ^= x >> 29; x *= 0xBF58476D1CE4E5B9ULL; x ^= x >> 32;
          c.score += 1e-3 * (double)(x & 0xFFFFF) / 0x100000;
        }
      size_t keep = std::min(cs.size(), (size_t)width * 8);
      std::partial_sort(cs.begin(), cs.begin() + keep, cs.end(), [](const Cand &x, const Cand &y) { return x.score < y.score; });
      frontier.clear();
      for (size_t i = 0; i < keep && (int)frontier.size() < width; i++) {
        const Cand &c = cs[i];
        if (c.g + (int)nodes[c.parent].stk.size() >= best) continue;
        Node n = nodes[c.parent];
        n.parent = c.parent;
        n.g = c.g;
        n.W = c.W;
        n.moves.clear();
        int s = n.stk.size();
        Perm24 P = WingGens::ident(), Pi = WingGens::ident();
        for (auto m : n.stk) P = WingGens::compose(P, G->face[m / 4][m % 4]);
        for (int q = s - 1; q >= 0; q--) Pi = WingGens::compose(Pi, G->face[n.stk[q] / 4][4 - n.stk[q] % 4]);
        auto applyV = [&](const Perm24 &cc, int k) {
          uint8_t *v = &n.V[k * 24], old[24];
          memcpy(old, v, 24);
          for (int j = 0; j < 24; j++) v[j] = old[cc[j]];
        };
        auto emit = [&](int dep, int f, int t) {
          n.moves.push_back({dep, (short)f, (short)t});
          n.tail.add(layerOf(N, f, dep, t));
        };
        if (c.kind == 0) { n.stk.push_back(c.a * 4 + c.b); emit(1, c.a, c.b); }
        else if (c.kind == 1) { int m = n.stk.back(); n.stk.pop_back(); emit(1, m / 4, 4 - m % 4); }
        else if (c.kind == 2) { applyV(conj(P, G->slice[c.a][c.b], Pi), c.c); emit(c.c + 2, c.a, c.b); }
        else {
          const Macro &M = macros[c.a];
          Perm24 cc = conj(P, M.m, Pi);
          std::vector<int> K;
          for (int k = 0; k < My; k++) {
            const uint8_t *v = &n.V[k * 24];
            int d = 0;
            for (int j = 0; j < 24; j++) if (cc[j] != j) d += (v[cc[j]] != j) - (v[j] != j);
            if (d < 0) K.push_back(k);
          }
          for (auto k : K) applyV(cc, k);
          auto emitS = [&](int tw) { for (auto k : K) emit(k + 2, M.sf, tw); };
          if (M.type == 0) { emitS(M.st); emit(1, M.xf, M.xt); emitS(4 - M.st); emit(1, M.xf, 4 - M.xt); }
          else { emit(1, M.xf, M.xt); emitS(M.st); emit(1, M.xf, 4 - M.xt); emitS(4 - M.st); }
        }
        if (!seen.insert(hashNode(n)).second) continue;
        nodes.push_back(std::move(n));
        frontier.push_back(nodes.size() - 1);
      }
      if (getenv("WDBG") && level % 10 == 0 && !frontier.empty())
        fprintf(stderr, "level %d frontier %zu best-W %d g %d\n", level, frontier.size(), nodes[frontier[0]].W, nodes[frontier[0]].g);
    }
    std::vector<mv> out = prefix;
    if (bestNode < 0) return {};
    std::vector<int> path;
    for (int n = bestNode; nodes[n].parent >= 0; n = nodes[n].parent) path.push_back(n);
    std::reverse(path.begin(), path.end());
    for (int n : path) out.insert(out.end(), nodes[n].moves.begin(), nodes[n].moves.end());
    out.insert(out.end(), bestFinish.begin(), bestFinish.end());
    return out;
  }
};

// Solve all wing orbits in independent chunks of `chunk` depths.  A chunk's
// sequence is the identity on every other depth's wings (setups are closed,
// and commutator face turns cancel on depths the commutator is not applied
// to), so chunks can be solved separately -- in parallel -- and concatenated.
// Time is linear in N; sharing of setups stays within a chunk.
#include <thread>
#include <atomic>
inline std::vector<mv> solveWingsChunked(const WingGens &G, const WingFinisher &F, int N, int My,
                                         const std::vector<Perm24> &start, int width, int chunk, int threads,
                                         bool merge = true, const std::atomic<bool> *stop = nullptr,
                                         uint64_t seed = 0, const WingTableBeam *tb = nullptr) {
  int nc = (My + chunk - 1) / chunk;
  std::vector<std::vector<mv>> out(nc);
  std::atomic<int> nx{0};
  std::atomic<bool> fail{false};
  std::atomic<int> nFallback{0};
  std::vector<std::thread> th;
  for (int t = 0; t < std::max(1, threads); t++)
    th.emplace_back([&] {
      for (int c; (c = nx++) < nc;) {
        if (stop && stop->load(std::memory_order_relaxed)) { fail = true; break; }
        int k0 = c * chunk, k1 = std::min(My, k0 + chunk);
        std::vector<Perm24> st(start.begin() + k0, start.begin() + k1);
        if (tb && k1 - k0 == 1) {  // table beam over the wing pool (one depth)
          bool ok;
          auto idx = tb->solve(st[0], width, ok, seed ? seed * 0x9E3779B97F4A7C15ULL + c + 1 : 0);
          if (ok) {
            std::vector<mv> m;
            for (auto r : idx)
              for (auto &x : tb->poolOf(r).movesOf(r & ~WingTableBeam::FIN))
                m.push_back({(short)(x[0] == 1 ? 1 : 2 + k0), (short)x[1], (short)x[2]});
            out[c] = std::move(m);
            continue;
          }
          // the table beam can dead-end at narrow widths: fall back to the
          // setup/commutator beam, which always completes (3-cycle finisher)
          if (stop && stop->load(std::memory_order_relaxed)) { fail = true; continue; }
          nFallback++;
        }
        WingBeam B{&G, N, k1 - k0};
        B.width = width;
        B.fin = &F;
        B.finishBelow = 1 << 20;
        B.stop = stop;
        B.seed = seed ? seed * 0x9E3779B97F4A7C15ULL + c + 1 : 0;
        auto m = B.solve(st);
        if (m.empty() && k1 > k0) { fail = true; continue; }
        for (auto &x : m)
          if (x.dep > 1) x.dep += k0;  // local orbit k -> depth k+2; shift to the chunk's depths
        out[c] = std::move(m);
      }
    });
  for (auto &t : th) t.join();
  std::vector<mv> all;
  if (fail) return all;
  if (!merge) {  // no cross-depth merging, but simplify each depth's sequence (linear)
    for (int c = 0; c < nc; c++) {
      Merger mg(N);
      mg.add(out[c], {});
      long r, m;
      auto s = mg.merge(r, m);
      // Check: this depth's sequence solves its own wing orbit, and its face
      // turns reduce to the empty word, so every other depth, the corners and
      // the middle edges are untouched.
      if (chunk == 1) {
        Perm24 v = start[c];
        // face-turn word reduced modulo same-axis commutation: axis groups of
        // (amount on the U/L/F face, amount on the opposite face)
        static const int AXf[6] = {0, 1, 2, 1, 2, 0}, SIDEf[6] = {0, 0, 0, 1, 1, 1};
        struct AG { int axis, amt[2]; };
        std::vector<AG> word;
        for (mv x : s) {
          if (x.dep > N + 1 - x.dep) { x.dep = N + 1 - x.dep; x.face = oppface[x.face]; x.twist = (4 - x.twist) & 3; }
          const Perm24 &g = x.dep == 1 ? G.face[x.face][x.twist] : G.slice[x.face][x.twist];
          if (x.dep != 1 && x.dep != c + 2) { if (getenv("WCHK")) fprintf(stderr, "chunk %d: slice at depth %d\n", c, x.dep); return {}; }
          Perm24 o = v;
          for (int j = 0; j < 24; j++) v[j] = o[g[j]];
          if (x.dep == 1) {
            int ax = AXf[x.face], sd = SIDEf[x.face];
            if (word.empty() || word.back().axis != ax) word.push_back({ax, {0, 0}});
            word.back().amt[sd] = (word.back().amt[sd] + x.twist) & 3;
            if (!word.back().amt[0] && !word.back().amt[1]) word.pop_back();
          }
        }
        for (int j = 0; j < 24; j++) if (v[j] != j) { if (getenv("WCHK")) fprintf(stderr, "chunk %d: orbit not solved\n", c); return {}; }
        if (!word.empty()) { if (getenv("WCHK")) fprintf(stderr, "chunk %d: face word not empty (%zu)\n", c, word.size()); return {}; }
      }
      all.insert(all.end(), s.begin(), s.end());
    }
    return all;
  }
  // Split each chunk's sequence into closed blocks (face turns reduce to the
  // empty word at the block end, so every setup is undone) and let the merger
  // interleave blocks of different chunks to cancel setups at the seams.
  // Blocks of one chunk stay in order.
  Merger mg(N);
  for (int c = 0; c < nc; c++) {
    std::vector<mv> blk;
    std::vector<std::pair<int, int>> word;  // (face, net twist) free reduction of face turns
    for (auto &m : out[c]) {
      blk.push_back(m);
      if (m.dep == 1) {
        if (!word.empty() && word.back().first == m.face) {
          word.back().second = (word.back().second + m.twist) & 3;
          if (!word.back().second) word.pop_back();
        } else word.push_back({m.face, m.twist});
      }
      if (word.empty()) { mg.add(blk, {{c, ~0ULL}}); blk.clear(); }
    }
    if (!blk.empty()) mg.add(blk, {{c, ~0ULL}});
  }
  long raw, merged;
  return mg.merge(raw, merged);
}
