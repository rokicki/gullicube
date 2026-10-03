// Portable orbit solvers (no SIMD intrinsics; plain loops the compiler may
// vectorise).
//
// ColourSolver: level beam over a pool of algorithms acting on up to 48
//   coloured positions (target colour (p%24)/4) -- used for oblique pairs, and
//   for diagonal / mid orbits whose pools are padded to 48 with the identity.
//   Optional endgame: at every beam node with few wrong pieces, enumerate one
//   small-support algorithm and look the result up in a one-algorithm finish
//   table (eg2.h).
// WingSolver:   level beam over a pool acting on 24 distinct wing pieces.
//
// Both rank candidates by  4*LAMBDA*correct - effcost  (quarter moves), where
// effcost charges face-turn setup moves (the X of X c X') at a discount, since
// those mostly cancel when many orbit solutions are merged.
#pragma once
#include <mutex>
#include <memory>
#include <cmath>
#include "mintree.h"
#include "../pair/pool.h"
#include "../pair/eg2.h"
#include <unordered_set>
#include <atomic>
#include <map>
#include <cstddef>
#include <thread>
#include <chrono>

struct Weights {
  int setupFace = 4, other = 4;  // quarter moves
};

inline int setupDepth(const std::vector<Move> &m) {
  int k = 0;
  while (2 * (k + 1) <= (int)m.size()) {
    const Move &a = m[k], &b = m[m.size() - 1 - k];
    if (a.depth != 1 || b.depth != 1 || a.face != b.face || (a.amt + b.amt) % 4) break;
    k++;
  }
  return k;
}
inline int effCost4(const std::vector<Move> &m, const Weights &w) {
  int sd = setupDepth(m), c = 0;
  for (int i = 0; i < (int)m.size(); i++) c += (i < sd || i >= (int)m.size() - sd) ? w.setupFace : w.other;
  return c;
}
// the algorithm's moves from its same-axis groups (the order groups_str prints)
inline std::vector<Move> poolMoves(const Pool &p, uint32_t a) {
  static const int FACEOF[3][2] = {{0, 5}, {1, 3}, {2, 4}};
  std::vector<Move> r;
  for (const uint16_t *g = p.gbegin(a); g != p.gend(a); g++)
    for (int slot = 0; slot < 6; slot++) {
      int amt = (gamt(*g) >> (2 * slot)) & 3;
      if (amt) r.push_back(Move{(uint8_t)FACEOF[gaxis(*g)][slot / 3], (uint8_t)(slot % 3 + 1), (uint8_t)amt});
    }
  return r;
}

// Algorithm reference in a solution: which pool, which index.
struct AlgRef {
  uint8_t pool;  // 0 = main pool, 1 = endgame enumeration pool, 2 = endgame finish pool
  uint32_t idx;
};

struct ColourSolver {
  const Pool *pool = nullptr;
  const EG2 *eg = nullptr;
  Weights wt;
  int lambda = 3, width = 16, maxMis = 8;
  const std::atomic<bool> *stop = nullptr;  // checked once per level
  // structure-of-arrays masks and scores
  size_t n = 0;
  std::vector<uint64_t> S[6];
  std::vector<int> eff4;
  // endgame algorithms' effective costs are only needed for the few that get
  // chosen, so compute them on demand rather than for all ~26M at startup
  int egEff4(const Pool &p, uint32_t a) const { return effCost4(poolMoves(p, a), wt); }

  void build(const Pool &p, const EG2 *e, const Weights &w) {
    pool = &p;
    eg = e;
    wt = w;
    n = p.size();
    for (int c = 0; c < 6; c++) {
      S[c].resize(n);
      for (size_t a = 0; a < n; a++) S[c][a] = p.S[a][c];
    }
    eff4.resize(n);
    for (size_t a = 0; a < n; a++) eff4[a] = effCost4(poolMoves(p, a), w);
  }

  // top-k algorithms by 4*lambda*correct - eff4
  void topK(const State &s, int k, std::vector<std::pair<int, uint32_t>> &out) const {
    out.clear();
    int thr = INT32_MIN;
    const uint64_t M0 = s.M[0], M1 = s.M[1], M2 = s.M[2], M3 = s.M[3], M4 = s.M[4], M5 = s.M[5];
    const int mul = 4 * lambda;
    const uint64_t *s0 = S[0].data(), *s1 = S[1].data(), *s2 = S[2].data(), *s3 = S[3].data(),
                   *s4 = S[4].data(), *s5 = S[5].data();
    const int *e4 = eff4.data();
    constexpr size_t B = 64;
    int sc[B];
    for (size_t base = 0; base < n; base += B) {
      size_t m = std::min(B, n - base);
      int mx = INT32_MIN;
      for (size_t j = 0; j < m; j++) {
        size_t a = base + j;
        // the S_c are disjoint (each source feeds one destination, of one
        // target colour) and so are the M_c, so the six popcounts collapse to one
        int c = __builtin_popcountll((s0[a] & M0) | (s1[a] & M1) | (s2[a] & M2) | (s3[a] & M3) |
                                     (s4[a] & M4) | (s5[a] & M5));
        sc[j] = mul * c - e4[a];
        mx = std::max(mx, sc[j]);
      }
      if (mx <= thr) continue;
      for (size_t j = 0; j < m; j++)
        if (sc[j] > thr) {
          out.push_back({sc[j], (uint32_t)(base + j)});
          if ((int)out.size() > 2 * k) {
            std::nth_element(out.begin(), out.begin() + k - 1, out.end(),
                             [](auto &x, auto &y) { return x.first > y.first; });
            out.resize(k);
            thr = std::min_element(out.begin(), out.end())->first;
          }
        }
    }
    std::sort(out.begin(), out.end(), [](auto &x, auto &y) { return x.first > y.first; });
    if ((int)out.size() > k) out.resize(k);
  }

  struct Node { State s; int g4, parent; uint32_t alg; };

  // Returns the chosen algorithms (empty and ok=false if unsolved).
  std::vector<AlgRef> solve(const State &start, bool &ok, int widthOverride = 0) const {
    ok = false;
    const int width = widthOverride ? widthOverride : this->width;
    std::vector<Node> nodes = {{start, 0, -1, 0}};
    std::vector<int> frontier = {0};
    std::unordered_set<uint64_t> seen = {start.hash()};
    int best4 = 1 << 30, bestNode = -1;
    EG2::Sol bestSol;
    std::vector<std::pair<int, uint32_t>> top;
    for (int depth = 0; depth < 60 && !frontier.empty(); depth++) {
      if (stop && stop->load(std::memory_order_relaxed)) return {};
      struct C { int key, g4, parent; uint32_t a; };
      std::vector<C> cs;
      for (int ni : frontier) {
        const Node &nd = nodes[ni];
        int m = 48 - nd.s.correct();
        if (m == 0) {
          if (nd.g4 < best4) { best4 = nd.g4; bestNode = ni; bestSol = EG2::Sol(); }
          continue;
        }
        if (eg && m <= maxMis) {
          auto sol = eg->solve(nd.s, 1, 200);
          if (sol.cost < 200) {
            int e = 0;
            for (auto x : sol.path) e += egEff4(eg->cand, x);
            if (sol.finAlg >= 0) e += egEff4(eg->fin, sol.finAlg);
            if (nd.g4 + e < best4) { best4 = nd.g4 + e; bestNode = ni; bestSol = sol; }
          }
        }
        topK(nd.s, 2 * width, top);
        for (auto &t : top) {
          int ng4 = nd.g4 + eff4[t.second];
          if (ng4 >= best4) continue;
          cs.push_back({t.first - nd.g4, ng4, ni, t.second});
        }
      }
      std::sort(cs.begin(), cs.end(), [](const C &x, const C &y) { return x.key > y.key; });
      frontier.clear();
      for (auto &c : cs) {
        State ns = nodes[c.parent].s.apply(*pool, c.a);
        if (!seen.insert(ns.hash()).second) continue;
        nodes.push_back({ns, c.g4, c.parent, c.a});
        frontier.push_back(nodes.size() - 1);
        if ((int)frontier.size() >= width) break;
      }
    }
    std::vector<AlgRef> r;
    if (bestNode < 0) return r;
    for (int nn = bestNode; nodes[nn].parent >= 0; nn = nodes[nn].parent) r.push_back({0, nodes[nn].alg});
    std::reverse(r.begin(), r.end());
    for (auto x : bestSol.path) r.push_back({1, x});
    if (bestSol.finAlg >= 0) r.push_back({2, (uint32_t)bestSol.finAlg});
    ok = true;
    return r;
  }
};

// Wings: 24 distinct pieces; state a[slot] = piece; solved when a[i] == i.
struct WingSolver {
  const Pool *pool = nullptr;
  Weights wt;
  int lambda = 3, width = 16;
  std::vector<int> eff4;
  std::vector<std::vector<uint8_t>> supp;
  void build(const Pool &p, const Weights &w) {
    pool = &p;
    wt = w;
    eff4.resize(p.size());
    supp.resize(p.size());
    for (size_t a = 0; a < p.size(); a++) {
      eff4[a] = effCost4(poolMoves(p, a), w);
      for (int d = 0; d < 24; d++)
        if (p.src[a][d] != d) supp[a].push_back(d);
    }
  }
  using WS = std::array<uint8_t, 24>;
  static int fixedCount(const WS &s) {
    int c = 0;
    for (int i = 0; i < 24; i++) c += s[i] == i;
    return c;
  }
  WS apply(const WS &s, uint32_t a) const {
    WS r;
    for (int d = 0; d < 24; d++) r[d] = s[pool->src[a][d]];
    return r;
  }
  static uint64_t hashOf(const WS &s) {
    uint64_t h = 1469598103934665603ULL;
    for (auto c : s) h = (h ^ c) * 1099511628211ULL;
    return h;
  }
  std::vector<AlgRef> solve(const WS &start, bool &ok) const {
    struct Node { WS s; int g4, parent; uint32_t alg; };
    ok = false;
    std::vector<Node> nodes = {{start, 0, -1, 0}};
    std::vector<int> frontier = {0};
    std::unordered_set<uint64_t> seen = {hashOf(start)};
    int best4 = 1 << 30, bestNode = -1;
    const int mul = 4 * lambda;
    for (int depth = 0; depth < 60 && !frontier.empty(); depth++) {
      struct C { int key, g4, parent; uint32_t a; };
      std::vector<C> cs;
      for (int ni : frontier) {
        const Node &nd = nodes[ni];
        int f0 = fixedCount(nd.s);
        if (f0 == 24) {
          if (nd.g4 < best4) { best4 = nd.g4; bestNode = ni; }
          continue;
        }
        std::vector<C> local;
        for (uint32_t a = 0; a < pool->size(); a++) {
          int f = f0;
          for (auto d : supp[a]) f += (nd.s[pool->src[a][d]] == d) - (nd.s[d] == d);
          int g4 = nd.g4 + eff4[a];
          if (g4 >= best4 || f <= f0 - 1) continue;  // wing algorithms are 3-cycles: require no loss
          local.push_back({mul * f - g4, g4, ni, a});
        }
        int k = std::min<int>(2 * width, local.size());
        std::partial_sort(local.begin(), local.begin() + k, local.end(), [](auto &x, auto &y) { return x.key > y.key; });
        cs.insert(cs.end(), local.begin(), local.begin() + k);
      }
      std::sort(cs.begin(), cs.end(), [](const C &x, const C &y) { return x.key > y.key; });
      frontier.clear();
      for (auto &c : cs) {
        WS ns = apply(nodes[c.parent].s, c.a);
        if (!seen.insert(hashOf(ns)).second) continue;
        nodes.push_back({ns, c.g4, c.parent, c.a});
        frontier.push_back(nodes.size() - 1);
        if ((int)frontier.size() >= width) break;
      }
    }
    std::vector<AlgRef> r;
    if (bestNode < 0) return r;
    for (int nn = bestNode; nodes[nn].parent >= 0; nn = nodes[nn].parent) r.push_back({0, nodes[nn].alg});
    std::reverse(r.begin(), r.end());
    ok = true;
    return r;
  }
};

// Guaranteed fallback for coloured orbits: decompose into 3-cycles from a
// complete table (every 3-cycle of each orbit present).  Each step picks a
// wrong position p, a wrong position q (same orbit) holding p's target colour,
// and a position r whose target is the colour now at p; the 3-cycle
// q->p->r->q fixes p and leaves r correct, so the number of wrong positions
// strictly decreases.  Such q and r always exist (4 identical pieces per
// colour), so this always terminates with the orbit(s) solved.
struct CycleFinisher {
  const Pool *pool = nullptr;
  std::vector<int32_t> table = std::vector<int32_t>(48 * 48 * 48, -1);  // key (q,p,r) -> alg
  static int key(int q, int p, int r) { return (q * 48 + p) * 48 + r; }
  void build(const Pool &p) {
    pool = &p;
    for (size_t a = 0; a < p.size(); a++) {
      int m[48], n = 0;
      for (int j = 0; j < 48; j++) if (p.src[a][j] != j) m[n++] = j;
      if (n != 3) continue;
      const auto &c = p.src[a];  // piece at c[j] goes to j
      int P = m[0], Q = c[P], R = -1;
      for (int t = 0; t < 3; t++) if (c[m[t]] == P) R = m[t];
      int ks[3] = {key(Q, P, R), key(P, R, Q), key(R, Q, P)};
      for (int k : ks)
        if (table[k] < 0 || p.cost[a] < p.cost[table[k]]) table[k] = a;
    }
  }
  // positions [lo, hi) form one orbit; returns false only if the table is incomplete
  bool finishOrbit(State &s, int lo, int hi, std::vector<AlgRef> &out, uint8_t poolId) const {
    for (int guard = 0; guard < 64; guard++) {
      int p = -1;
      for (int j = lo; j < hi && p < 0; j++) if (s.col[j] != target(j)) p = j;
      if (p < 0) return true;
      int q = -1, r = -1;
      for (int j = lo; j < hi && q < 0; j++) if (j != p && s.col[j] != target(j) && s.col[j] == target(p)) q = j;
      if (q < 0) return false;
      // best r: wrong, wants p's colour, holds q's target; then wrong wanting p's colour; then any wanting it
      for (int pass = 0; pass < 3 && r < 0; pass++)
        for (int j = lo; j < hi && r < 0; j++) {
          if (j == p || j == q || target(j) != s.col[p]) continue;
          bool wrong = s.col[j] != target(j);
          if (pass == 0 && wrong && s.col[j] == target(q)) r = j;
          else if (pass == 1 && wrong) r = j;
          else if (pass == 2) r = j;
        }
      if (r < 0) return false;
      int a = table[key(q, p, r)];
      if (a < 0) return false;
      s = s.apply(*pool, a);
      out.push_back({poolId, (uint32_t)a});
    }
    return false;
  }
};

// Level-by-level beam with a direct-mapped table (no top-k, no sort, no
// collision resolution).  Every (state, algorithm) child is scored; a child
// whose score beats the current weakest in the table has its hash computed
// incrementally (Zobrist over the algorithm's moved positions, without
// building the child) and is dropped into slot hash & (T-1), replacing the
// occupant only if better.  The occupied slots form the next level.  A second
// direct-mapped table of earlier-level hashes blocks revisits.  T = 1 is plain
// greedy.
// Compact finish table: Zobrist hash (the beam's keys) of every colour state a
// single algorithm of `pool` solves -> the cheapest such algorithm and its
// effective cost.  Lets the beam see, for every child it admits (whose hash it
// already has), whether one more algorithm finishes it and at exactly what
// cost -- "perfect" within the algorithm set.
#ifndef FINISH_RING
#define FINISH_RING 32
#endif
struct FinishTable {
  struct Slot { uint64_t key; uint32_t alg; int16_t e4; uint8_t pid, cost; };  // one probe, one line
  static_assert(sizeof(Slot) == 16);
  std::vector<Slot> own;        // the table (built, or read from a saved file)
  const Slot *slots = nullptr;  // own.data()
  uint64_t nslots = 0;          // any size below 2^32: the home slot is a multiply-shift of the key's high half
  size_t used = 0;
  size_t home(uint64_t k) const { return (size_t)(((k >> 32) * nslots) >> 32); }
  size_t step(size_t i) const { return i + 1 == nslots ? 0 : i + 1; }
  const FinishTable *next = nullptr;  // a second table probed as well (cheapest hit wins)
  FinishTable() = default;
  FinishTable(const FinishTable &) = delete;
  void init(size_t capacity) {  // capacity: the most states it will hold
    initSlots(capacity + capacity / 4 + 16);
  }
  void initSlots(size_t sz) {
    own.assign(sz, Slot{0, 0, 0, 0, 0});
    slots = own.data();
    nslots = sz;
    used = 0;
  }
  // add every algorithm of p (effective costs eff4) under pool id `id`
  void add(const Pool &p, const uint64_t (&Z)[48][6], const std::vector<int> &eff4, uint8_t id) {
    // the state b solves: position src[d] holds colour target(d)
    for (size_t b = 0; b < p.size(); b++) insert(p.src[b], Z, eff4[b], p.cost[b], id, b);
  }
  static uint64_t keyOf(const std::array<uint8_t, 48> &src, const uint64_t (&Z)[48][6]) {
    uint64_t k = 0;  // the state the algorithm solves: position src[d] holds colour target(d)
    for (int d = 0; d < 48; d++) k ^= Z[src[d]][target(d)];
    return k;
  }
  void insert(const std::array<uint8_t, 48> &src, const uint64_t (&Z)[48][6], int e, int c, uint8_t id, uint32_t b) {
    insertKey(keyOf(src, Z), e, c, id, b);
  }
  // Bulk inserts: each key's slot is prefetched and the insert itself is done
  // RING inserts later, when the cache line has arrived (random writes into a
  // table of hundreds of MB are otherwise one cache miss each).
  struct Batch {
    static constexpr int RING = FINISH_RING;
    FinishTable &t;
    Slot ring[RING];
    int head = 0, n = 0;  // circular queue of pending inserts
    bool atomic = false;
    size_t claimed = 0;  // new slots (atomic mode)
    explicit Batch(FinishTable &t, bool atomic = false) : t(t), atomic(atomic) {}
    void put(const Slot &o) {
      if (atomic) claimed += t.insertKeyAtomic(o.key, o.e4, o.cost, o.pid, o.alg);
      else t.insertKey(o.key, o.e4, o.cost, o.pid, o.alg);
    }
    ~Batch() { flush(); }
    void add(uint64_t k, int e, int c, uint8_t id, uint32_t b) {
      __builtin_prefetch(&t.own[t.home(k)], 1);
      if (n == RING) {  // the oldest pending insert's line has had RING inserts' time to arrive
        put(ring[head]);
        head = (head + 1) % RING;
        n--;
      }
      ring[(head + n) % RING] = {k, b, (int16_t)e, id, (uint8_t)c};
      n++;
    }
    void flush() {
      for (; n > 0; n--, head = (head + 1) % RING) put(ring[head]);
    }
  };
  // Thread-safe insert (parallel builds): the key is claimed with a CAS, and
  // the value word (alg, e4, pid, cost) is replaced with a CAS while cheaper.
  // returns 1 if a new slot was claimed
  int insertKeyAtomic(uint64_t k, int e, int c, uint8_t id, uint32_t b) {
    int claimed = 0;
    static_assert(sizeof(Slot) == 16 && offsetof(Slot, alg) == 8);
    Slot nv{k, b, (int16_t)e, id, (uint8_t)c};
    uint64_t nval;
    memcpy(&nval, (const char *)&nv + 8, 8);
    for (size_t i = home(k);; i = step(i)) {
      std::atomic_ref<uint64_t> key(own[i].key);
      uint64_t cur = key.load(std::memory_order_acquire);
      if (cur == 0) {
        uint64_t zero = 0;
        if (!key.compare_exchange_strong(zero, k)) { cur = zero; if (cur != k) continue; }
        else claimed = 1;
      } else if (cur != k) continue;
      std::atomic_ref<uint64_t> val(*(uint64_t *)((char *)&own[i] + 8));
      uint64_t old = val.load(std::memory_order_relaxed);
      for (;;) {
        Slot os;
        memcpy((char *)&os + 8, &old, 8);
        // keep the cheapest, ties to the lower (pool, algorithm): the result does not depend on thread timing
        if (old != 0 && (os.e4 < e || (os.e4 == e && (os.pid < id || (os.pid == id && os.alg <= b))))) return claimed;
        if (val.compare_exchange_weak(old, nval)) return claimed;
      }
    }
  }
  void insertKey(uint64_t k, int e, int c, uint8_t id, uint32_t b) {
    for (size_t i = home(k);; i = step(i)) {
      Slot &s = own[i];
      if (s.key == k) {
        if (e < s.e4 || (e == s.e4 && (id < s.pid || (id == s.pid && b < s.alg)))) s = {k, b, (int16_t)e, id, (uint8_t)c};
        return;
      }
      if (s.key == 0) { s = {k, b, (int16_t)e, id, (uint8_t)c}; used++; return; }
    }
  }
  const Slot *find1(uint64_t k) const {
    for (size_t i = home(k);; i = step(i)) {
      if (slots[i].key == k) return &slots[i];
      if (slots[i].key == 0) return nullptr;
    }
  }
  // cheapest finishing algorithm for the state with hash k (this table and next), or nullptr
  const Slot *find(uint64_t k) const {
    const Slot *a = slots ? find1(k) : nullptr, *b = next ? next->find(k) : nullptr;
    return !a ? b : !b ? a : (b->e4 < a->e4 ? b : a);
  }

};


// Finish tables from a symmetry-class pool (classpool.h): `ext` gets every
// image of every class (pool id `pid`, algorithm = class * 96 + image, effective
// cost 4 * moves); `own` holds the beam's pool (pool id 0) and chains to ext.
// Built in parallel on every load (0.3 s for the 49M-image diagonal classes
// on 16 threads: keys straight from the class representatives, prefetched
// batched inserts), so there is no cache file.
template <class CP, class OM>
inline void buildFinishFromClasses(FinishTable &own, FinishTable &ext, const Pool &pool, const uint64_t (&Z)[48][6],
                                   const std::vector<int> &eff4, const CP &cp, const OM &model, uint8_t pid) {
  buildFinishFromClassList(own, ext, pool, Z, eff4, std::vector<std::pair<const CP *, uint8_t>>{{&cp, pid}}, model);
}
template <class CP, class OM>
inline void buildFinishFromClassList(FinishTable &own, FinishTable &ext, const Pool &pool, const uint64_t (&Z)[48][6],
                                     const std::vector<int> &eff4, const std::vector<std::pair<const CP *, uint8_t>> &list,
                                     const OM &model) {
  own.init(pool.size());
  own.add(pool, Z, eff4, 0);
  own.next = &ext;
  {
    // Keys straight from the class representative: image = sigma o q o sigma^-1
    // (q = the representative's perm or its inverse), so its key
    // XOR_d Z[image[d]][target(d)] = XOR_j Z[sigma[q[j]]][target(sigma[j])]
    // = XOR_j KT_s[q[j]][j]; the image perm is never built.  Positions P..47
    // are the identity: their part is a constant.
    const int P = model.P;
    uint64_t k0 = 0;
    for (int d = P; d < 48; d++) k0 ^= Z[d][target(d)];
    std::vector<uint64_t> KT(48 * P * P);
    for (int sy = 0; sy < 48; sy++)
      for (int x = 0; x < P; x++)
        for (int j = 0; j < P; j++) KT[(sy * P + x) * P + j] = Z[model.sigma[sy][x]][target(model.sigma[sy][j])];
    const int threads = std::max(1u, std::thread::hardware_concurrency());
    // every image's key, on all threads: emit(threadState, key, cost, pid, alg)
    auto forEachKey = [&](auto makeState, auto emit, auto finish) {
      for (auto &e : list) {
        const CP &cp = *e.first;
        const uint8_t pid = e.second;
        std::atomic<size_t> next{0};
        std::vector<std::thread> th;
        for (int t = 0; t < threads; t++)
          th.emplace_back([&] {
            auto st = makeState();
            uint8_t q[2][48];
            for (size_t lo; (lo = next.fetch_add(4096)) < cp.size();)
              for (size_t r = lo; r < std::min(cp.size(), lo + 4096); r++) {
                auto base = model.permOf(cp.rep[r]);
                for (int j = 0; j < P; j++) { q[0][j] = base[j]; q[1][base[j]] = j; }
                const int c = cp.cost[r];
                for (int img = 0; img < 96; img++) {
                  const uint8_t *qq = q[img & 1];
                  const uint64_t *kt = &KT[(size_t)(img >> 1) * P * P];
                  uint64_t k = k0;
                  for (int j = 0; j < P; j++) k ^= kt[qq[j] * P + j];
                  emit(*st, k, c, pid, (uint32_t)(r * 96 + img));
                }
              }
            finish(*st);
          });
        for (auto &t : th) t.join();
      }
    };
    // Pass 1: the number of distinct states (many images finish the same one),
    // estimated with a HyperLogLog sketch (2^14 registers, about 1% error), so
    // the table is sized for the states it will hold.  The keys are Zobrist
    // XORs, so their bits are already uniform.
    constexpr int HB = 14;
    std::vector<uint8_t> reg(1 << HB, 0);
    std::mutex regMu;
    forEachKey([] { return std::make_unique<std::vector<uint8_t>>(1 << HB, 0); },
               [](std::vector<uint8_t> &r, uint64_t k, int, uint8_t, uint32_t) {
                 const uint64_t rest = k << HB;
                 const uint8_t rho = rest ? (uint8_t)(__builtin_clzll(rest) + 1) : (uint8_t)(64 - HB + 1);
                 uint8_t &x = r[k >> (64 - HB)];
                 if (rho > x) x = rho;
               },
               [&](std::vector<uint8_t> &r) {
                 std::lock_guard<std::mutex> g(regMu);
                 for (size_t i = 0; i < reg.size(); i++) reg[i] = std::max(reg[i], r[i]);
               });
    double sum = 0;
    size_t zeros = 0;
    for (auto x : reg) { sum += std::ldexp(1.0, -x); zeros += x == 0; }
    const double m = reg.size(), alpha = 0.7213 / (1 + 1.079 / m);
    double est = alpha * m * m / sum;
    if (est < 2.5 * m && zeros) est = m * std::log(m / zeros);  // small range: linear counting
    // Pass 2: insert, the table about 74% full (estimate x 1.35, with 3% margin for the estimate)
    ext.initSlots((size_t)(est * 1.03 * 1.35) + 64);
    std::atomic<size_t> claimedAll{0};
    forEachKey([&] { return std::make_unique<FinishTable::Batch>(ext, true); },
               [](FinishTable::Batch &b, uint64_t k, int c, uint8_t pid, uint32_t alg) { b.add(k, 4 * c, c, pid, alg); },
               [&](FinishTable::Batch &b) { b.flush(); claimedAll += b.claimed; });
    ext.used = claimedAll;
  }
}

// Optional refinement of the beam score: the base score is mul * correct - g4
// (an estimate of 12/4 = 3 moves per misplaced piece); a heuristic subtracts a
// non-negative state-dependent extra (quarter-move units) from it.  The base
// score is then an upper bound, so the cheap pre-filters stay valid and the
// extra is only evaluated (on the built child) for candidates that pass them.
struct PairHeur {
  virtual ~PairHeur() = default;
  virtual int extra4(const State &s) const = 0;
};

struct TableBeam {
  const ColourSolver *cs = nullptr;  // pool, SoA masks, eff4, endgame, weights
  uint64_t Z[48][6];
  std::vector<std::vector<uint8_t>> supp;  // moved positions per algorithm
  // movement index: for each movement s -> d (piece at s lands on d), the
  // algorithms containing it (CSR: moff[d*48+s] .. moff[d*48+s+1] in mdat)
  std::vector<uint32_t> moff, mdat;
  bool useIndex = false;  // movement-index scoring (measured slower than the bit-sliced scan)
  // Bit-sliced algorithms: for each source position s, the target colour of the
  // destination it lands on (3 bits, as planes T0..T2 over the 48 positions).
  // With the state's colours as planes C0..C2, the correct count after the
  // algorithm is popcount(~((C0^T0)|(C1^T1)|(C2^T2)) & 48 bits).
  std::vector<uint64_t> T0, T1, T2;
  static constexpr uint64_t M48 = (1ULL << 48) - 1;
  // Face-leave masks: positions whose piece the algorithm moves onto another
  // face.  With keepSolved, an algorithm is rejected (one AND) if it would
  // carry any currently-correct piece off its face.
  std::vector<uint64_t> X;
  bool keepSolved = true;
  bool strictKeep = false;
  // experiment: with keepSolved, levels 1..rootLevels score only the rootTop
  // algorithms that scored best at the start (0 = off)
  uint32_t rootTop = 0;
  int rootLevels = 0;
  int firstWidth = 0;  // experiment: the start keeps this many children (0 = the width)
  double widthGrowth = 1;  // experiment: the table at depth d holds width * widthGrowth^d nodes ...
  int widthCap = 1 << 16;  // ... up to this
  const std::atomic<bool> *stop = nullptr;  // checked once per level  // reject any algorithm that moves a solved piece at all
  bool timing = false;
  bool noIncr = getenv("NOINCR") != nullptr, noAoS = getenv("NOAOS") != nullptr, noList = getenv("NOLIST") != nullptr;
  // Transposed: bitset over algorithms, per position p, of the algorithms that
  // carry p's piece off its face.  Rejected = OR over correct positions.
  size_t nw = 0;
  std::vector<uint64_t> BP;  // 48 * nw words
  struct alignas(32) Rec { uint64_t t0, t1, t2; int32_t e4, pad; };
  std::vector<Rec> rec;  // per algorithm: bit-planes + cost in one cache line
  mutable std::atomic<long> statScanned{0}, statSurvived{0}, statAdmit{0}, statLevels{0}, statSolves{0};
  mutable std::atomic<long> nsFilter{0}, nsScore{0}, nsNodePre{0}, nsSeam{0}, nsExpand{0}, nWeakScan{0}, nsAdmitFin{0}, nAdmitFin{0}, nFromList{0}, nRestricted{0}, nBitset{0}, nRebuild{0};
  // with timing, per node by its number of correct positions: nodes, algorithms
  // surviving the keep-solved filter (each fully scored), and admitted (score >= weakest)
  mutable std::atomic<long> histNodes[49] = {}, histSurv[49] = {}, histAdmit[49] = {};
  // with timing, per level (depth): nodes, scored, correct positions summed, and wall time
  mutable std::atomic<long> lvNodes[60] = {}, lvSurv[60] = {}, lvCorrect[60] = {}, lvNs[60] = {};
  static long nsNow() { return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
  void build(const ColourSolver &c) {
    cs = &c;
    std::mt19937_64 r(4242);
    for (auto &z : Z) for (auto &x : z) x = r();
    const Pool &p = *c.pool;
    supp.resize(p.size());
    for (size_t a = 0; a < p.size(); a++)
      for (int d = 0; d < 48; d++)
        if (p.src[a][d] != d) supp[a].push_back(d);
    T0.assign(p.size(), 0); T1.assign(p.size(), 0); T2.assign(p.size(), 0);
    X.assign(p.size(), 0);
    for (size_t a = 0; a < p.size(); a++)
      for (int d = 0; d < 48; d++) {
        int sp = p.src[a][d];  // the piece at sp lands on d
        if (strictKeep ? sp != d : (sp % 24) / 4 != (d % 24) / 4) X[a] |= 1ULL << sp;
      }
    nw = (p.size() + 63) / 64;
    BP.assign(48 * nw, 0);
    for (size_t a = 0; a < p.size(); a++)
      for (int q = 0; q < 48; q++)
        if (X[a] >> q & 1) BP[q * nw + a / 64] |= 1ULL << (a % 64);
    for (size_t a = 0; a < p.size(); a++)
      for (int d = 0; d < 48; d++) {
        int sp = p.src[a][d], t = target(d);
        if (t & 1) T0[a] |= 1ULL << sp;
        if (t & 2) T1[a] |= 1ULL << sp;
        if (t & 4) T2[a] |= 1ULL << sp;
      }
    // seam data: first / last group per algorithm, distinct first groups
    {
      std::map<uint16_t, uint16_t> hmap;
      heads.clear();
      headIdx.resize(p.size());
      lastG.resize(p.size());
      nGroups.resize(p.size());
      zeroIdx.assign(p.size(), 0);
      for (size_t a = 0; a < p.size(); a++) {
        uint16_t h = p.head[a];
        auto it = hmap.find(h);
        if (it == hmap.end()) { it = hmap.emplace(h, heads.size()).first; heads.push_back(h); }
        headIdx[a] = it->second;
        nGroups[a] = p.gend(a) - p.gbegin(a);
        lastG[a] = nGroups[a] ? *(p.gend(a) - 1) : 0xFFFF;
      }
    }
    rec.resize(p.size());  // after the planes are filled
    for (size_t a = 0; a < p.size(); a++) rec[a] = {T0[a], T1[a], T2[a], c.eff4[a], 0};
    std::vector<uint32_t> cnt(48 * 48, 0);
    for (size_t a = 0; a < p.size(); a++)
      for (auto d : supp[a]) cnt[d * 48 + p.src[a][d]]++;
    moff.assign(48 * 48 + 1, 0);
    for (int m = 0; m < 48 * 48; m++) moff[m + 1] = moff[m] + cnt[m];
    mdat.resize(moff.back());
    std::vector<uint32_t> fill(moff.begin(), moff.end() - 1);
    for (size_t a = 0; a < p.size(); a++)
      for (auto d : supp[a]) mdat[fill[d * 48 + p.src[a][d]]++] = a;
  }
  uint64_t zob(const State &s) const {
    uint64_t h = 0;
    for (int d = 0; d < 48; d++) h ^= Z[d][s.col[d]];
    return h;
  }
  struct Node { State s; uint64_t h; int g4, parent; uint32_t alg; uint16_t tail = 0xFFFF; };
  // Seam-aware costs: an algorithm's cost minus the moves that cancel against
  // the last same-axis group of the sequence so far (only that group can
  // interact with the algorithm's first group).  First groups come from a
  // small set, so per node the saving is computed once per distinct first
  // group (bonus4[headIdx[a]]).
  bool seamAware = false;
  std::vector<uint16_t> zeroIdx;  // all zero (seam-unaware lookups)
  std::vector<uint16_t> heads;    // distinct first groups
  std::vector<uint16_t> headIdx;  // per algorithm: index into heads
  std::vector<uint16_t> lastG;    // per algorithm: last group
  std::vector<uint8_t> nGroups;
  // moves saved (>= 0) when a sequence ending in group t is followed by one
  // starting with group f; merged = the combined group (0xFFFF if it vanished)
  static int seamSave(uint16_t t, uint16_t f, uint16_t &merged) {
    merged = 0xFFFF;
    if (t == 0xFFFF || f == 0xFFFF || gaxis(t) != gaxis(f)) return 0;
    uint16_t sum = gadd(gamt(t), gamt(f));
    merged = sum ? (uint16_t)((gaxis(t) << 12) | sum) : 0xFFFF;
    return gnz(gamt(t)) + gnz(gamt(f)) - gnz(sum);
  }
  const PairHeur *heur = nullptr;
  const FinishTable *finish = nullptr;  // built with this beam's Z
  // finish-table steering: a child with <= finishLBMis misplaced pieces that no
  // single algorithm finishes needs at least two, so its score is capped at
  // "solved minus two cheapest algorithms" (0: off); a child one algorithm
  // from solved is recorded as a solution and not kept in the beam (finishSkip)
  int finishLBMis = 0;
  bool finishSkip = false;
  int mulOverride = 0;    // quarter-moves per correct piece in the base score (0: 4 * lambda)
  // optional base score by correct count (quarter-moves, increasing in c):
  // replaces mul * c, i.e. a non-linear estimate of the remaining cost
  std::vector<int> scoreTab;
  // optional base score by correct count per orbit (cA * 25 + cB), quarter-moves;
  // overrides scoreTab.  Pair algorithms fix both orbits at once, so a state
  // with one orbit nearly solved and the other not is much dearer than its
  // total count suggests.
  std::vector<int> scoreTab2;
  // optional per-face term: faceExtra[k] (quarter-moves, >= 0) is subtracted
  // for every (orbit, face) with k of its 4 pieces correct after the algorithm.
  // A piece of colour c that lands correctly lands on face c of its orbit, so
  // k = popcount(correct-after & M_c & orbit): no child state is built.
  std::vector<int> faceExtra;
  uint64_t slotSeed = 0;  // varies which children collide in the table (diversity across runs)
  struct Slot { int score = INT32_MIN; int parent; uint32_t alg; int g4; uint64_t h; uint16_t tail; };

  // seeds (testbed): start the beam from these algorithm sets (each applied to the
  // start in order) instead of from the start; startList: levels 0 .. startListLevels-1
  // score only these algorithms
  std::vector<AlgRef> solve(const State &start, bool &ok, int T, const std::vector<std::vector<uint32_t>> *seeds = nullptr,
                            const std::vector<uint32_t> *startList = nullptr, int startListLevels = 0) const {
    ok = false;
    const Pool &pool = *cs->pool;
    const int mul = mulOverride ? mulOverride : 4 * cs->lambda;
    int tab[49];
    for (int c = 0; c <= 48; c++) tab[c] = scoreTab.empty() ? mul * c : scoreTab[c];
    const bool two = !scoreTab2.empty();
    const int *tab2 = two ? scoreTab2.data() : nullptr;
    constexpr uint64_t M24 = (1ULL << 24) - 1;
    // base score of the correct-after mask x
    auto baseScore = [&](uint64_t x) {
      return two ? tab2[__builtin_popcountll(x & M24) * 25 + __builtin_popcountll((x >> 24) & M24)]
                 : tab[__builtin_popcountll(x & M48)];
    };
    // exact table size T: slot = (hash * T) >> 64 (multiply-shift range reduction)
    const int TS = std::max(1, T);
    int curTS = TS;  // this level's table size (firstWidth experiment: the start's level differs)
    const uint64_t seedMix = slotSeed * 0x9E3779B97F4A7C15ULL;
    auto tieKey = [seedMix](uint64_t h) {
      uint64_t x = (h + seedMix) * 0xD6E8FEB86659FD93ULL;
      return x ^ (x >> 32);
    };
    auto slotOf = [&curTS, seedMix](uint64_t h) {
      uint64_t x = h ^ seedMix;
      x ^= x >> 31;
      x *= 0xBF58476D1CE4E5B9ULL;
      x ^= x >> 29;
      return (size_t)(((unsigned __int128)x * (uint64_t)curTS) >> 64);
    };
    std::vector<Node> nodes;
    nodes.push_back({start, zob(start), 0, -1, 0, 0xFFFF});
    std::vector<int> level = {0};
    int SEEN = 1 << 12;  // duplicate filter, sized to the beam
    while (SEEN < 64 * TS && SEEN < (1 << 20)) SEEN <<= 1;
    std::vector<uint64_t> seen(SEEN, 0);
    seen[nodes[0].h & (SEEN - 1)] = nodes[0].h;
    if (seeds && !seeds->empty()) {
      // each seed a chain of nodes from the start; the chains' last nodes (contiguous,
      // as a level must be) form level 0
      std::vector<int> lastInner;
      for (auto &sd : *seeds) {
        int par = 0;
        for (size_t i = 0; i + 1 < sd.size(); i++) {
          const Node &pn = nodes[par];
          State ns = pn.s.apply(pool, sd[i]);
          nodes.push_back({ns, zob(ns), pn.g4 + cs->eff4[sd[i]], par, sd[i], 0xFFFF});
          par = nodes.size() - 1;
        }
        lastInner.push_back(par);
      }
      level.clear();
      for (size_t k = 0; k < seeds->size(); k++) {
        const auto &sd = (*seeds)[k];
        if (sd.empty()) continue;
        const int par = lastInner[k];
        State ns = nodes[par].s.apply(pool, sd.back());
        nodes.push_back({ns, zob(ns), nodes[par].g4 + cs->eff4[sd.back()], par, sd.back(), 0xFFFF});
        level.push_back(nodes.size() - 1);
        seen[nodes.back().h & (SEEN - 1)] = nodes.back().h;
      }
      if (level.empty()) level = {0};
    }
    std::vector<Slot> table(TS);
    MinTree minTree;
    int best4 = 1 << 30, bestNode = -1;
    long bestAlg = -1;
    const FinishTable::Slot *bestFin = nullptr;
    const int minE4 = finish ? *std::min_element(cs->eff4.begin(), cs->eff4.end()) : 0;  // finish-table solution: bestNode, then bestAlg, then the table's algorithm
    EG2::Sol bestSol;
    const size_t n = pool.size();
    const int *e4 = cs->eff4.data();
    std::vector<uint8_t> hits(useIndex ? n : 0, 0);
    // per node rejected set (nw words) and correct-position mask, for incremental
    // updates; only the current level and its parent level are kept (each level
    // is a contiguous range of nodes)
    std::vector<uint64_t> rejArena, rejPrev, rejCorr, rejCorrPrev;
    int levelBase = 0, prevBase = 0;
    std::vector<uint32_t> touched;
    std::vector<int> bonus4, zeroBonus(1, 0);
    std::vector<std::pair<int, uint32_t>> rootScores;  // rootTop experiment: scores at the start
    std::vector<uint32_t> rootList;                     // ... and the best rootTop, in pool order
    // survivor lists: a node whose keep-solved survivors are few keeps them for
    // its children (a child's survivors are a subset of its parent's whenever no
    // correct position became incorrect); CSR per level, current and previous
    constexpr uint32_t NOLIST = UINT32_MAX;
    std::vector<uint32_t> survDat, survDatPrev, survBeg, survBegPrev, survEnd, survEndPrev;
    std::vector<uint8_t> rejOk, rejOkPrev;  // the node's incremental reject set was built

    long nAdmit = 0, nLevels = 0;
    for (int depth = 0; depth < 60 && !level.empty(); depth++) {
      const long tLevel = timing ? nsNow() : 0;
      if (stop && stop->load(std::memory_order_relaxed)) return {};
      nLevels++;
      rejArena.swap(rejPrev);
      rejCorr.swap(rejCorrPrev);
      survDat.swap(survDatPrev); survBeg.swap(survBegPrev); survEnd.swap(survEndPrev); rejOk.swap(rejOkPrev);
      survDat.clear();
      survBeg.assign(level.size(), NOLIST);
      survEnd.assign(level.size(), 0);
      rejOk.assign(level.size(), 0);
      // lists pay off below about two bitset walks, and are capped at 64MB per level
      const size_t listMax = noList ? 0 : std::min<size_t>(2 * nw, (16u << 20) / level.size());
      prevBase = levelBase;
      levelBase = level[0];
      // the incremental sets cost width * pool/8 bytes per level: above 64MB,
      // rebuild each node's set from scratch instead
      const bool incr = !noIncr && keepSolved && !useIndex && level.size() * nw * 8 <= (64u << 20);
      const bool parentIncr = incr && !rejPrev.empty();
      rejArena.resize(incr ? level.size() * nw : 0);
      rejCorr.resize(level.size());  // also used by the survivor lists
      if (!incr) std::vector<uint64_t>().swap(rejArena);
      std::vector<uint64_t> rejScratch(incr ? 0 : nw);
      curTS = depth == 0 && firstWidth > 0 ? firstWidth
              : widthGrowth != 1 ? (int)std::min<double>(widthCap, std::max(1.0, TS * std::pow(widthGrowth, depth))) : TS;
      table.resize(curTS);
      for (auto &sl : table) sl.score = INT32_MIN;
      minTree.reset(curTS);
      int filled = 0, weakest = INT32_MIN;  // weakest score in the table once every slot is filled
      size_t nodeCount = 0;
      for (int ni : level) {
        if ((++nodeCount & 63) == 0 && stop && stop->load(std::memory_order_relaxed)) return {};  // Ctrl-C within a level
        const Node &nd = nodes[ni];
        int m = 48 - nd.s.correct();
        if (m == 0) {
          if (nd.g4 < best4) { best4 = nd.g4; bestNode = ni; bestSol = EG2::Sol(); bestAlg = -1; bestFin = nullptr; }
          continue;
        }
        const long tPre = timing ? nsNow() : 0;
        if (cs->eg && m <= cs->maxMis) {
          auto sol = cs->eg->solve(nd.s, 1, 200);
          if (sol.cost < 200) {
            int e = 0;
            for (auto x : sol.path) e += cs->egEff4(cs->eg->cand, x);
            if (sol.finAlg >= 0) e += cs->egEff4(cs->eg->fin, sol.finAlg);
            if (nd.g4 + e < best4) { best4 = nd.g4 + e; bestNode = ni; bestSol = sol; bestAlg = -1; bestFin = nullptr; }
          }
        }
        if (finish) {  // the node itself one algorithm from solved
          const FinishTable::Slot *fi = finish->find(nd.h);
          if (fi && nd.g4 + fi->e4 < best4) {
            best4 = nd.g4 + fi->e4;
            bestNode = nd.parent >= 0 ? nd.parent : ni;
            bestAlg = nd.parent >= 0 ? (long)nd.alg : -2; bestFin = fi; bestSol = EG2::Sol();
          }
        }
        const long tSeam0 = timing ? nsNow() : 0;
        if (timing) nsNodePre += tSeam0 - tPre;
        if (seamAware) {
          bonus4.resize(heads.size());
          uint16_t mg;
          for (size_t hI = 0; hI < heads.size(); hI++) bonus4[hI] = 4 * seamSave(nd.tail, heads[hI], mg);
        }
        if (timing) nsSeam += nsNow() - tSeam0;
        const int *bon = seamAware ? bonus4.data() : zeroBonus.data();
        const uint16_t *hix = seamAware ? headIdx.data() : zeroIdx.data();
        const uint64_t M1 = nd.s.M[1], M2 = nd.s.M[2], M3 = nd.s.M[3], M4 = nd.s.M[4], M5 = nd.s.M[5];
        // colour bit-planes of the state
        const uint64_t C0 = M1 | M3 | M5, C1 = M2 | M3, C2 = M4 | M5;
        const uint64_t *t0 = T0.data(), *t1 = T1.data(), *t2 = T2.data();
        auto admit = [&](size_t a, int sc) {
          nAdmit++;
          int g4 = nd.g4 + e4[a] - bon[hix[a]];
          uint16_t ntail = 0xFFFF;
          if (seamAware) {
            uint16_t mg;
            int sv = seamSave(nd.tail, pool.head[a], mg);
            ntail = nGroups[a] == 1 ? (sv || (nd.tail != 0xFFFF && gaxis(nd.tail) == gaxis(pool.head[a])) ? mg : lastG[a]) : lastG[a];
          }
          if (g4 >= best4) return;
          uint64_t h = nd.h;
          const auto &src = pool.src[a];
          for (auto d : supp[a]) h ^= Z[d][nd.s.col[d]] ^ Z[d][nd.s.col[src[d]]];
          if (seen[h & (SEEN - 1)] == h) return;
          if (finish) {
            const long tF = timing ? nsNow() : 0;
            const FinishTable::Slot *fi = finish->find(h);
            if (timing) { nsAdmitFin += nsNow() - tF; nAdmitFin++; }
            if (fi) {
              if (g4 + fi->e4 < best4) {  // one more algorithm solves it
                best4 = g4 + fi->e4;
                bestNode = ni; bestAlg = a; bestFin = fi; bestSol = EG2::Sol();
              }
              if (finishSkip) return;
            } else if (finishLBMis) {
              const int ca = __builtin_popcountll(~((C0 ^ t0[a]) | (C1 ^ t1[a]) | (C2 ^ t2[a])) & M48);
              if (ca < 48 && 48 - ca <= finishLBMis) {
                sc = std::min(sc, tab[48] - g4 - 2 * minE4);
                if (sc < weakest) return;
              }
            }
          }
          if (!faceExtra.empty()) {
            const uint64_t x = ~((C0 ^ t0[a]) | (C1 ^ t1[a]) | (C2 ^ t2[a]));
            int ex = 0;
            for (int c = 0; c < 6; c++) {
              const uint64_t xc = x & nd.s.M[c];
              ex += faceExtra[__builtin_popcountll(xc & 0xFFFFFFULL)] + faceExtra[__builtin_popcountll(xc & (0xFFFFFFULL << 24))];
            }
            sc -= ex;
            if (sc < weakest) return;
          }
          if (heur) {
            sc -= heur->extra4(nd.s.apply(pool, a));
            if (sc < weakest) return;
          }
          const size_t si = slotOf(h);
          Slot &sl = table[si];
          if (sl.score == INT32_MIN) filled++;
          else if (sc < sl.score) return;
          else if (sc == sl.score && tieKey(h) >= tieKey(sl.h)) return;  // seeded tie-break: runs differ even at width 1
          sl = {sc, ni, (uint32_t)a, g4, h, ntail};
          minTree.set(si, sc);
          if (filled == curTS) {  // admission threshold: weakest slot
            if (timing) nWeakScan++;
            weakest = minTree.min();
          }
        };
        if (useIndex) {
          // Only algorithms containing a beneficial movement (a piece of colour
          // target(d) moving onto a wrong d) can gain.  Count such movements
          // per algorithm via the movement index; that count bounds the gain.
          const int c0 = 48 - m;
          touched.clear();
          for (int d = 0; d < 48; d++) {
            if (nd.s.col[d] == target(d)) continue;
            const int lo = d < 24 ? 0 : 24, t = target(d);
            for (int sp = lo; sp < lo + 24; sp++) {
              if (sp == d || nd.s.col[sp] != t) continue;
              for (uint32_t k = moff[d * 48 + sp], e = moff[d * 48 + sp + 1]; k < e; k++) {
                uint32_t a = mdat[k];
                if (hits[a]++ == 0) touched.push_back(a);
              }
            }
          }
          for (uint32_t a : touched) {
            int h = hits[a];
            hits[a] = 0;
            if (mul * (c0 + h) - (nd.g4 + e4[a]) < weakest) continue;  // gain < hits (ties reach the seeded tie-break)
            int c = __builtin_popcountll(~((C0 ^ t0[a]) | (C1 ^ t1[a]) | (C2 ^ t2[a])) & M48);
            int sc = tab[c] - (nd.g4 + e4[a]);
            if (sc >= weakest) admit(a, sc);
          }
        } else if (keepSolved) {
          // pass 1: reject (one AND) every algorithm that moves a solved piece off
          // its face, compacting the survivors; pass 2: score only those
          uint64_t corr = 0;
          for (int p = 0; p < 48; p++) corr |= (uint64_t)(nd.s.col[p] == target(p)) << p;
          // rejected = OR of the per-position bitsets over correct positions;
          // incrementally from the parent's set when no position became incorrect
          long tA = timing ? nsNow() : 0;
          const int li = ni - levelBase;
          int pi = nd.parent - prevBase;  // parent's index within the previous level
          const std::vector<uint32_t> *rList = startList && depth < startListLevels ? startList
                                               : depth >= 1 && depth <= rootLevels && !rootList.empty() ? &rootList : nullptr;
          const bool restricted = rList != nullptr;
          // a child's survivors are its parent's minus those now carrying a correct piece
          // off its face, unless some position lost its correct piece (shuffled within
          // its face): then rebuild
          const bool lost = depth > 0 && nd.parent >= 0 && (rejCorrPrev[pi] & ~corr);
          const bool fromList = !restricted && !lost && nd.parent >= 0 && !survBegPrev.empty() && survBegPrev[pi] != NOLIST;
          rejCorr[li] = corr;
          uint64_t *rej = incr ? &rejArena[(size_t)li * nw] : rejScratch.data();
          if (timing) (fromList ? nFromList : restricted ? nRestricted : nBitset)++;
          if (!fromList) {
            uint64_t base = 0;
            if (parentIncr && nd.parent >= 0 && rejOkPrev[pi] && !lost) {
              memcpy(rej, &rejPrev[(size_t)pi * nw], nw * 8);
              base = rejCorrPrev[pi];
            } else {
              if (timing && depth > 0) nRebuild++;
              memset(rej, 0, nw * 8);
              if (n % 64) rej[nw - 1] = ~0ULL << (n % 64);  // padding never survives
            }
            for (uint64_t cb = corr & ~base; cb; cb &= cb - 1) {
              const uint64_t *b = &BP[__builtin_ctzll(cb) * nw];
              for (size_t w = 0; w < nw; w++) rej[w] |= b[w];
            }
            if (incr) rejOk[li] = 1;
          }
          long tB = timing ? nsNow() : 0;
          const long admit0 = nAdmit;
          size_t k = 0;
          auto scoreOne = [&](size_t a) {
            k++;
            int c, e;
            if (!noAoS) { const Rec &R = rec[a]; c = baseScore(~((C0 ^ R.t0) | (C1 ^ R.t1) | (C2 ^ R.t2))); e = R.e4; }
            else { c = baseScore(~((C0 ^ t0[a]) | (C1 ^ t1[a]) | (C2 ^ t2[a]))); e = e4[a]; }
            int sc = c - (nd.g4 + e) + bon[hix[a]];
            if (depth == 0 && rootTop) rootScores.push_back({sc, (uint32_t)a});
            if (sc >= weakest) admit(a, sc);
          };
          if (restricted) {
            for (uint32_t a : *rList)
              if (!(rej[a >> 6] >> (a & 63) & 1)) scoreOne(a);
          } else if (fromList) {
            const uint64_t *Xa = X.data();
            const size_t b0 = survDat.size();
            for (uint32_t j = survBegPrev[pi], je = survEndPrev[pi]; j < je; j++) {
              const uint32_t a = survDatPrev[j];
              if (Xa[a] & corr) continue;
              survDat.push_back(a);
              scoreOne(a);
            }
            survBeg[li] = b0; survEnd[li] = survDat.size();
          } else {
            const size_t b0 = survDat.size();
            bool keep = listMax > 0;
            for (size_t w = 0; w < nw; w++)
              for (uint64_t bits = ~rej[w]; bits; bits &= bits - 1) {
                const size_t a = w * 64 + __builtin_ctzll(bits);
                if (keep) {
                  if (survDat.size() - b0 < listMax) survDat.push_back(a);
                  else { keep = false; survDat.resize(b0); }
                }
                scoreOne(a);
              }
            if (keep) { survBeg[li] = b0; survEnd[li] = survDat.size(); }
          }
          if (timing) {
            long tC = nsNow();
            nsFilter += tB - tA; nsScore += tC - tB;
            histNodes[48 - m]++; histSurv[48 - m] += k; histAdmit[48 - m] += nAdmit - admit0;
            lvNodes[depth]++; lvSurv[depth] += k; lvCorrect[depth] += 48 - m;
          }
          statScanned += n;
          statSurvived += k;
        } else {
        // score blocks of algorithms branch-free (vectorisable), then admit only
        // blocks whose best can beat the weakest slot.  The S_c are disjoint
        // (each source feeds one destination of one target colour), as are
        // the M_c, so the six popcounts collapse into one.
        constexpr size_t B = 64;
        int sc[B];
        for (size_t base = 0; base < n; base += B) {
          size_t mB = std::min(B, n - base);
          int mx = INT32_MIN;
          for (size_t j = 0; j < mB; j++) {
            size_t a = base + j;
            sc[j] = baseScore(~((C0 ^ t0[a]) | (C1 ^ t1[a]) | (C2 ^ t2[a]))) - (nd.g4 + e4[a]) + bon[hix[a]];
            mx = std::max(mx, sc[j]);
          }
          if (mx < weakest) continue;
          for (size_t j = 0; j < mB; j++) {
            if (sc[j] < weakest) continue;
            size_t a = base + j;
            admit(a, sc[j]);
          }
        }
      }
        }
      const long tExp = timing ? nsNow() : 0;
      level.clear();
      for (auto &sl : table) {
        if (sl.score == INT32_MIN) continue;
        State ns = nodes[sl.parent].s.apply(pool, sl.alg);
        seen[sl.h & (SEEN - 1)] = sl.h;
        nodes.push_back({ns, sl.h, sl.g4, sl.parent, sl.alg, sl.tail});
        level.push_back(nodes.size() - 1);
      }
      if (timing) nsExpand += nsNow() - tExp;
      if (depth == 0 && rootTop && rootScores.size() > 0) {
        const size_t M = std::min<size_t>(rootTop, rootScores.size());
        std::nth_element(rootScores.begin(), rootScores.begin() + (M - 1), rootScores.end(),
                         [](auto &x, auto &y) { return x.first > y.first; });
        for (size_t i = 0; i < M; i++) rootList.push_back(rootScores[i].second);
        std::sort(rootList.begin(), rootList.end());
        std::vector<std::pair<int, uint32_t>>().swap(rootScores);
      }
      if (timing) lvNs[depth] += nsNow() - tLevel;
    }
    statAdmit += nAdmit;
    statLevels += nLevels;
    statSolves++;
    std::vector<AlgRef> r;
    if (bestNode < 0) return r;
    for (int nn = bestNode; nodes[nn].parent >= 0; nn = nodes[nn].parent) r.push_back({0, nodes[nn].alg});
    std::reverse(r.begin(), r.end());
    if (bestFin) {
      if (bestAlg >= 0) r.push_back({0, (uint32_t)bestAlg});
      r.push_back({bestFin->pid, bestFin->alg});
    }
    for (auto x : bestSol.path) r.push_back({1, x});
    if (bestSol.finAlg >= 0) r.push_back({2, (uint32_t)bestSol.finAlg});
    ok = true;
    return r;
  }
};
