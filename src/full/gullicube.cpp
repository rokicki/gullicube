// Full NxNxN solve built on brobdicube's xcube:
//
//   0. corners (+ middle edges and centre centres on odd cubes): brobdicube
//   1. wings, with centres free (wings.h): one beam per wing depth over
//      face-turn setups and [slice, face] commutators, finished with
//      brobdicube's wing 3-cycles; corners and middle edges are preserved
//   2. oblique centre pairs (x,y)/(y,x): pair beam solver (merged, or
//      concatenated with --fast)
//   3. diagonal orbits (a,a) and mid-centre orbits (Mx,b): single-orbit beams
//
//   gullicube [-i | -f | -r n | -R] [-F] [-o] [-b WIDTH] [-2] [--fast] [-E] [-t THREADS]
//             [--seed S] [--setupcost Q] [--verify] N
//
// Options follow brobdicube's where they overlap.  -i reads scramble moves from
// stdin, -f a Kociemba facelet string (facelets.h), -r n / -R scramble with
// n / 50*N random moves; without any of these the state is uniformly random
// (brobdicube 'R' mode).  -F prints the start state as a facelet string.
// --verify replays the scramble (or loads the facelets) and the solution on
// an independent sticker simulator.
//
// -2: solve the whole cube at beam width 1, 2, 4, 8, ... (growth --grow, --reps
// runs per width, each with its own seed) from the same scrambled cube,
// printing each result (width, moves split into edges / corners+middle edges /
// centres, time), and keep the best.  With -b W the growth stops at width W.
// Ctrl-C stops the run (the width-1 solve is always completed) and the best
// solution is reported / written.  --perpair is the older variant (wings once,
// best per pair combined).
//
// Structure: everything loaded once (pools, solvers, tables) lives in Res and is
// read-only during solves; solveOnce() works on a copy of the scrambled cube
// and keeps all of its state local, so solves can be repeated.
#include <deque>
#include "layout.h"
#include "merge.h"
#include "solvers.h"
#include "poolcache.h"
#include "classpool.h"
#include "wings.h"
#include "replay.h"
#include "facelets.h"
#include "cornersolver.h"
#include "../pair/cube.h"
#include <atomic>
#include <chrono>
#include <csignal>
#include <map>
#include <memory>
#include <thread>


extern std::vector<mv> *g_prep_recorder;

static double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// move list of a pool algorithm with slice placeholders '2'/'3' mapped
// Merge with the beam or greedily: the beam's cost grows about as
// instances^2 x width while its gain shrinks for large merges (measured: -3%
// at 10^3 and 16^3 with width 256; -1% at 64^3 at 30x the time).
static bool gPairMerge = false;  // --pairmerge: plain pairwise merge
static bool gOldMerge = false;   // --oldmerge: merge beam / greedy merge (the previous default)
static int gMergeRestarts = 256; // --mergerestarts
static int gMergeThreads = 1;    // the merge's restarts run on the solve's threads (-t)
// restarts: at most this many restarts of the default merge (-1: --mergerestarts)
static std::vector<mv> mergeWith(Merger &mg, int beamOpt, long &raw, long &merged, std::vector<int> *order = nullptr,
                                 double budget = 0, const char *what = "", int restarts = -1) {
  size_t n = mg.insts.size();
  if (const char *dump = getenv("MERGEDUMP")) {  // experiment: write the merge input (devtools/mergebound)
    static std::atomic<int> seq{0};
    std::string fn = std::string(dump) + "-" + std::to_string(seq++) + "-" + what + ".txt";
    if (FILE *f = fopen(fn.c_str(), "w")) {
      fprintf(f, "# merge %s N %d instances %zu\n", what, mg.N, n);
      for (auto &in : mg.insts) {
        fprintf(f, "%zu", in.deps.size());
        for (auto &d : in.deps) fprintf(f, " %d %llx", d.first, (unsigned long long)d.second);
        fprintf(f, " %zu", in.moves.size());
        for (auto &m : in.moves) fprintf(f, " %d %d %d", m.axis, m.pos, m.amt);
        fprintf(f, "\n");
      }
      fclose(f);
    }
  }
  mg.budget = budget;
  if (gPairMerge) {
    mg.budget = budget;
    double t0 = now();
    auto r = mg.mergePairwise(raw, merged, order);
    if (getenv("MERGETIME")) fprintf(stderr, "  merge %s: %zu instances, pairwise: %.3fs\n", what, n, now() - t0);
    return r;
  }
  if (!gOldMerge && beamOpt < 0) {  // pairwise (or levels) + improvement + restarts within the budget
    mg.budget = budget;
    mg.restarts = restarts >= 0 ? std::min(restarts, gMergeRestarts) : gMergeRestarts;
    mg.threads = gMergeThreads;
    double t0 = now();
    auto r = mg.mergeBest(raw, merged, order);
    if (getenv("MERGETIME"))
      fprintf(stderr, "  merge %s: %zu instances, best of %ld restarts: %.3fs\n", what, n, mg.statRestarts + 1, now() - t0);
    return r;
  }
  int w = beamOpt >= 0 ? beamOpt : n <= 400 ? 256 : n <= 1500 ? 64 : 0;
  // the merge beam costs about n^2 x width x 9e-8 s (measured): fit the width to the budget
  if (beamOpt < 0 && budget > 0)
    while (w >= 4 && 9e-8 * (double)n * n * w > budget) w /= 2;
  if (w < 4) w = 0;
  double t0 = now();
  auto r = w ? mg.mergeBeam(w, raw, merged, order) : mg.merge(raw, merged, order);
  if (getenv("MERGETIME")) fprintf(stderr, "  merge: %zu instances, width %d: %.2fs\n", n, w, now() - t0);
  return r;
}

// same-axis groups -> moves with the placeholder depths substituted
static std::vector<mv> groupsMapped(const std::vector<uint16_t> &g, int d2, int d3) {
  static const int FACEOF[3][2] = {{0, 5}, {1, 3}, {2, 4}};
  std::vector<mv> r;
  for (auto k : g)
    for (int slot = 0; slot < 6; slot++) {
      int amt = (gamt(k) >> (2 * slot)) & 3, d = slot % 3 + 1;
      if (amt) r.push_back({(short)(d == 1 ? 1 : d == 2 ? d2 : d3), (short)FACEOF[gaxis(k)][slot / 3], (short)amt});
    }
  return r;
}
static std::vector<mv> mapped(const Pool &p, uint32_t a, int d2, int d3) {
  std::vector<mv> r;
  for (auto &m : poolMoves(p, a)) r.push_back({m.depth == 1 ? 1 : (m.depth == 2 ? d2 : d3), (short)m.face, (short)m.amt});
  return r;
}

// ------------------------------------------------------------------ options
struct Opts {
  int N = 0, width = 16, threads = (int)std::thread::hardware_concurrency();
  long scrLen = 0;
  bool readFacelets = false, printFacelets = false;
  bool readScramble = false, writeMoves = false, useEG = false, verify = false, fast = false, doubling = false;
  bool fallbackOnly = false, useTopK = false, useIndex = false, keepSolved = true, strictKeep = false;
  uint64_t seed = 0;  // 0: pick one from entropy (printed, so the run can be repeated)
  const char *usageFile = nullptr, *pairPoolFile = nullptr, *diagPoolFile = nullptr, *midPoolFile = nullptr;
  const char *diagFinFile = nullptr, *midFinFile = nullptr, *wingPoolFile = nullptr, *wingFinFile = nullptr;
  // data files (*.cls): data/ beside the executable as named on the command line
  // (argv[0]); ./data when that has no directory part.  No environment lookups.
  std::string dataDir = "data/";
  std::string data(const char *name) const { return dataDir + name; }
  std::string wingDefault;  // storage for the default wingPoolFile
  Weights wt;
  double grow = 2.0;   // -2: width growth factor per step
  int reps = 1;        // -2: runs per width, each with a different table seed
  int wingWidth = 16;  // -2: beam width for the (once-only) wing phase
  int wingFactor = 4;  // wing beam width = this x the -b width
  int wingBeamCost = 7;  // wing class file: beam pool = images up to this many moves
  bool oneMerge = true;  // pairs + diagonals + mids in one merge (quality path)
  int mergeBeam = -1;    // merge beam width (0: greedy merge; -1: by merge size, see mergeWith)
  bool pairMerge = false;  // pairwise (greedy-edge) merge instead (merge.h mergePairwise)
  double mergeShare = 0.5;  // greedy merge time budget = share x pair beam time
  bool seam = true;      // pair/diag/mid beams: seam-aware costs (measured: -1% quality, -4.8% --fast)
  bool whole = true;  // -2: re-solve the whole cube (wings too) at every width, keep the best whole solve
                     // (false = --perpair: wings once, best solution per pair combined across runs)
  bool finish = true;  // pair finish table (default except --fast; --finish / --nofinish)
  bool finishGiven = false;
};

// ------------------------------------------------------------------ loaded once, read-only while solving
struct Res {
  Layout L;
  Pool pairPool, fbPool, diagPool, midPool, wingPool, pair3Pool, diag3Pool, mid3Pool;
  EG2 eg;
  ColourSolver pairSolver, fbSolver, diagSolver, midSolver;
  TableBeam pairBeam, diagBeam, midBeam;
  CycleFinisher pairFin, diagFin, midFin;
  WingGens wingGens;
  WingFinisher wingFin;
  WingPool minedWings, minedWingsFin;  // --wingpool / --wingfin: mined wing algorithms (beam / finish table)
  ClassPool wingClasses;
  OrbitModel wingModelR;
  WingTableBeam wingTable;
  // pair finish table: the pair pool (AlgRef pool 0) and the endgame pools
  // streamed from their binaries (pool ids 10+i, read on demand); algorithms
  // used from those are copied into extPool (pool id 6) after each solve
  FinishTable finish, finishExt;
  std::vector<ClassPool> pairClasses;  // pair finish algorithms (symmetry classes), pool ids 10..12
  OrbitModel pairModelR;
  bool fbSame = false;  // the fallback solver uses pairPool
  mutable Pool extPool;
  mutable std::map<std::pair<int, uint32_t>, uint32_t> extIndex;
  // diagonal finish table (same scheme; external algorithms copied into diagExt, pool id 6)
  FinishTable diagFinish, diagFinishExt;
  ClassPool diagClasses;  // diagonal finish algorithms: symmetry classes, expanded on load
  OrbitModel diagModelR;
  mutable Pool diagExt;
  mutable std::map<std::pair<int, uint32_t>, uint32_t> diagExtIndex;
  FinishTable midFinish, midFinishExt;  // mid orbits (odd N): the same, into midExt
  ClassPool midClasses;
  OrbitModel midModelR;
  mutable Pool midExt;
  mutable std::map<std::pair<int, uint32_t>, uint32_t> midExtIndex;
  bool load(const Opts &o) {
    double tl = now(), tl0 = tl;
    std::deque<std::string> names;  // default data file paths (stable storage for the const char * below)
    auto D = [&](const char *name) { names.push_back(o.data(name)); return names.back().c_str(); };
    auto mark = [&](const char *what) { if (getenv("LOADTIME")) fprintf(stderr, "  load %-28s %.2fs\n", what, now() - tl); tl = now(); };
    bool odd = o.N & 1;
    if (!L.build()) { fprintf(stderr, "layout mapping failed\n"); return false; }
    pairModelR = pairModel();
    diagModelR = diagModel();
    if (odd) midModelR = midModel();
    // pools: symmetry-class files (.cls, expanded here) or expanded pool files
    auto loadAny = [&](Pool &p, const char *fn, const OrbitModel &m) {
      if (endsWith(fn, ".cls")) {
        ClassPool cp;
        if (!cp.load(fn)) { fprintf(stderr, "cannot read %s\n", fn); return false; }
        poolFromClasses(p, cp, m);
      } else loadPoolCached(p, fn);
      p.sortByCost();
      return true;
    };
    const char *pairFile = o.pairPoolFile ? o.pairPoolFile : o.fast ? D("fast6.cls") : D("p8s10.cls");
    if (!loadAny(pairPool, pairFile, pairModelR)) return false;
    mark("pair pool");
    // fallback for the rare pair the main beam fails on: p8s10 (shared when it is the pair pool)
    fbSame = std::string(pairFile) == D("p8s10.cls");
    if (!fbSame && !loadAny(fbPool, D("p8s10.cls"), pairModelR)) return false;
    mark("fallback pool");
    if (!loadAny(diagPool, o.diagPoolFile ? o.diagPoolFile : D("diag2.cls"), diagModelR)) return false;
    if (odd && !loadAny(midPool, o.midPoolFile ? o.midPoolFile : D("mid2.cls"), midModelR)) return false;
    {  // brobdicube's wing 3-cycles (the old wing finisher)
      ClassPool cp;
      if (!cp.load(D("wing3.cls"))) { fprintf(stderr, "cannot read %s\n", names.back().c_str()); return false; }
      static OrbitModel wm = wingModel();
      poolFromClasses(wingPool, cp, wm);
      wingPool.sortByCost();
    }
    mark("diag/mid/wing pools");
    if (o.useEG) {
      // older endgame search (superseded by the finish tables): expands several GB
      if (!loadAny(eg.cand, D("eg8.cls"), pairModelR) || !loadAny(eg.fin, D("eg8.cls"), pairModelR) ||
          !loadAny(eg.fin, D("p9.cls"), pairModelR) || !loadAny(eg.fin, D("c12x.cls"), pairModelR))
        return false;
      eg.build();
    }
    pairSolver.build(pairPool, o.useEG ? &eg : nullptr, o.wt);
    fbSolver.build(fbSame ? pairPool : fbPool, nullptr, o.wt);
    pairBeam.strictKeep = o.strictKeep;
    pairBeam.build(pairSolver);
    pairBeam.useIndex = o.useIndex;
    pairBeam.keepSolved = o.keepSolved && !o.useIndex;
    pairBeam.timing = getenv("BEAMTIME") != nullptr;
    pairBeam.seamAware = o.seam;
    mark("pair solvers + beam");
    if (o.finish) {
      const char *files[] = {D("c12x.cls"), D("p9.cls"), D("eg8.cls")};
      pairClasses.resize(3);
      std::vector<std::pair<const ClassPool *, uint8_t>> list;
      for (int i = 0; i < 3; i++) {
        if (!pairClasses[i].load(files[i])) { fprintf(stderr, "cannot read %s\n", files[i]); return false; }
        list.push_back({&pairClasses[i], (uint8_t)(10 + i)});
      }
      buildFinishFromClassList(finish, finishExt, pairPool, pairBeam.Z, pairSolver.eff4, list, pairModelR);
      pairBeam.finish = &finish;
      mark("pair finish table");
    }
    // guaranteed fallbacks: complete 3-cycle tables
    if (!loadAny(pair3Pool, D("pair3.cls"), pairModelR) || !loadAny(diag3Pool, D("diag3.cls"), diagModelR) ||
        (odd && !loadAny(mid3Pool, D("mid3.cls"), midModelR)))
      return false;
    pairFin.build(pair3Pool);
    diagFin.build(diag3Pool);
    if (odd) midFin.build(mid3Pool);
    mark("3-cycle fallbacks");
    diagSolver.build(diagPool, nullptr, o.wt);
    diagBeam.build(diagSolver);
    diagBeam.seamAware = o.seam;
    diagBeam.timing = pairBeam.timing;
    mark("diag solver + beam");
    if (o.finish) {
      const char *fn = o.diagFinFile ? o.diagFinFile : D("diag_x10.cls");
      if (!diagClasses.load(fn)) { fprintf(stderr, "cannot read %s\n", fn); return false; }
      buildFinishFromClasses(diagFinish, diagFinishExt, diagPool, diagBeam.Z, diagSolver.eff4, diagClasses, diagModelR, 10);
      diagBeam.finish = &diagFinish;
      mark("diag finish table");
    }
    if (odd) {
      midSolver.build(midPool, nullptr, o.wt);
      midBeam.build(midSolver);
      midBeam.seamAware = o.seam;
      if (o.finish) {
        const char *fn = o.midFinFile ? o.midFinFile : D("mid_x10.cls");
        if (!midClasses.load(fn)) { fprintf(stderr, "cannot read %s\n", fn); return false; }
        buildFinishFromClasses(midFinish, midFinishExt, midPool, midBeam.Z, midSolver.eff4, midClasses, midModelR, 10);
        midBeam.finish = &midFinish;
        mark("mid finish table");
      }
    }
    wingGens.build();
    std::vector<std::vector<std::array<int, 3>>> seqs;
    for (size_t a = 0; a < wingPool.size(); a++) {
      std::vector<std::array<int, 3>> q;
      for (auto &m : poolMoves(wingPool, a)) q.push_back({m.depth > 1 ? 1 : 0, m.face, m.amt});
      seqs.push_back(q);
    }
    wingFin.buildFromPool(wingPool, seqs);
    mark("wing finisher");
    if (o.wingPoolFile) {
      std::string wp = o.wingPoolFile;
      if (wp.size() > 4 && wp.compare(wp.size() - 4, 4, ".cls") == 0) {
        // symmetry classes: beam pool = images up to wingBeamCost moves, finish table = all
        if (!wingClasses.load(o.wingPoolFile)) { fprintf(stderr, "cannot read %s\n", o.wingPoolFile); return false; }
        wingModelR = wingModel();
        loadWingPoolFromClasses(minedWings, wingClasses, wingModelR, o.wingBeamCost);
        loadWingPoolFromClasses(minedWingsFin, wingClasses, wingModelR);
        wingTable.build(minedWings, &minedWingsFin);
      } else {
        if (!minedWings.load(o.wingPoolFile)) { fprintf(stderr, "cannot read %s\n", o.wingPoolFile); return false; }
        if (o.wingFinFile && !minedWingsFin.load(o.wingFinFile)) { fprintf(stderr, "cannot read %s\n", o.wingFinFile); return false; }
        wingTable.build(minedWings, o.wingFinFile ? &minedWingsFin : nullptr);
      }
    }
    return true;
    mark("wing table beam"); (void)tl0;
  }
  // seed every beam's tie-breaks / table slots for one run (never 0: 0 disables the wing tie-break)
  void setSeed(uint64_t seed) { pairBeam.slotSeed = diagBeam.slotSeed = midBeam.slotSeed = seed | 1; }
  void setStop(const std::atomic<bool> *s) {
    pairBeam.stop = diagBeam.stop = midBeam.stop = pairSolver.stop = fbSolver.stop = diagSolver.stop = midSolver.stop = s;
  }
  const Pool &pairPoolOf(uint8_t id) const {
    return id == 0 ? pairPool : id == 1 ? eg.cand : id == 2 ? eg.fin : id == 3 ? (fbSame ? pairPool : fbPool)
           : id == 6 ? extPool : pair3Pool;
  }
};

// ------------------------------------------------------------------ one solve
struct Result {
  bool ok = false, aborted = false, solved = false;
  long corners = 0, wings = 0, pairs = 0, diagmid = 0;  // moves per phase
  long pairsRaw = 0, diagmidRaw = 0, rawKind[2] = {0, 0}, nKind[2] = {0, 0};
  long retries = 0, retryFailed = 0, sRetries = 0, sFallback = 0;
  double tCorners = 0, tWings = 0, tPairs = 0, tPairSolve = 0, tDiag = 0, tDiagSolve = 0, tTotal = 0;
  std::vector<mv> moves;  // the whole solution, only when kept (-o / --verify)
  bool combined = false;  // pairs, diagonals and mids merged together: `pairs` holds the merged total
  long total() const { return corners + wings + pairs + diagmid; }
};

// The solve is split into three phases so that -2 can freeze the edges and
// combine the best oblique-pair solutions found at different beam widths:
//   A. runEdges:      corners + wings (once); the cube after A fixes every centre state
//   B. solvePairs:    every oblique pair from the post-A cube (per width)
//   C. finishCentres: apply chosen pair solutions to a copy of the post-A cube,
//                     then solve diagonals and mids
struct EdgesOut {
  bool ok = false, aborted = false;
  xcube cube{1, 1, 'C'};  // after corners and wings
  long corners = 0, wings = 0;
  double tCorners = 0, tWings = 0;
  std::vector<mv> moves;  // corner + wing moves (only when kept)
};
struct PairSol {
  int a, b;
  std::vector<AlgRef> sol;
  long cost = 0;
  bool ok = false;
  State start{};  // the pair's state in the pair layout
};

// fast path: pieces are checked in their own orbit layouts and moves are not
// replayed on the cube (except the centre side effects later phases read)
static bool quickPath(const Opts &o) { return o.fast && !o.verify; }

template <class F>
static void parallelFor(int n, int threads, const std::atomic<bool> *stop, F fn) {
  std::atomic<int> nx{0};
  std::vector<std::thread> th;
  for (int t = 0; t < threads; t++)
    th.emplace_back([&] { for (int i; (i = nx++) < n && !(stop && stop->load());) fn(i); });
  for (auto &t : th) t.join();
}

static EdgesOut runEdges(const Res &R, const Opts &o, const xcube &initial, int wingWidth, const std::atomic<bool> *stop,
                         bool keepMoves, uint64_t seed) {
  EdgesOut E;
  const int N = o.N, My = (N - 2) / 2, Mx = N - 2 - My;
  const bool odd = N & 1;
  auto stopped = [&] { return stop && stop->load(); };
  E.cube = initial;
  xcube &cube = E.cube;
  auto emit = [&](const std::vector<mv> &s) { if (keepMoves) E.moves.insert(E.moves.end(), s.begin(), s.end()); };

  // ---------------------------------------------------------------- 0. corners
  double ts = now();
  if (odd) {
    // Orient the centre centres with middle-layer turns, judged by xcube::ccs
    // (which matches physical cubes).  threecube::centercenters rotates the
    // other way for middle-layer turns, so brobdicube's own corner phase would
    // misorient them; after this it sees them solved and leaves them alone.
    auto bad = [&] { int w = 0; for (int f = 0; f < 6; f++) w += cube.ccs[f] != f; return w; };
    for (int step = 0; step < 2 && bad(); step++) {
      int w0 = bad();
      bool done = false;
      for (short f = 0; f < 3 && !done; f++)
        for (short tw = 1; tw < 4 && !done; tw++) {
          domove({Mx + 1, f, tw}, cube, true);
          if (bad() < w0) { emit({{Mx + 1, f, tw}}); E.corners++; done = true; }
          else domove({Mx + 1, f, (short)(4 - tw)}, cube, true);
        }
    }
    if (bad()) { fprintf(stderr, "could not orient centre centres\n"); return E; }
    for (int f = 0; f < 6; f++) cube.tc.centercenters[f] = f;
  }
  {
    std::vector<mv> rec;
    g_prep_recorder = &rec;
    std::streambuf *old = std::cout.rdbuf(nullptr);  // silence brobdicube's progress output
    solvecorners(cube);
    std::cout.rdbuf(old);
    g_prep_recorder = nullptr;
    emit(rec);
    E.corners += rec.size();
  }
  E.tCorners = now() - ts;

  // ---------------------------------------------------------------- 1. wings (centres free)
  ts = now();
  std::vector<Perm24> st(My);
  for (int k = 0; k < My; k++)
    for (int j = 0; j < 24; j++) st[k][j] = cube.wingorbits[k].a[j];
  // one beam per wing depth (in parallel); each depth's sequence is checked to
  // solve its own orbit with its face turns cancelling (solveWingsChunked)
  auto seq = solveWingsChunked(R.wingGens, R.wingFin, N, My, st, wingWidth, 1, o.threads, !o.fast, stop, seed,
                               R.minedWings.size() ? &R.wingTable : nullptr);
  if (stopped()) { E.aborted = true; return E; }
  if (seq.empty() && My > 0) { fprintf(stderr, "wing solve failed\n"); return E; }
  if (quickPath(o)) {
    // cache-friendly replay of the centres only; wings are solved by
    // construction and corners / middle edges are untouched
    CentreReplay::apply(cube, seq, o.threads);
    for (auto &wo : cube.wingorbits) uclear(wo);
  } else
    for (auto &m : seq) domove(m, cube, true);
  emit(seq);
  E.wings = seq.size();
  for (int k = 0; k < My; k++)
    for (int j = 0; j < 24; j++)
      if (cube.wingorbits[k].a[j] != j) { fprintf(stderr, "wing orbit %d not solved\n", k); return E; }
  E.tWings = now() - ts;
  E.ok = true;
  return E;
}

// Solve every oblique pair of the post-edges cube at the given width.
static bool solvePairs(const Res &R, const Opts &o, const xcube &cube, int width, const std::atomic<bool> *stop,
                       std::vector<PairSol> &pj, long &retriesOut, long &fallbackOut) {
  const int My = (o.N - 2) / 2;
  pj.clear();
  for (int a = 1; a <= My; a++)
    for (int b = a + 1; b <= My; b++) pj.push_back({a, b, {}, 0, false});
  std::atomic<long> retries{0}, retryFailed{0};
  parallelFor(pj.size(), o.threads, stop, [&](int i) {
    auto &j = pj[i];
    uint8_t X0[24], X1[24];
    readOrbit(cube, j.a, j.b, X0);
    readOrbit(cube, j.b, j.a, X1);
    State s;
    for (int p = 0; p < 48; p++) s.col[p] = R.L.phi[p] < 24 ? X0[R.L.phi[p]] : X1[R.L.phi[p] - 24];
    s.setMasks();
    if (!o.fallbackOnly) j.sol = o.useTopK ? R.pairSolver.solve(s, j.ok, width) : R.pairBeam.solve(s, j.ok, width);
    if (stop && stop->load()) return;  // interrupted: abandon (reported as aborted)
    // a beam can dead-end (every successor already visited): retry wider, bigger pool
    for (int w = 4; !j.ok && !o.fallbackOnly && w <= 256; w *= 4) {
      retries++;
      j.sol = R.fbSolver.solve(s, j.ok, w);
      for (auto &r : j.sol) r.pool = 3;
    }
    if (!j.ok) {  // guaranteed fallback: 3-cycle decomposition, orbit A then B
      retryFailed++;
      State t = s;
      j.sol.clear();
      j.ok = R.pairFin.finishOrbit(t, 0, 24, j.sol, 4) && R.pairFin.finishOrbit(t, 24, 48, j.sol, 4);
      if (!j.ok) return;
    }
    j.start = s;
  });
  // copy finish-table algorithms into extPool (sequentially: it grows), then
  // check every solution in the pair layout
  for (auto &j : pj)
    for (auto &r : j.sol)
      if (r.pool >= 10) {
        auto key = std::make_pair((int)r.pool, r.idx);
        auto it = R.extIndex.find(key);
        if (it == R.extIndex.end())
          it = R.extIndex.emplace(key, appendImage(R.pairClasses[r.pool - 10], R.pairModelR, r.idx, R.extPool)).first;
        r = {6, it->second};
      }
  parallelFor(pj.size(), o.threads, nullptr, [&](int i) {
    auto &j = pj[i];
    if (!j.ok) return;
    State t = j.start;
    j.cost = 0;
    for (auto &r : j.sol) { t = t.apply(R.pairPoolOf(r.pool), r.idx); j.cost += R.pairPoolOf(r.pool).cost[r.idx]; }
    if (t.correct() != 48) j.ok = false;
  });
  retriesOut = retries;
  fallbackOut = retryFailed;
  if (stop && stop->load()) return false;
  for (auto &j : pj)
    if (!j.ok) { fprintf(stderr, "pair (%d,%d) unsolved\n", j.a, j.b); return false; }
  if (o.usageFile) {
    std::vector<long> use(R.pairPool.size(), 0);
    for (auto &j : pj) for (auto &r : j.sol) if (r.pool == 0) use[r.idx]++;
    FILE *uf = fopen(o.usageFile, "w");
    for (size_t a = 0; a < use.size(); a++) {
      fprintf(uf, "%ld %d", use[a], (int)R.pairPool.cost[a]);
      for (int d = 0; d < 48; d++) fprintf(uf, " %d", R.pairPool.src[a][d]);
      fprintf(uf, "  %s D%d\n", groups_str(std::vector<uint16_t>(R.pairPool.gbegin(a), R.pairPool.gend(a))).c_str(),
              R.pairPool.diag[a]);
    }
    fclose(uf);
  }
  return true;
}

// Apply the chosen pair solutions to a copy of the post-edges cube, then solve
// the diagonals and mids.  Returns the complete result (including phase A).
static Result finishCentres(const Res &R, const Opts &o, const EdgesOut &E, const std::vector<PairSol> &pj, int width,
                            const std::atomic<bool> *stop, bool keepMoves, double tPairSolve = 0) {
  Result res;
  const int N = o.N, My = (N - 2) / 2, Mx = N - 2 - My;
  const bool odd = N & 1, quick = quickPath(o);
  auto stopped = [&] { return stop && stop->load(); };
  double t0 = now();
  xcube cube = E.cube;
  res.corners = E.corners;
  res.wings = E.wings;
  res.tCorners = E.tCorners;
  res.tWings = E.tWings;
  if (keepMoves) res.moves = E.moves;
  auto emit = [&](const std::vector<mv> &s) { if (keepMoves) res.moves.insert(res.moves.end(), s.begin(), s.end()); };
  auto apply = [&](const std::vector<mv> &s) {
    for (auto &m : s) domove(m, cube, true);
    emit(s);
  };

  // ---------------------------------------------------------------- 2. apply the pair solutions
  double ts = now();
  // one merge (quality path): pairs are merged first as usual; after the
  // diagonals and mids are solved, everything is merged again together, with
  // the pair algorithms that disturb a diagonal chained in the first merge's
  // order (the order the diagonal states came from); the shorter result wins
  // (the combined merge re-merges every pair instance; with many pairs the few
  // diagonal instances do not repay a second full merge)
  size_t nPairInst = 0;
  for (auto &j : pj) nPairInst += j.sol.size();
  const bool one = o.oneMerge && !quick && !o.fast && nPairInst <= 20000;
  // merge time budget: a share of the pair beam's time, at least 1 s (small
  // cubes: the beam takes milliseconds, and the merge beam is worth ~3% there)
  const double mergeBudget = o.mergeShare > 0 ? std::max(1.0, o.mergeShare * tPairSolve) : 0;
  // replay the merged sequences on the cube (and check it) only when verifying:
  // replaying millions of moves costs ~N per move and dominated large quality solves
  const bool replay = !quick && (o.verify || getenv("REPLAY"));
  std::vector<int> pairOrder;        // first merge: instance order
  std::vector<std::pair<int, int>> pairInst;  // instance -> (pair, algorithm index)
  std::unique_ptr<xcube> cubeBeforePairs;  // for the combined merge
  size_t movesBeforePairs = 0;
  long pairsMerged1 = 0;
  // Pair algorithms' side effects on the diagonal orbits, applied as
  // permutations in solution order (pair by pair, algorithm by algorithm).
  // touched[i][k]: bit 0 / 1 = algorithm k of pair i moves diagonal a / b.
  std::vector<std::vector<uint8_t>> touched(pj.size());
  // ord: (pair, algorithm index) in the order the algorithms are applied
  auto applyDiagEffects = [&](const std::vector<std::pair<int, int>> &ord) {
    std::map<std::pair<int, uint32_t>, std::pair<std::array<uint8_t, 24>, std::array<uint8_t, 24>>> cache;
    for (size_t i = 0; i < pj.size(); i++) touched[i].assign(pj[i].sol.size(), 0);
    for (auto [i, k] : ord) {
      auto &j = pj[i];
      {
        const AlgRef &r = j.sol[k];
        const Pool &p = R.pairPoolOf(r.pool);
        uint8_t t = 0;
        if (p.diag[r.idx]) {
          auto key = std::make_pair((int)r.pool, r.idx);
          auto it = cache.find(key);
          if (it == cache.end()) {
            std::vector<mv> alg;
            for (auto &m : poolMoves(p, r.idx)) alg.push_back({m.depth, (short)m.face, (short)m.amt});
            std::array<uint8_t, 24> D2, D3;
            diagEffect(alg, D2, D3);
            it = cache.emplace(key, std::make_pair(D2, D3)).first;
          }
          for (int w = 0; w < 2; w++) {
            int d = w ? j.b : j.a;
            const auto &D = w ? it->second.second : it->second.first;
            bool id = true;
            for (int k = 0; k < 24; k++) id &= D[k] == k;
            if (id) continue;
            t |= 1 << w;
            uint8_t cur[24], nw[24];
            readOrbit(cube, d, d, cur);
            for (int k = 0; k < 24; k++) nw[k] = cur[D[k]];
            writeOrbit(cube, d, d, nw);
          }
        }
        touched[i][k] = t;
      }
    }
  };
  auto solutionOrder = [&] {
    std::vector<std::pair<int, int>> ord;
    for (size_t i = 0; i < pj.size(); i++)
      for (size_t k = 0; k < pj[i].sol.size(); k++) ord.push_back({(int)i, (int)k});
    return ord;
  };
    if (quick) {
    // No merging, and the pair moves are not replayed on the cube (each pair is
    // solved by construction).  Their only side effect that matters later is on
    // the diagonal orbits, applied directly as permutations in solution order.
    // Each pair's own sequence is simplified (same-axis groups combined at the
    // seams between its algorithms; linear, and exact in placeholder notation
    // since 2 and 3 are distinct layers for every pair).
    for (auto &j : pj) {
      std::vector<uint16_t> g;
      for (auto &r : j.sol) {
        const Pool &p = R.pairPoolOf(r.pool);
        res.pairsRaw += p.cost[r.idx];
        merge_into(g, std::vector<uint16_t>(p.gbegin(r.idx), p.gend(r.idx)));
      }
      res.pairs += groups_cost(g);
      if (keepMoves) emit(groupsMapped(g, j.a + 1, j.b + 1));
    }
    applyDiagEffects(solutionOrder());
  } else {
    Merger mg(N);
    for (size_t i = 0; i < pj.size(); i++)
      for (auto &r : pj[i].sol) {
        const Pool &p = R.pairPoolOf(r.pool);
        mg.add(mapped(p, r.idx, pj[i].a + 1, pj[i].b + 1), {{(int)i, p.support[r.idx]}});
      }
    std::vector<mv> seq;
    if (o.fast) {  // no merging: plain concatenation
      for (auto &in : mg.insts)
        for (auto &l : in.moves) seq.push_back(mg.toMove(l));
      res.pairsRaw = res.pairs = seq.size();
    } else {
      if (getenv("MERGESTAT")) {  // face turns in the raw pair sequences (the only moves that can cancel between algorithms)
        long faces = 0, all = 0;
        for (auto &in : mg.insts) for (auto &l : in.moves) { all++; faces += l.pos == 1 || l.pos == N; }
        fprintf(stderr, "MERGESTAT pairs: %zu algorithm instances, %ld raw moves, %ld face turns\n", mg.insts.size(), all, faces);
      }
      mg.stop = stop;
      // with the combined merge to come, this one only orders the pairs (and is a fallback): no restarts
      seq = mergeWith(mg, o.mergeBeam, res.pairsRaw, res.pairs, &pairOrder, mergeBudget, "pairs", one ? 0 : -1);
      if (stopped()) { res.aborted = true; return res; }
      pairInst = solutionOrder();  // instance index -> (pair, algorithm)
      pairsMerged1 = res.pairs;
    }
    movesBeforePairs = res.moves.size();
    if (replay) {
      if (one) cubeBeforePairs = std::make_unique<xcube>(cube);
      apply(seq);
      for (int y = 1; y <= My; y++)  // check: every oblique orbit solved
        for (int x = 1; x <= My; x++) {
          if (x == y) continue;
          uint8_t c[24];
          readOrbit(cube, x, y, c);
          for (int k = 0; k < 24; k++)
            if (c[k] != k / 4) { fprintf(stderr, "oblique orbit (%d,%d) not solved\n", x, y); return res; }
        }
    } else {
      // no replay: every pair was checked in its own layout, and the merged
      // sequence is exactly the product of its instances in merge order, so
      // the diagonals get the pair algorithms' effects in that order
      emit(seq);
      std::vector<std::pair<int, int>> ord;
      if (o.fast) ord = solutionOrder();
      else for (int q : pairOrder) ord.push_back(pairInst[q]);
      applyDiagEffects(ord);
    }
  }
  res.tPairs = now() - ts;
  if (stopped()) { res.aborted = true; return res; }

  // ---------------------------------------------------------------- 3. diagonals, mids
  ts = now();
  struct SJob { int kind, a; std::vector<AlgRef> sol; bool ok; State start{}; };  // kind 0 diag, 1 mid
  std::vector<SJob> sj;
  for (int a = 1; a <= My; a++) sj.push_back({0, a, {}, false, {}});
  if (odd)
    for (int b = 1; b <= My; b++) sj.push_back({1, b, {}, false, {}});
  auto sPool = [&](const SJob &j, const AlgRef &r) -> const Pool & {
    if (r.pool == 5) return j.kind == 0 ? R.diag3Pool : R.mid3Pool;
    if (r.pool == 6) return j.kind == 0 ? R.diagExt : R.midExt;
    return j.kind == 0 ? R.diagPool : R.midPool;
  };
  std::atomic<long> sRetries{0}, sFallback{0};
  parallelFor(sj.size(), o.threads, stop, [&](int i) {
    auto &j = sj[i];
    uint8_t c[24];
    readOrbit(cube, j.kind == 0 ? j.a : Mx, j.a, c);
    State s;
    for (int p = 0; p < 24; p++) s.col[p] = c[p];
    for (int p = 24; p < 48; p++) s.col[p] = (p % 24) / 4;
    s.setMasks();
    const ColourSolver &sv = j.kind == 0 ? R.diagSolver : R.midSolver;
    if (!o.fallbackOnly)
      j.sol = o.useTopK ? sv.solve(s, j.ok, width) : (j.kind == 0 ? R.diagBeam : R.midBeam).solve(s, j.ok, width);
    if (stopped()) return;
    for (int w = 4 * std::max(1, width); !j.ok && !o.fallbackOnly && w <= 256; w *= 4) {
      sRetries++;
      j.sol = sv.solve(s, j.ok, w);
    }
    if (!j.ok) {  // guaranteed fallback: 3-cycle decomposition (pool id 5)
      sFallback++;
      State t = s;
      j.sol.clear();
      j.ok = (j.kind == 0 ? R.diagFin : R.midFin).finishOrbit(t, 0, 24, j.sol, 5);
    }
    j.start = s;
  });
  // copy finish-table algorithms into diagExt (sequentially), then check every
  // solution in the orbit's own layout
  for (auto &j : sj)
    for (auto &r : j.sol)
      if (r.pool >= 10) {
        auto key = std::make_pair((int)r.pool, r.idx);
        auto &index = j.kind == 0 ? R.diagExtIndex : R.midExtIndex;
        auto it = index.find(key);
        if (it == index.end())
          it = index.emplace(key, j.kind == 0 ? appendImage(R.diagClasses, R.diagModelR, r.idx, R.diagExt)
                                              : appendImage(R.midClasses, R.midModelR, r.idx, R.midExt)).first;
        r = {6, it->second};
      }
  parallelFor(sj.size(), o.threads, nullptr, [&](int i) {
    auto &j = sj[i];
    if (!j.ok) return;
    State t = j.start;
    for (auto &r : j.sol) t = t.apply(sPool(j, r), r.idx);
    if (t.correct() != 48) j.ok = false;
  });
  res.sRetries = sRetries;
  res.sFallback = sFallback;
  res.tDiagSolve = now() - ts;
  if (stopped()) { res.aborted = true; return res; }
  {
    Merger mg(N);
    const int P = pj.size();
    for (size_t i = 0; i < sj.size(); i++) {
      auto &j = sj[i];
      if (!j.ok) { fprintf(stderr, "orbit kind %d #%d unsolved\n", j.kind, j.a); return res; }
      res.nKind[j.kind]++;
      for (auto &r : j.sol) {
        const Pool &p = sPool(j, r);
        res.rawKind[j.kind] += p.cost[r.idx];
        int d2 = j.a + 1, d3 = j.kind == 0 ? j.a + 1 : Mx + 1;  // diagonal: both at a+1; mid: '3' = middle layer
        mg.add(mapped(p, r.idx, d2, d3), {{(int)i, p.support[r.idx]}});
      }
    }
    if (quick) {  // each orbit was checked in its own layout; count, don't replay
      for (auto &in : mg.insts) {
        res.diagmidRaw += in.moves.size();
        if (keepMoves)
          for (auto &l : in.moves) res.moves.push_back(mg.toMove(l));
      }
      res.diagmid = res.diagmidRaw;
    } else {
      std::vector<mv> seq;
      if (o.fast) {
        for (auto &in : mg.insts)
          for (auto &l : in.moves) seq.push_back(mg.toMove(l));
        res.diagmidRaw = res.diagmid = seq.size();
      } else {
        mg.stop = stop;
        seq = mergeWith(mg, o.mergeBeam, res.diagmidRaw, res.diagmid, nullptr, one ? mergeBudget : 0.2 * mergeBudget, "diagmid",
                        one ? 0 : -1);
        if (stopped()) { res.aborted = true; return res; }
      }
      if (one) {
        // components: pair i -> i; diagonal orbit a -> P + a; mid orbit b -> P + My + 1 + b
        std::map<std::pair<int, uint32_t>, uint8_t> tcache;  // (pool, alg) -> diagonals disturbed (bit 0: a, 1: b)
        auto touches = [&](const AlgRef &r) {
          const Pool &p = R.pairPoolOf(r.pool);
          if (!p.diag[r.idx]) return (uint8_t)0;
          auto key = std::make_pair((int)r.pool, r.idx);
          auto it = tcache.find(key);
          if (it != tcache.end()) return it->second;
          std::vector<mv> alg;
          for (auto &m : poolMoves(p, r.idx)) alg.push_back({m.depth, (short)m.face, (short)m.amt});
          std::array<uint8_t, 24> D2, D3;
          diagEffect(alg, D2, D3);
          uint8_t t = 0;
          for (int k = 0; k < 24; k++) { if (D2[k] != k) t |= 1; if (D3[k] != k) t |= 2; }
          return tcache[key] = t;
        };
        Merger all(N);
        for (int q : pairOrder) {  // pair instances in the first merge's order
          auto [i, k] = pairInst[q];
          const auto &r = pj[i].sol[k];
          const Pool &p = R.pairPoolOf(r.pool);
          std::vector<std::pair<int, uint64_t>> deps = {{i, p.support[r.idx]}};
          uint8_t t = touches(r);
          if (t & 1) deps.push_back({P + pj[i].a, ~0ULL});
          if (t & 2) deps.push_back({P + pj[i].b, ~0ULL});
          all.add(mapped(p, r.idx, pj[i].a + 1, pj[i].b + 1), deps);
        }
        for (size_t i = 0; i < sj.size(); i++)
          for (auto &r : sj[i].sol) {
            const Pool &p = sPool(sj[i], r);
            int d2 = sj[i].a + 1, d3 = sj[i].kind == 0 ? sj[i].a + 1 : Mx + 1;
            int comp = sj[i].kind == 0 ? P + sj[i].a : P + My + 1 + sj[i].a;
            all.add(mapped(p, r.idx, d2, d3), {{comp, p.support[r.idx]}});
          }
        all.stop = stop;
        long raw, merged;
        auto seqAll = mergeWith(all, o.mergeBeam, raw, merged, nullptr, mergeBudget, "all");
        if (stopped()) { res.aborted = true; return res; }
        if (merged < pairsMerged1 + res.diagmid) {  // the combined merge is shorter: use it
          if (replay) cube = *cubeBeforePairs;
          if (keepMoves) res.moves.resize(movesBeforePairs);  // drop the first merge's pair moves
          res.pairs = merged;
          res.diagmid = 0;
          res.combined = true;
          seq = std::move(seqAll);
        }
      }
      if (replay) apply(seq);
      else emit(seq);
    }
  }
  res.tDiag = now() - ts;

  // ---------------------------------------------------------------- final check
  bool solved = true;
  if (replay)  // otherwise every centre orbit was checked in its own layout
    for (int y = 1; y <= My; y++)
      for (int x = 1; x <= Mx; x++) {
        uint8_t c[24];
        readOrbit(cube, x, y, c);
        for (int k = 0; k < 24; k++) solved &= c[k] == k / 4;
      }
  for (auto &w : cube.wingorbits)
    for (int k = 0; k < 24; k++) solved &= w.a[k] == k;
  for (int k = 0; k < 24; k++) solved &= cube.tc.corner.a[k] == k / 4;
  if (odd)
    for (int k = 0; k < 24; k++) solved &= cube.tc.edge.a[k] == k / 4;
  for (int f = 0; f < 6; f++) solved &= cube.ccs[f] == f;  // ccs, not tc.centercenters (see above)
  res.solved = solved;
  res.ok = true;
  res.tTotal = now() - t0;
  return res;
}

// A single complete solve at one width: A, B, C.
static Result solveOnce(const Res &R, const Opts &o, const xcube &initial, int width, const std::atomic<bool> *stop,
                        bool keepMoves, uint64_t seed) {
  double t0 = now();
  EdgesOut E = runEdges(R, o, initial, width * o.wingFactor, stop, keepMoves, seed);
  if (!E.ok) { Result r; r.aborted = E.aborted; return r; }
  double ts = now();
  std::vector<PairSol> pj;
  long retries = 0, fallback = 0;
  if (!solvePairs(R, o, E.cube, width, stop, pj, retries, fallback)) { Result r; r.aborted = stop && stop->load(); return r; }
  double tPairSolve = now() - ts;
  Result r = finishCentres(R, o, E, pj, width, stop, keepMoves, tPairSolve);
  r.retries = retries;
  r.retryFailed = fallback;
  r.tPairSolve = tPairSolve;
  r.tPairs += tPairSolve;
  r.tTotal = now() - t0;
  return r;
}

// ------------------------------------------------------------------ reporting
static void printDetail(const Opts &o, const Result &r) {
  int My = (o.N - 2) / 2;
  long nPairs = (long)My * (My - 1) / 2;
  printf("N=%d  %s  pairs %ld, diagonals %d, mids %d, wing orbits %d\n", o.N, r.solved ? "SOLVED" : "NOT SOLVED", nPairs,
         My, (o.N & 1) ? My : 0, My);
  printf("  corners        %6ld moves\n", r.corners);
  printf("  wings          %6ld moves, %.1f per wing orbit; %.2fs\n", r.wings, (double)r.wings / std::max(1, My), r.tWings);
  if (r.combined) {
    long raw = r.pairsRaw + r.diagmidRaw;
    printf("  pairs+diag+mid %6ld moves, merged together (raw: pairs %ld + diag/mid %ld = %ld, %.1f%% merged away);\n"
           "                 pair solve %.2fs, diag/mid solve %.2fs, total %.2fs\n", r.pairs, r.pairsRaw, r.diagmidRaw, raw,
           100.0 * (raw - r.pairs) / std::max(1L, raw), r.tPairSolve, r.tDiagSolve, r.tPairs + r.tDiag);
  } else {
    printf("  oblique pairs  %6ld moves (raw %ld, %.1f%% merged away), %.2f per pair; solve %.2fs, total %.2fs\n", r.pairs,
           r.pairsRaw, 100.0 * (r.pairsRaw - r.pairs) / std::max(1L, r.pairsRaw), (double)r.pairs / std::max(1L, nPairs),
           r.tPairSolve, r.tPairs);
    printf("  diag+mid       %6ld moves (raw %ld, %.1f%% merged away); solve %.2fs, total %.2fs\n", r.diagmid, r.diagmidRaw,
           100.0 * (r.diagmidRaw - r.diagmid) / std::max(1L, r.diagmidRaw), r.tDiagSolve, r.tDiag);
  }
  printf("    raw per orbit: diagonal %.1f, mid %.1f\n", (double)r.rawKind[0] / std::max(1L, r.nKind[0]),
         (double)r.rawKind[1] / std::max(1L, r.nKind[1]));
  if (r.retries || r.retryFailed)
    printf("  pair retries   %ld (wider beam, bigger pool); 3-cycle fallback used %ld\n", r.retries, r.retryFailed);
  if (r.sRetries || r.sFallback) printf("  diag/mid retries %ld; 3-cycle fallback used %ld\n", r.sRetries, r.sFallback);
}

// Conventional move names: no layer deeper than the middle.  Layer d of face f
// is layer N+1-d of the opposite face turned the other way, so a move deeper
// than N/2 (the middle slice of an odd cube is kept) is renamed that way.  The
// permutation is the same; the move count does not change.
static void conventionalMoves(std::vector<mv> &ms, int N) {
  for (auto &m : ms)
    if (2 * m.dep > N + 1) {
      m.dep = N + 1 - m.dep;
      m.face = oppface[m.face];
      m.twist = 4 - m.twist;
    }
}

static std::atomic<bool> gStop{false};  // set by Ctrl-C during -2
static void onSigint(int) {
  gStop = true;
  std::signal(SIGINT, SIG_DFL);  // a second Ctrl-C kills the process
}

static void usage() {
  fputs(
      "usage: gullicube N [options]      solve a scrambled N x N x N cube (N >= 2)\n"
      "  scramble:\n"
      "    -r n            scramble with n random moves\n"
      "    -R              scramble with 50*N random moves\n"
      "    -i              read the scramble moves from stdin\n"
      "    -f              read the cube state from stdin as a Kociemba facelet string: 6*N*N\n"
      "                    letters, faces U R F D L B, each row by row as seen in the usual net\n"
      "                    (U with B at the top, D with F at the top, the others with U at the\n"
      "                    top); each letter names the face whose colour the sticker has when\n"
      "                    solved; whitespace is ignored\n"
      "    -F              print the start state as a Kociemba facelet string\n"
      "    --seed s        random seed (default: from entropy; the seed used is printed)\n"
      "  search:\n"
      "    -b w            beam width, all phases (with -2: the maximum width)\n"
      "    -2              growing widths until Ctrl-C: the whole cube is re-solved at each width\n"
      "                    (wings at 4x the width, a fresh seed per run); the best solve is kept\n"
      "    --grow f        -2: width growth factor (default 2)\n"
      "    --reps n        -2: runs per width, each with a different table seed (default 1)\n"
      "    --perpair       -2: solve wings once (at --wingwidth) and keep the best solution per\n"
      "                    oblique pair across runs (older mode; measured worse than a flat beam)\n"
      "    --wingwidth w   -2 --perpair: wing beam width (default 16)\n"
      "    --fast          speed over moves: small pair pool, no merging, width 1 unless -b\n"
      "    -E              older pair endgame search (superseded by the finish tables; several GB)\n"
      "    --onemerge/--twomerge  one merge for pairs + diagonals + mids (default), or pairs first\n"
      "    --mergeshare f  greedy merge time budget = f x the pair beam time (default 0.5; 0 = unlimited)\n"
      "    --mergebeam w   merge with a beam of width w (0: the greedy merge); with --oldmerge the\n"
      "                    default is 256 up to 400 algorithm instances, 64 up to 1500, else greedy\n"
      "    --mergerestarts k  merge: at most k seeded restarts of the pairwise merge within the\n"
      "                    merge time budget, on the -t threads (default 256; -2: at least 4 x width)\n"
      "    --oldmerge      merge with the merge beam / greedy merge (the previous default)\n"
      "    --pairmerge     merge with the plain pairwise merge (no improvement, no restarts)\n"
      "    --seam/--noseam pair/diagonal/mid beams: cost minus moves cancelling with the previous\n"
      "                    algorithm (default on)\n"
      "    --finish        pair finish table: every beam child one algorithm from solved is seen\n"
      "                    (from c12x/p9/eg8.cls, built on load; default unless --fast)\n"
      "    --nofinish      no finish table\n"
      "    --setupcost q   quarter moves charged per setup face turn (default 3)\n"
      "    -t n            threads (default: all cores)\n"
      "  data files (*.cls) are read from data/ beside the executable as named on the\n"
      "  command line (else ./data)\n"
      "  output and checks:\n"
      "    -o              write the solution (lines start with a space)\n"
      "    --verify        replay scramble + solution on the independent sticker simulator\n"
      "    --usage file    write how often each pair-pool algorithm was used\n"
      "  testing / comparison:\n"
      "    --fallback-only every orbit via the guaranteed 3-cycle fallback\n"
      "    --nokeep        allow algorithms that carry solved pieces off their face\n"
      "    --strictkeep    never move a solved piece at all\n"
      "    --topk          older top-k beam\n"
      "    --index         movement-index scoring\n"
      "    --pairpool f    override the pair pool file\n"
      "    --diagpool f    override the diagonal pool file (default diag2.pool)\n"
      "    --midpool f     override the mid-orbit pool file (default mid2.pool)\n"
      "    --diagfin f     diagonal finish-table class file (default data/diag_x10.cls)\n"
      "    --midfin f      mid finish-table class file (default data/mid_x10.cls)\n"
      "    --wingpool f    wings: mined wing algorithms for the table beam: a symmetry-class file\n"
      "                    (default data/wing_d6c.cls; beam pool = images up to --wingbeamcost moves,\n"
      "                    finish table = all) or an expanded pool (with --wingfin for the finish\n"
      "                    table).  --fast uses the older wing beam unless --wingpool is given\n"
      "    --wingbeamcost c  wing beam pool: class images up to c moves (default 7)\n"
      "    --wingfin f     wings: finish-table pool for an expanded --wingpool\n"
      "    --oldwings      wings: the older setup/commutator beam\n"
      "    --wingfactor f  wing beam width = f x the -b width (default 4)\n",
      stderr);
}

int main(int argc, char **argv) {
  if (argc < 2) { usage(); return 1; }
  printf("#");
  for (int i = 0; i < argc; i++) printf(" %s", argv[i]);
  printf("\n");
  Opts o;
  {  // data/ beside the executable, as named on the command line; else ./data
    const std::string a0 = argv[0];
    const size_t k = a0.find_last_of("/\\");
    if (k != std::string::npos) o.dataDir = a0.substr(0, k + 1) + "data/";
    if (FILE *f = fopen(o.data("p8s10.cls").c_str(), "rb")) fclose(f);
    else o.dataDir = "data/";
    if (FILE *f = fopen(o.data("p8s10.cls").c_str(), "rb")) fclose(f);
    else fprintf(stderr, "gullicube: data files not found: run it by a path (e.g. ./gullicube, or\n"
                         "  path/to/gullicube/gullicube) or from its directory, which holds data/\n");
  }
  o.wt.setupFace = 3;  // setup face turns at 3/4 move (best with seam-aware costs)
  bool widthGiven = false, autoLen = false, oldWings = false;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a == "-i") o.readScramble = true;                 // read scramble moves from stdin
    else if (a == "-f") o.readFacelets = true;            // read a Kociemba facelet string from stdin
    else if (a == "-F") o.printFacelets = true;           // print the start state as a facelet string
    else if (a == "-r") o.scrLen = atol(argv[++i]);       // scramble with n random moves
    else if (a == "-R") autoLen = true;                   // scramble with 50*N random moves
    else if (a == "-o") o.writeMoves = true;              // write solution (lines start with a space)
    else if (a == "-t") o.threads = atoi(argv[++i]);
    else if (a == "-b") { o.width = atoi(argv[++i]); widthGiven = true; }  // beam width, all phases
    else if (a == "-2") o.doubling = true;                // growing widths until Ctrl-C; keep the best per pair
    else if (a == "--grow") o.grow = atof(argv[++i]);     // -2: width growth factor (default 2)
    else if (a == "--reps") o.reps = atoi(argv[++i]);     // -2: runs per width with different table seeds
    else if (a == "--wingwidth") o.wingWidth = atoi(argv[++i]);  // -2: wing beam width (default 16)
    else if (a == "--whole") o.whole = true;              // -2: re-solve everything at each width (default)
    else if (a == "--perpair") o.whole = false;           // -2: wings once; best per pair across runs
    else if (a == "--fast") o.fast = true;                // speed over moves: small pool, no merging
    else if (a == "-E") o.useEG = true;                   // pair endgame search
    else if (a == "--finish") { o.finish = true; o.finishGiven = true; }     // pair finish table (default unless --fast)
    else if (a == "--nofinish") { o.finish = false; o.finishGiven = true; }
    else if (a == "--seed") o.seed = atoll(argv[++i]);
    else if (a == "--setupcost") o.wt.setupFace = atoi(argv[++i]);  // quarter moves per setup face turn
    else if (a == "--verify") o.verify = true;            // replay on the independent sticker simulator
    else if (a == "--fallback-only") o.fallbackOnly = true;  // test: every orbit via the 3-cycle fallback
    else if (a == "--nokeep") o.keepSolved = false;       // allow algorithms that carry solved pieces off their face
    else if (a == "--strictkeep") o.strictKeep = true;    // never move a solved piece at all
    else if (a == "--topk") o.useTopK = true;             // older top-k beam, for comparison
    else if (a == "--index") o.useIndex = true;           // movement-index scoring, for comparison
    else if (a == "--pairpool") o.pairPoolFile = argv[++i];  // override the pair pool
    else if (a == "--diagpool") o.diagPoolFile = argv[++i];  // override the diagonal pool
    else if (a == "--midpool") o.midPoolFile = argv[++i];    // override the mid-orbit pool
    else if (a == "--diagfin") o.diagFinFile = argv[++i];    // diagonal finish-table pool
    else if (a == "--midfin") o.midFinFile = argv[++i];      // mid finish-table pool
    else if (a == "--wingpool") o.wingPoolFile = argv[++i];  // wings: table beam over this mined pool
    else if (a == "--wingfin") o.wingFinFile = argv[++i];    // wings: finish table from this pool
    else if (a == "--oldwings") oldWings = true;             // wings: the older WingBeam
    else if (a == "--wingfactor") o.wingFactor = atoi(argv[++i]);  // wing width = factor x -b
    else if (a == "--wingbeamcost") o.wingBeamCost = atoi(argv[++i]);  // wing beam pool: images up to c moves
    else if (a == "--seam") o.seam = true;                   // pair/diag/mid beams: seam-aware costs
    else if (a == "--onemerge") o.oneMerge = true;           // pairs + diagonals + mids merged together
    else if (a == "--twomerge") o.oneMerge = false;          // pairs merged, then diagonals + mids
    else if (a == "--pairmerge") o.pairMerge = gPairMerge = true;  // plain pairwise merge (comparison)
    else if (a == "--oldmerge") gOldMerge = true;                   // previous default merge (comparison)
    else if (a == "--mergerestarts") gMergeRestarts = atoi(argv[++i]);
    else if (a == "--mergebeam") o.mergeBeam = atoi(argv[++i]);  // merge beam width (0: greedy; default by size)
    else if (a == "--mergeshare") o.mergeShare = atof(argv[++i]);  // merge time budget / pair beam time
    else if (a == "--noseam") o.seam = false;
    else if (a == "--usage") o.usageFile = argv[++i];     // write how often each pair-pool algorithm was used
    else if (a[0] == '-') { fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    else if (isdigit((unsigned char)a[0])) o.N = atoi(a.c_str());
    else { fprintf(stderr, "unknown argument %s\n", a.c_str()); return 2; }
  }
  if (autoLen) o.scrLen = 50L * o.N;
  if (o.fast && !widthGiven) o.width = 1;
  gMergeThreads = std::max(1, o.threads);
  if (o.fast && !o.finishGiven) o.finish = false;
  // wings: the mined-pool table beam by default, except --fast (thousands of
  // wing orbits at huge N) or --oldwings
  if (!oldWings && !o.fast && !o.wingPoolFile) { o.wingDefault = o.data("wing_d6c.cls"); o.wingPoolFile = o.wingDefault.c_str(); }
  if (oldWings) o.wingPoolFile = o.wingFinFile = nullptr;
  if (o.N < 2) { usage(); return 1; }

  if (!o.seed) {
    std::random_device rd;
    o.seed = ((uint64_t)rd() << 32 ^ rd() ^ (uint64_t)time(nullptr) * 0x9E3779B97F4A7C15ULL ^ ((uint64_t)(uintptr_t)&o << 16)) % 1000000000000ULL;
    if (!o.seed) o.seed = 1;
  }
  printf("# seed %llu\n", (unsigned long long)o.seed);
  rng.seed(o.seed);
  const int N = o.N, My = (N - 2) / 2, Mx = N - 2 - My;

  double t0 = now();
  auto R = std::make_unique<Res>();
  if (!R->load(o)) return 1;
  double tLoad = now() - t0;

  // the scrambled cube; every solve works on a copy
  std::vector<mv> scramble;
  std::string facelets;  // with -f
  if (o.readFacelets) {
    std::string line;
    while (std::getline(std::cin, line)) facelets += line;
  } else if (o.readScramble) {
    std::string line, all;
    while (std::getline(std::cin, line)) all += line + " ";
    scramble = parse_moves(all);
  } else
    for (long i = 0; i < o.scrLen; i++)
      scramble.push_back({1 + (int)(rng() % ((N + 1) / 2)), (short)(rng() % 6), (short)(1 + rng() % 3)});
  const bool fromMoves = o.readScramble || o.scrLen;
  moveset = scramble;  // brobdicube's 'I' construction applies this global
  xcube initial(Mx, My, o.readFacelets ? 'C' : fromMoves ? 'I' : 'R');
  moveset.clear();
  if (o.readFacelets) {
    std::string err;
    if (!xcube_from_facelets(initial, facelets, err)) { fprintf(stderr, "bad facelets: %s\n", err.c_str()); return 1; }
  }
  if (o.printFacelets) printf("facelets %s\n", xcube_to_facelets(initial).c_str());

  const bool keep = o.writeMoves || o.verify;
  const int baseRestarts = gMergeRestarts;
  Result best;
  if (!o.doubling) {
    R->setSeed(o.seed);
    best = solveOnce(*R, o, initial, o.width, nullptr, keep, o.seed | 1);
    if (!best.ok) return 1;
    printDetail(o, best);
    printf("  TOTAL          %6ld moves   (load %.1fs, corners %.2fs, wall %.3fs)\n", best.total(), tLoad, best.tCorners,
           best.tTotal);
  } else {
    std::signal(SIGINT, onSigint);
    // Edges are solved once (phase A) so every width sees the same centre
    // states; each oblique pair keeps its best solution across all widths.
    // With -b, widths grow only up to that value (to bound the time).
    const int maxW = widthGiven ? o.width : INT32_MAX;
    auto runSeed = [&](long k) { uint64_t x = (o.seed + k * 0x9E3779B97F4A7C15ULL) * 0xBF58476D1CE4E5B9ULL; return (x ^ (x >> 31)) | 1; };
    if (o.whole) {
      // --whole: every run re-solves the whole cube (wings included) at the
      // step's width (wings at 4x) with its own seed; the best complete
      // solution is kept, nothing is combined across runs.
      printf("N=%d: whole-cube solves at widths 1, then x%.2f each step%s, %d run(s) per width (Ctrl-C to stop)\n", N,
             o.grow, widthGiven ? (", up to " + std::to_string(maxW)).c_str() : "", o.reps);
      long runNo = 0;
      for (int w = 1; w <= maxW;) {
        for (int rep = 0; rep < o.reps && !gStop; rep++) {
          const std::atomic<bool> *stop = best.ok ? &gStop : nullptr;  // the first run always completes
          R->setStop(stop);
          uint64_t sd = runSeed(runNo++);
          R->setSeed(sd);
          gMergeRestarts = std::max(baseRestarts, 4 * w);  // -2: wider steps, bigger budgets, more restarts
          Result r = solveOnce(*R, o, initial, w, stop, keep, sd);
          if (r.aborted || !r.ok) { gStop = true; break; }
          bool better = !best.ok || r.total() < best.total();
          printf("  width %5d: total %8ld  (edges %ld, corners+middle edges %ld, centres %ld)  %.2fs%s\n", w, r.total(),
                 r.wings, r.corners, r.pairs + r.diagmid, r.tTotal, better ? "  best" : "");
          fflush(stdout);
          if (better) best = std::move(r);
        }
        if (gStop) break;
        int nw = (int)std::lround(w * o.grow);
        if (nw <= w) nw = w + 1;
        if (w < maxW && nw > maxW) nw = maxW;
        else if (w >= maxW) break;
        w = nw;
      }
      printf("best: %ld moves\n", best.total());
      printDetail(o, best);
    } else {
    printf("N=%d: widths 1, then x%.2f each step%s; best solution per pair kept across widths (Ctrl-C to stop)\n", N,
           o.grow, widthGiven ? (", up to " + std::to_string(maxW)).c_str() : "");
    EdgesOut E = runEdges(*R, o, initial, o.wingWidth, nullptr, keep, runSeed(-1));
    if (!E.ok) return 1;
    printf("  edges: corners+middle edges %ld, wings %ld (wing width %d)  %.2fs\n", E.corners, E.wings, o.wingWidth,
           E.tCorners + E.tWings);
    std::vector<PairSol> bestPairs;
    long runNo = 0;
    for (int w = 1; w <= maxW;) {
     for (int rep = 0; rep < o.reps && !gStop; rep++) {
      const std::atomic<bool> *stop = best.ok ? &gStop : nullptr;  // the first run always completes
      R->setStop(stop);
      R->setSeed(runSeed(runNo++));  // a different seed per run: different collisions and tie-breaks, different paths
      double ts = now();
      std::vector<PairSol> pj;
      long retries = 0, fallback = 0;
      if (!solvePairs(*R, o, E.cube, w, stop, pj, retries, fallback)) { gStop = true; break; }
      long thisW = 0, improved = 0;
      for (auto &j : pj) thisW += j.cost;
      if (bestPairs.empty()) bestPairs = pj;
      else
        for (size_t i = 0; i < pj.size(); i++)
          if (pj[i].cost < bestPairs[i].cost) { bestPairs[i] = std::move(pj[i]); improved++; }
      long combined = 0;
      for (auto &j : bestPairs) combined += j.cost;
      Result r = finishCentres(*R, o, E, bestPairs, w, stop, keep);
      if (r.aborted || !r.ok) { gStop = true; break; }
      bool better = !best.ok || r.total() < best.total();
      printf("  width %5d: pairs this width %10ld, best-of %10ld (%ld pairs improved)  total %10ld  (edges %ld, "
             "corners+middle edges %ld, centres %ld)  %.2fs%s\n",
             w, thisW, combined, improved, r.total(), r.wings, r.corners, r.pairs + r.diagmid, now() - ts,
             better ? "  best" : "");
      fflush(stdout);
      if (better) best = std::move(r);
     }
      if (gStop) break;
      int nw = (int)std::lround(w * o.grow);
      if (nw <= w) nw = w + 1;
      if (w < maxW && nw > maxW) nw = maxW;  // finish exactly at the cap
      else if (w >= maxW) break;
      w = nw;
    }
    printf("best: %ld moves\n", best.total());
    printDetail(o, best);
    }
  }

  conventionalMoves(best.moves, N);  // what -o prints and --verify replays
  if (o.writeMoves) {
    std::string line;
    for (size_t i = 0; i < best.moves.size(); i++) {
      const mv &m = best.moves[i];
      line += ' ';
      if (m.dep > 1) line += std::to_string(m.dep);
      line += "ULFRBD"[m.face];
      if (m.twist == 2) line += '2';
      if (m.twist == 3) line += '\'';
      if (line.size() > 100 || i + 1 == best.moves.size()) { printf("%s\n", line.c_str()); line.clear(); }
    }
  }
  if (o.verify) {
    if (!fromMoves && !o.readFacelets) { printf("  (--verify needs a move scramble or facelets: -i, -f, -r or -R)\n"); return !best.solved; }
    // independent check: replay scramble (or load the facelets) + solution on the
    // pair-project sticker simulator
    Cube sim(N);
    std::vector<uint8_t> col(sim.pos.size());
    for (size_t i = 0; i < col.size(); i++) col[i] = sim.face[i];
    if (o.readFacelets) {  // the simulator's own facelet geometry (Cube::kociemba)
      std::string s;
      for (char c : facelets) if (c > ' ') s += c;
      for (int kf = 0, k = 0; kf < 6; kf++)
        for (int r = 0; r < N; r++)
          for (int c = 0; c < N; c++, k++) col[sim.kociemba(kf, r, c)] = face_of_char(s[k]);
    }
    auto play = [&](const std::vector<mv> &ms) {
      for (auto m : ms) {
        const auto &P = sim.movePerm[m.face][m.dep - 1];
        for (int k = 0; k < m.twist; k++) {
          std::vector<uint8_t> nc(col.size());
          for (size_t i = 0; i < col.size(); i++) nc[P[i]] = col[i];
          col.swap(nc);
        }
      }
    };
    play(scramble);
    play(best.moves);
    bool ok = true;
    std::map<std::pair<int, int>, int> wrong;  // (row depth, col depth) -> wrong stickers
    for (size_t i = 0; i < col.size(); i++)
      if (col[i] != sim.face[i]) { ok = false; wrong[{sim.drow(i), sim.dcol(i)}]++; }
    printf("  independent sticker-simulator check: %s (%zu moves replayed)\n", ok ? "SOLVED" : "NOT SOLVED",
           best.moves.size());
    for (auto &w : wrong) printf("    wrong stickers at depth (%d,%d): %d\n", w.first.first, w.first.second, w.second);
  }
  return best.solved ? 0 : 1;
}
