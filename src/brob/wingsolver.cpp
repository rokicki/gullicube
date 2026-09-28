#include "wingsolver.h"
#include <algorithm>
using namespace std ;

#if WRITEFRAMES
extern ll writemoveinc, writemovetarget ;
extern void writeframe(xcube &) ;
#endif

int wing_cycle_index[24][24][24] ;
int wing_cycle_triple[4048][3] ;

void build_wing_cycle_index() {
   for (int a=0; a<24; a++)
      for (int b=0; b<24; b++)
         for (int c=0; c<24; c++)
            wing_cycle_index[a][b][c] = -1 ;
   int idx = 0 ;
   for (int a=0; a<24; a++)
      for (int b=0; b<24; b++)
         for (int c=0; c<24; c++) {
            if (a==b || b==c || a==c)
               continue ;
            if (a < b && a < c) {
               wing_cycle_index[a][b][c] = idx ;
               wing_cycle_index[b][c][a] = idx ;
               wing_cycle_index[c][a][b] = idx ;
               wing_cycle_triple[idx][0] = a ;
               wing_cycle_triple[idx][1] = b ;
               wing_cycle_triple[idx][2] = c ;
               idx++ ;
            }
         }
   if (idx != 4048)
      error("error, not 4048") ;
}

static vector<int> solve_wings(vector<int> orb) {
   vector<int> algs ;
   for (int i=0; i<24; i++) {
      while (orb[i] != i) {
         int t1 = find(orb.begin(), orb.end(), i) - orb.begin() ;
         int t2 = i ;
         int t3 = orb[i] ;
         if (t3 == t1) {
            bool seen[24] = {} ;
            int cyclebreak = -1 ;
            for (int s=i+1; s<24 && cyclebreak < 0; s++) {
               if (orb[s] == s || s == t1 || seen[s])
                  continue ;
               int len = 0 ;
               int x = s ;
               do { seen[x] = true ; len++ ; x = orb[x] ; } while (x != s) ;
               if ((len % 2) == 0)
                  cyclebreak = s ;
            }
            t3 = cyclebreak ;
         }
         algs.push_back(wing_cycle_index[t1][t2][t3]) ;
         int tmp = orb[t3] ;
         orb[t3] = orb[t2] ;
         orb[t2] = orb[t1] ;
         orb[t1] = tmp ;
      }
   }
   return algs ;
}

static vector<mv> wing_alg_moves[4048] ;

static int wing_cycle_of(const vector<mv> &moves) {
   xcube cube(4, 5, 'U') ;
   for (auto mv : moves)
      domove(mv, cube) ;
   xorbit &o = cube.wingorbits[0] ;
   int dest[24] ;
   for (int t=0; t<24; t++)
      dest[o.a[t]] = t ;
   int t1 = -1 ;
   int moved = 0 ;
   for (int s=0; s<24; s++)
      if (dest[s] != s) {
         moved++ ;
         if (t1 < 0)
            t1 = s ;
      }
   if (moved != 3)
      return -1 ;
   return wing_cycle_index[t1][dest[t1]][dest[dest[t1]]] ;
}

static vector<mv> rotate_moves(const vector<mv> &in, int m, bool inv) {
   vector<mv> r = in ;
   if (inv) {
      reverse(r.begin(), r.end()) ;
      for (auto &mv : r)
         mv.twist = (4 - mv.twist) & 3 ;
   }
   for (auto &mv : r)
      mv.face = rots[m][mv.face] ;
   return r ;
}

static const char *wing_base_algs[] = {
   "F D 2B' D' F' D 2B D'",
   "F' L' 2B L F L' 2B' L",
   "F D 2B2 D' F' D 2B2 D'",
   "F' L' 2B2 L F L' 2B2 L",
   "F D 2B D' F' D 2B' D'",
   "F' L' 2B' L F L' 2B L",
   "F D 2F' D' F' D 2F D'",
   "F' L' 2F L F L' 2F' L",
   "F D 2F2 D' F' D 2F2 D'",
   "F' L' 2F2 L F L' 2F2 L",
   "F D2 F' 2D F D2 F' 2D'",
   "F' L2 F 2L' F' L2 F 2L",
   "F D2 F' 2D2 F D2 F' 2D2",
   "F' L2 F 2L2 F' L2 F 2L2",
   "F D2 F' 2D' F D2 F' 2D",
   "F' L2 F 2L F' L2 F 2L'",
   "F D2 F' 2U F D2 F' 2U'",
   "F' L2 F 2R' F' L2 F 2R",
   "F D' 2B2 D F' D' 2B2 D",
   "F' L 2B2 L' F L 2B2 L'",
   "F D2 F' 2U2 F D2 F' 2U2",
   "F' L2 F 2R2 F' L2 F 2R2",
   "F D' 2F2 D F' D' 2F2 D",
   "F' L 2F2 L' F L 2F2 L'",
   "F D2 F' 2U' F D2 F' 2U",
   "F' L2 F 2R F' L2 F 2R'",
   "F D' F' 2D2 F D F' 2D2",
   "F' L F 2L2 F' L' F 2L2",
   "F D' F' 2D' F D F' 2D",
   "F' L F 2L F' L' F 2L'",
   "F D' F' 2U F D F' 2U'",
   "F' L F 2R' F' L' F 2R",
   "F D' F' 2U2 F D F' 2U2",
   "F' L F 2R2 F' L' F 2R2",
   "F 2D F' D2 F 2D' F' D2",
   "F' 2L' F L2 F' 2L F L2",
   "F 2D F' U2 F 2D' F' U2",
   "F' 2L' F R2 F' 2L F R2",
   "F 2D2 F' D2 F 2D2 F' D2",
   "F' 2L2 F L2 F' 2L2 F L2",
   "F 2D2 F' U2 F 2D2 F' U2",
   "F' 2L2 F R2 F' 2L2 F R2",
   "F 2D' F' U2 F 2D F' U2",
   "F' 2L F R2 F' 2L' F R2",
   "F D 2B D' F D 2B' D' F2",
   "F' L' 2B' L F' L' 2B L F2",
   "F D 2B' D' F D 2B D' F2",
   "F D 2F D' F D 2F' D' F2",
   "F D 2B D' F2 D 2B' D' F",
   "F' L' 2B' L F2 L' 2B L F'",
   "F D F 2D F' D' F 2D' F2",
   "F' L' F' 2L' F L F' 2L F2",
   "F D2 F 2D F' D2 F 2D' F2",
   "F' L2 F' 2L' F L2 F' 2L F2",
   "F D' F 2D F' D F 2D' F2",
   "F' L F' 2L' F L' F' 2L F2",
   "F D' F 2D2 F' D F 2D2 F2",
   "F' L F' 2L2 F L' F' 2L2 F2",
   "F 2D F D2 F' 2D' F D2 F2",
   "F 2D F U F' 2D' F U' F2",
   "F' 2L' F' R' F 2L F' R F2",
   "F 2D2 F D F' 2D2 F D' F2",
   "F 2D2 F U F' 2D2 F U' F2",
   "F B D F' 2D' F D' F' 2D B'",
   "F B D F' 2U F D' F' 2U' B'",
   "F B D' F' 2D' F D F' 2D B'",
   "F' B' L F 2L F' L' F 2L' B",
   "F B 2D2 F' D F 2D2 F' D' B'",
   "F' B' 2L2 F L' F' 2L2 F L B",
   "F B' D F' 2D' F D' F' 2D B",
   "F' B L' F 2L F' L F 2L' B'",
   "F D R' 2D2 R D' R' 2D2 R F'",
   "F D R' 2U2 R D' R' 2U2 R F'",
   "F D2 L' 2D' L D2 L' 2D L F'",
   "F' L2 D 2L D' L2 D 2L' D' F",
   "F D B2 L' 2F L B2 L' 2F' L D' F'",
   "F D' L F' D 2L D' F L' D F' 2L'",
   "F' L D' F L' 2D' L F' D L' F 2D",
   "F D' L F' D 2L2 D' F L' D F' 2L2",
   "F D' L 2F2 L' D F' L D' 2F2 D L'",
   "F' L D' 2F2 D L' F D' L 2F2 L' D",
   "F D' L 2F' L' D F' L D' 2F D L'",
   "F' L D' 2F D L' F D' L 2F' L' D",
   "F D L' F 2D' F' L D' F L' 2D L F2",
   "F' L' D F' 2L F D' L F' D 2L' D' F2",
   "F D L' F 2U F' L D' F L' 2U' L F2",
   "F' L' D F' 2R' F D' L F' D 2R D' F2",
} ;

void loadwingalgs() {
   for (int idx=0; idx<4048; idx++)
      wing_alg_moves[idx].clear() ;
   for (int i=0; i<87; i++) {
      vector<mv> m0 = parse_moves(wing_base_algs[i]) ;
      for (int m=0; m<24; m++)
         for (int inv=0; inv<2; inv++) {
            vector<mv> v = rotate_moves(m0, m, inv) ;
            int idx = wing_cycle_of(v) ;
            if (idx >= 0 && wing_alg_moves[idx].empty())
               wing_alg_moves[idx] = v ;
         }
   }
   for (int idx=0; idx<4048; idx++)
      if (wing_alg_moves[idx].empty())
         error("! wing base set does not cover all 4048 cycles") ;
}

ll doblock_wings(xcube &cube, int idx, const vector<int> &depths) {
#if WRITEFRAMES
   if (writemoveinc && runningmovecount >= writemovetarget)
      writeframe(cube) ;
#endif
   blocks++ ;
   auto &alg = wing_alg_moves[idx] ;
   int nslice = 0 ;
   for (auto &mv : alg)
      if (mv.dep == 2)
         nslice++ ;
   if (writemoves) {
      for (auto mv : alg) {
         if (mv.dep == 2) {
            for (int k : depths) {
               mv.dep = k + 2 ;
               curblock.push_back(mv) ;
            }
         } else {
            curblock.push_back(mv) ;
         }
      }
      sendblock() ;
   }
   for (auto mv : alg) {
      if (mv.dep == 2) {
         for (int k : depths) {
            mv.dep = k + 2 ;
            domove(mv, cube) ;
         }
      } else {
         domove(mv, cube) ;
      }
   }
   ll mvs = (ll)alg.size() + (ll)nslice * ((ll)depths.size() - 1) ;
   runningmovecount += mvs ;
   return mvs ;
}

ll do_wing_cycles(xcube &cube) {
   ll wr = 0 ;
   vector<vector<int>> seq((size_t)cube.My) ;
   for (ll k=0; k<cube.My; k++)
      seq[(size_t)k] = solve_wings(cube.get_wings((ll)k)[0]) ;
   vector<vector<int>> bucket(4048) ;
   for (ll k=0; k<cube.My; k++)
      for (int idx : seq[k])
         bucket[idx].push_back(k) ;
   for (int v=0; v<4048; v++)
      if (!bucket[v].empty()) {
         ll wm = doblock_wings(cube, v, bucket[v]) ;
         wr += wm ;
         wingmoves += wm ;
      }
   return wr ;
}
