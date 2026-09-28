#include "cornersolver.h"
using namespace std ;

void solvecorners(xcube &cube) {
   int odd = (int)((cube.Mx + cube.My) & 1) ;
   threecube cc = cube.tc ;
   static int inited = 0 ;
   vector<mv> cornersolution ;
   if (odd) {
      ll midslice = (ll)(cube.Mx + 1) ;
      for (int mvs=0; mvs<2; mvs++) {
         int centerwrongcount = cc.countbadcenters() ;
         if (centerwrongcount) {
            int found = 0 ;
            for (int f=0; found==0 && f<3; f++)
               for (int tw=1; found==0 && tw<4; tw++) {
                  cc.domove({2, (short)f, (short)tw}) ;
                  if (cc.countbadcenters() < centerwrongcount) {
                     found = 1 ;
                     cornersolution.push_back({(int)midslice, (short)f, (short)tw}) ;
                  } else {
                     cc.domove({2, (short)f, (short)(4-tw)}) ;
                  }
               }
         }
      }
      if (cc.countbadcenters())
         error("! impossible centercenter state") ;
      if (cornersolution.size()) {
         cout << "Oriented center centers in " << cornersolution.size() << " moves." << endl ;
         for (auto &mv: cornersolution)
            domove(mv, cube) ;
      }
      static struct prunetable<7, 10000019, 12> phase1prune("phase1") ;
      static struct prunetable<8, 50000017, 18, 1> phase2prune("phase2") ;
      if (!inited) {
         inited = 1 ;
         threecube c1('1', 1) ;
         phase1prune.fillit(c1) ;
         threecube c2('C', 1) ;
         phase2prune.fillit(c2) ;
      }
      threecube c1 = cc.tophase1() ;
      int r1 = phase1prune.solve(c1) ;
      if (r1 < 0)
         error("! failed to solve phase 1") ;
      cornersolution.insert(cornersolution.end(),
                            phase1prune.mvs.begin(), phase1prune.mvs.end()) ;
      for (auto &mv: phase1prune.mvs)
         domove(mv, cube) ;
      cc = cube.tc ;
      int r2 = phase2prune.solve(cc) ;
      if (r2 < 0)
         error("! failed to solve phase 2") ;
      for (auto &mv: phase2prune.mvs)
         domove(mv, cube) ;
      cornersolution.insert(cornersolution.end(),
                            phase2prune.mvs.begin(), phase2prune.mvs.end()) ;
      cout << "Solved corners plus middle edges in " << r1 << " + " << r2
                               << " = " << (r1 + r2) << " moves." << endl ;
   } else {
      static struct prunetable<6, 2000003, 11> cornerprune("corner") ;
      if (!inited) {
         inited = 1 ;
         threecube c1('C', 0) ;
         cornerprune.fillit(c1) ;
      }
      int r = cornerprune.solve(cc) ;
      if (r < 0)
         error("! failed to solve corners") ;
      for (auto &mv: cornerprune.mvs)
         domove(mv, cube) ;
      cornersolution = cornerprune.mvs ;
      cout << "Solved corners in " << cornerprune.mvs.size() << " moves." << endl ;
   }
   runningmovecount += cornersolution.size() ;
   if (writemoves) {
      for (auto mv: cornersolution) {
         cout << " " ;
         showmove(mv) ;
      }
      cout << endl ;
   }
}
