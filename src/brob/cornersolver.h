#pragma once
#include <cstring>
#include <iostream>
#include "xcube.h"

double duration() ;
extern ll runningmovecount ;
extern int writemoves ;

template<int> int goodmove(int, int) { return 1 ; }
template<> inline int goodmove<1>(int f, int tw) {
   return f==0 || f==5 || tw == 2 ;
}

template<int dep, int size, int maxdepth, int limitmoves=0> struct prunetable {
   prunetable(const char *name_) : solved('C', 0), name(name_) {
      memset(prune, dep, size) ;
   }
   void filltable(threecube &cc, int togo, int prevf, int d) {
      if (togo == 0) {
         int h = cc.hash() % size ;
         if (prune[h] == dep) {
            prune[h] = d ;
         }
      } else {
         for (int f=0; f<6; f++)
            if (f != prevf) {
               for (int tw=1; tw<4; tw++) {
                  cc.domove({1, (short)f, 1}) ;
                  if (goodmove<limitmoves>(f, tw) &&
                      (d != togo || cc != solved))
                     filltable(cc, togo-1, f, d) ;
               }
               cc.domove({1, (short)f, 1}) ;
            }
      }
   }
   void fillit(threecube &cc) {
      duration() ;
      solved = cc ;
      for (int d=0; d<dep; d++)
         filltable(cc, d, -1, d) ;
      std::cout << "Created pruning table " << name << " in " << duration() << std::endl ;
   }
   int lookup(const threecube &cc) {
      return prune[cc.hash() % size] ;
   }
   int solve(threecube &cc, int togo, int prevf) {
      if (togo == 0)
         return (cc == solved) ;
      if (lookup(cc) > togo)
         return 0 ;
      for (int f=0; f<6; f++)
         if (f != prevf) {
            for (int tw=1; tw<4; tw++) {
               mvs.push_back({1, (short)f, (short)tw}) ;
               cc.domove({1, (short)f, 1}) ;
               if (goodmove<limitmoves>(f, tw) && solve(cc, togo-1, f))
                  return 1 ;
               mvs.pop_back() ;
            }
            cc.domove({1, (short)f, 1}) ;
         }
      return 0 ;
   }
   int solve(threecube &cc) {
      duration() ;
      mvs.clear() ;
      for (int curdep=lookup(cc); curdep<=maxdepth; curdep++) {
         if (solve(cc, curdep, -1)) {
            std::cout << "Solved " << name << " in " << duration() << std::endl ;
            return curdep ;
         }
      }
      return -1 ;
   }
   unsigned char prune[size] ;
   threecube solved ;
   std::vector<mv> mvs ;
   const char *name ;
} ;

void solvecorners(xcube &cube) ;
