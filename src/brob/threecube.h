#pragma once
#include <cassert>
#include <string_view>
#include "xorbit.h"
#include "mv.h"
#include "util.h"

/*
 *   This does 333 work, and 222 work with vestigial edges.
 */
struct threecube {
   static constexpr int cor[8][3] = {
      {0, 4, 17}, {1, 16,  13},  {2, 12, 9}, {3, 8, 5},
      {20, 6, 11}, {23, 18, 7}, {21, 10, 15}, {22, 14, 19}
   } ;
   static constexpr int edg[12][2] = {
      {0, 16}, {1, 12}, {2, 8}, {3, 4}, {11, 5}, {9, 15}, {19, 13},
      {17, 7}, {20, 10}, {21, 14}, {22, 18}, {23, 6}
   } ;
   threecube(int src_, int odd) {
      for (int i=0; i<6; i++)
         centercenters[i] = i ;
      if (src_ == 'I' || src_ == 'C') { // inputmoves, or clear
         oclear(corner) ;
         oclear(edge) ;
         if (!odd)
            for (int i=0; i<24; i++)
               edge.a[i] = 0 ;
      } else if (src_ == 'U') { // orbits have numbers 0..23
         uclear(corner) ;
         uclear(edge) ;
      } else if (src_ == 'R') { // random
         oclear(corner) ;
         if (odd)
            oclear(edge) ;
         else
            for (int i=0; i<24; i++)
               edge.a[i] = 0 ;
         int parity = 0 ;
         for (int i=0; i<8; i++) {
            std::uniform_int_distribution<int> distrib(i, 7) ;
            int j = distrib(rng) ;
            if (j != i) {
               parity++ ;
               std::swap(corner.a[cor[i][0]], corner.a[cor[j][0]]) ;
               std::swap(corner.a[cor[i][1]], corner.a[cor[j][1]]) ;
               std::swap(corner.a[cor[i][2]], corner.a[cor[j][2]]) ;
            }
         }
         if (odd) {
            for (int i=0; i<11; i++) {
               std::uniform_int_distribution<int> distrib(i, 11) ;
               int j = distrib(rng) ;
               if (j != i) {
                  parity++ ;
                  std::swap(edge.a[edg[i][0]], edge.a[edg[j][0]]) ;
                  std::swap(edge.a[edg[i][1]], edge.a[edg[j][1]]) ;
               }
            }
            if (parity & 1) {
               std::swap(edge.a[edg[10][0]], edge.a[edg[11][0]]) ;
               std::swap(edge.a[edg[10][1]], edge.a[edg[11][1]]) ;
            }
         }
         int twist = 0 ;
         std::uniform_int_distribution<int> tw(0, 2) ;
         for (int i=0; i<8; i++) {
            int j ;
            if (i<7) {
               j = tw(rng) ;
               twist += j ;
            } else
               j = (3 * 8 - twist) % 3 ;
            switch (j) {
case 0:        break ;
case 1:        std::swap(corner.a[cor[i][0]], corner.a[cor[i][1]]) ;
               std::swap(corner.a[cor[i][1]], corner.a[cor[i][2]]) ;
               break ;
case 2:        std::swap(corner.a[cor[i][1]], corner.a[cor[i][2]]) ;
               std::swap(corner.a[cor[i][0]], corner.a[cor[i][1]]) ;
               break ;
            }
         }
         if (odd) {
            twist = 0 ;
            std::uniform_int_distribution<int> tw2(0, 1) ;
            for (int i=0; i<12; i++) {
               int j ;
               if (i<11) {
                  j = tw2(rng) ;
                  twist += j ;
               } else
                  j = (twist & 1) ;
               if (j)
                  std::swap(edge.a[edg[i][0]], edge.a[edg[i][1]]) ;
            }
         }
      } else if (src_ == '1') { // phase 1
         for (int i=0; i<24; i++) {
            corner.a[i] = 0 ;
            edge.a[i] = 0 ;
         }
         for (int i=0; i<4; i++) {
            corner.a[i] = 1 ;
            edge.a[i] = 1 ;
            corner.a[20+i] = 1 ;
            edge.a[20+i] = 1 ;
         }
         edge.a[11] = edge.a[9] = edge.a[19] = edge.a[17] = 2 ;
      }
   }
   void rot4c(int a, int b, int c, int d) {
      auto t = centercenters[a] ;
      centercenters[a] = centercenters[b] ;
      centercenters[b] = centercenters[c] ;
      centercenters[c] = centercenters[d] ;
      centercenters[d] = t ;
   }
   int countbadcenters() {
      int centerwrongcount = 0 ;
      for (int i=0; i<6; i++)
         if (centercenters[i] != i)
            centerwrongcount++ ;
      return centerwrongcount ;
   }
   void domove(mv mv) {
      assert(mv.dep == 1 || mv.dep == 2) ;
      assert(mv.twist > 0) ;
      assert(mv.twist < 4) ;
      for (int i=0; i<mv.twist; i++) {
         switch (mv.face + 6 * mv.dep - 6) {
case 0: // U
            corner.rot4(0, 1, 2, 3) ;
            corner.rot4(16, 12, 8, 4) ;
            corner.rot4(17, 13, 9, 5) ;
            edge.rot4(0, 1, 2, 3) ;
            edge.rot4(16, 12, 8, 4) ;
            break ;
case 1: // L
            corner.rot4(4, 5, 6, 7) ;
            corner.rot4(0, 8, 20, 18) ;
            corner.rot4(3, 11, 23, 17) ;
            edge.rot4(4, 5, 6, 7) ;
            edge.rot4(3, 11, 23, 17) ;
            break ;
case 2: // F
            corner.rot4(8, 9, 10, 11) ;
            corner.rot4(3, 12, 21,  6) ;
            corner.rot4(2, 15, 20, 5) ;
            edge.rot4(8, 9, 10, 11) ;
            edge.rot4(2, 15, 20, 5) ;
            break ;
case 3: // R
            corner.rot4(12, 13, 14, 15) ;
            corner.rot4(22, 10, 2, 16) ;
            corner.rot4(21, 9, 1, 19) ;
            edge.rot4(12, 13, 14, 15) ;
            edge.rot4(1, 19, 21, 9) ;
            break ;
case 4: // B
            corner.rot4(16, 17, 18, 19) ;
            corner.rot4(1, 4, 23, 14) ;
            corner.rot4(0, 7, 22, 13) ;
            edge.rot4(16, 17, 18, 19) ;
            edge.rot4(0, 7, 22, 13) ;
            break ;
case 5: // D
            corner.rot4(20, 21, 22, 23) ;
            corner.rot4(7, 11, 15, 19) ;
            corner.rot4(6, 10, 14, 18) ;
            edge.rot4(20, 21, 22, 23) ;
            edge.rot4(10, 14, 18, 6) ;
            break ;
case 6: // 2U
            edge.rot4(19, 15, 11, 7) ;
            edge.rot4(17, 13, 9, 5) ;
            rot4c(1, 2, 3, 4) ;
            break ;
case 7: // 2L
            edge.rot4(0, 8, 20, 18) ;
            edge.rot4(2, 10, 22, 16) ;
            rot4c(4, 5, 2, 0) ;
            break ;
case 8: // 2F
            edge.rot4(1, 14, 23, 4) ;
            edge.rot4(3, 12, 21, 6) ;
            rot4c(1, 5, 3, 0) ;
            break ;
case 9: // 2R
            edge.rot4(18, 20, 8, 0) ;
            edge.rot4(16, 22, 10, 2) ;
            rot4c(0, 2, 5, 4) ;
            break ;
case 10: // 2B
            edge.rot4(4, 23, 14, 1) ;
            edge.rot4(6, 21, 12, 3) ;
            rot4c(0, 3, 5, 1) ;
            break ;
case 11: // 2D
            edge.rot4(7, 11, 15, 19) ;
            edge.rot4(5, 9, 13, 17) ;
            rot4c(4, 3, 2, 1) ;
            break ;
         }
      }
   }
   threecube tophase1() const {
      threecube r('C', 1) ;
      for (int i=0; i<24; i++) {
         r.edge.a[i] = 0 ;
         r.corner.a[i] = 0 ;
      }
      for (int i=0; i<12; i++) {
         int a = edg[i][0] ;
         int b = edg[i][1] ;
         if (edge.a[a] == 0 || edge.a[a] == 5) {
            r.edge.a[a] = 1 ;
         } else if (edge.a[b] == 0 || edge.a[b] == 5) {
            r.edge.a[b] = 1 ;
         } else if (edge.a[a] == 2 || edge.a[a] == 4) {
            r.edge.a[a] = 2 ;
         } else if (edge.a[b] == 2 || edge.a[b] == 4) {
            r.edge.a[b] = 2 ;
         } else {
            error("! Unable to convert cube edges to phase1") ;
         }
      }
      for (int i=0; i<24; i++)
         if (corner.a[i] == 0 || corner.a[i] == 5)
            r.corner.a[i] = 1 ;
      return r ;
   }
   // we don't currently hash centercenters because we solve those
   // first and then don't use slice moves for corner/middle edge
   // solves.
   size_t hash() const {
      std::string_view sv((const char *)this, sizeof(*this)) ;
      return std::hash<std::string_view>{}(sv) ;
   }
   bool operator==(const threecube &b) const {
      for (int i=0; i<24; i++)
         if (corner.a[i] != b.corner.a[i] || edge.a[i] != b.edge.a[i])
            return 0 ;
      return 1 ;
   }
   xorbit corner, edge ;
   uchar centercenters[6] ;
} ;
