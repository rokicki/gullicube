#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <vector>
#include <string>
#include <iostream>
#include "xorbit.h"
#include "mv.h"
#include "threecube.h"
#include "util.h"

using ll = long long ;
using cubesizetype = ll ;

extern int rots[24][6] ;
extern const char *faceorder ;
extern const int oppface[] ;
extern std::vector<mv> moveset ;
extern std::vector<mv> prevmvs ;
extern std::vector<mv> *g_prep_recorder ;
extern int alwaysinmem ;
extern ll cancellations ;

void coordtospeffz(int Mx, int My, int x, int y, int &ox, int &oy, int &oz) ;
void border() ;
int encoderot(int *p) ;
void expandrotations() ;

struct xcube ;
void domove(mv mv, xcube &cube, bool suppress_cancellations=false) ;

struct xcube {
   xcube(cubesizetype Mx_, cubesizetype My_, int src_) :
  Mx(Mx_), My(My_), rows(My), tc(src_, (Mx_ + My_) & 1), src(src_) {
      for (int i=0; i<6; i++)
         facerots[i] = 0 ;
      for (int i=0; i<6; i++)
         ccs[i] = i ;
      for (int a=0; a<My_; a++){
         xorbit xo ;
         uclear(xo) ;
         if (src_ == 'R')
            std::shuffle(xo.a, xo.a+24, rng) ;
         wingorbits.push_back(xo) ;
      }
      if (src_ == 'U' || alwaysinmem) {
         for (int i=0; i<My; i++) {
            rows[i].resize(Mx) ;
            if (src == 'I' || src == 'C') {
               for (auto &xo: rows[i])
                  oclear(xo) ;
            } else if (src == 'U') {
               for (auto &xo: rows[i])
                  uclear(xo) ;
            } else if (src == 'R') {
               for (auto &xo: rows[i])
                  scramble(xo) ;
            }
         }
      }
      if (src == 'I')
         for (auto mv: moveset)
            domove(mv, *this) ;
   }
   void normalizefacerots() {
      for (int f=0; f<6; f++) {
         if (facerots[f] == 0)
            continue ;
         for (auto &row: rows)
            for (auto &xo: row) {
               switch (facerots[f]) {
                  case 1: xo.rot4i(4*f, 4*f+1, 4*f+2, 4*f+3) ; break ;
                  case 2: xo.rot2(4*f, 4*f+1, 4*f+2, 4*f+3) ; break ;
                  case 3: xo.rot4(4*f, 4*f+1, 4*f+2, 4*f+3) ; break ;
               }
            }
         facerots[f] = 0 ;
      }
   }
   char get(int f, int x, int y) const {
      int ox, oy, oz ;
      coordtospeffz(Mx, My, x, y, ox, oy, oz) ;
      if (oz == 4)
         return '0'+ccs[f] ;
      if (rows[oy-1].size() == 0)
         return '?' ;
      return (char)('0' + rows[oy-1][ox-1].a[4*f+((oz+facerots[f])&3)]) ;
   }
   void wings_cycle(int p1, int p2, int p3, int p4, int dep) {
      auto t = wingorbits[dep].a[p1] ;
      wingorbits[dep].a[p1] = wingorbits[dep].a[p4] ;
      wingorbits[dep].a[p4] = wingorbits[dep].a[p3] ;
      wingorbits[dep].a[p3] = wingorbits[dep].a[p2] ;
      wingorbits[dep].a[p2] = t ;
   }
   void wings_cycle2(int p1, int p2, int p3, int p4, int dep) {
      std::swap(wingorbits[dep].a[p1], wingorbits[dep].a[p3]) ;
      std::swap(wingorbits[dep].a[p2], wingorbits[dep].a[p4]) ;
   }
   void rot4hor(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      if (rows[dep].empty())
         return ;
      for (int i=0; i<Mx; i++) {
         uchar &c1 = rows[dep][i].a[4*f1+((o1+facerots[f1])&3)] ;
         uchar &c2 = rows[dep][i].a[4*f2+((o2+facerots[f2])&3)] ;
         uchar &c3 = rows[dep][i].a[4*f3+((o3+facerots[f3])&3)] ;
         uchar &c4 = rows[dep][i].a[4*f4+((o4+facerots[f4])&3)] ;
         auto t = c1 ; c1 = c2 ; c2 = c3 ; c3 = c4 ; c4 = t ;
      }
   }
   void rot4ver(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      for (int i=0; i<My; i++) {
         if (rows[i].empty())
            continue ;
         uchar &c1 = rows[i][dep].a[4*f1+((o1+facerots[f1])&3)] ;
         uchar &c2 = rows[i][dep].a[4*f2+((o2+facerots[f2])&3)] ;
         uchar &c3 = rows[i][dep].a[4*f3+((o3+facerots[f3])&3)] ;
         uchar &c4 = rows[i][dep].a[4*f4+((o4+facerots[f4])&3)] ;
         auto t = c1 ; c1 = c2 ; c2 = c3 ; c3 = c4 ; c4 = t ;
      }
   }
   void rot4ver2(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      rot4ver(f1, o1, f2, o2, f3, o3, f4, o4, dep) ;
      rot4ver(f1, o1^2, f2, o2^2, f3, o3^2, f4, o4^2, dep) ;
      uchar &c1 = ccs[f1] ;
      uchar &c2 = ccs[f2] ;
      uchar &c3 = ccs[f3] ;
      uchar &c4 = ccs[f4] ;
      auto t = c1 ; c1 = c2 ; c2 = c3 ; c3 = c4 ; c4 = t ;
   }
   void rot2hor(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      if (rows[dep].empty())
         return ;
      for (int i=0; i<Mx; i++) {
         uchar &c1 = rows[dep][i].a[4*f1+((o1+facerots[f1])&3)] ;
         uchar &c2 = rows[dep][i].a[4*f2+((o2+facerots[f2])&3)] ;
         uchar &c3 = rows[dep][i].a[4*f3+((o3+facerots[f3])&3)] ;
         uchar &c4 = rows[dep][i].a[4*f4+((o4+facerots[f4])&3)] ;
         std::swap(c1, c3) ;
         std::swap(c2, c4) ;
      }
   }
   void rot2ver(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      for (int i=0; i<My; i++) {
         if (rows[i].empty())
            continue ;
         uchar &c1 = rows[i][dep].a[4*f1+((o1+facerots[f1])&3)] ;
         uchar &c2 = rows[i][dep].a[4*f2+((o2+facerots[f2])&3)] ;
         uchar &c3 = rows[i][dep].a[4*f3+((o3+facerots[f3])&3)] ;
         uchar &c4 = rows[i][dep].a[4*f4+((o4+facerots[f4])&3)] ;
         std::swap(c1, c3) ;
         std::swap(c2, c4) ;
      }
   }
   void rot2ver2(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      rot2ver(f1, o1, f2, o2, f3, o3, f4, o4, dep) ;
      rot2ver(f1, o1^2, f2, o2^2, f3, o3^2, f4, o4^2, dep) ;
      uchar &c1 = ccs[f1] ;
      uchar &c2 = ccs[f2] ;
      uchar &c3 = ccs[f3] ;
      uchar &c4 = ccs[f4] ;
      std::swap(c1, c3) ;
      std::swap(c2, c4) ;
   }
   void getrow(int i) {
      if (alwaysinmem)
         return ;
      if (src != 'R')
         error("! only support random right now for not in memory") ;
      rows[i].resize(Mx) ;
      for (auto &orb: rows[i])
         scramble(orb) ;
   }
   void checkrow(int i) const {
      if ((ll)rows[i].size() != Mx)
         error("! checking a nonexistent row?") ;
      for (auto &orb: rows[i])
         for (int j=0; j<24; j++)
            if (orb.a[j] != (j>>2)) {
               show() ;
               error("! cube is not solved") ;
            }
   }
   void checkrelease(int i) {
      if (alwaysinmem)
         return ;
      checkrow(i) ;
      std::vector<xorbit> empty ;
      std::swap(empty, rows[i]) ;
   }
   int wing_color(int oz_wings, int f, int j, int N) {
      static constexpr int colors_of_wings[][2] = {
         {0,4}, {0,3}, {0,2}, {0,1}, {1,0}, {1,2}, {1,5}, {1,4},
         {2,0}, {2,3}, {2,5}, {2,1}, {3,0}, {3,4}, {3,5}, {3,2},
         {4,0}, {4,1}, {4,5}, {4,3}, {5,2}, {5,3}, {5,4}, {5,1}
      } ;
      static constexpr int PARTNER[24] = {16,12,8,4,  3,11,23,17, 2,15,20,5,  1,19,21,9, 0,7,22,13, 10,14,18,6} ;
      int k = (j <= N/2) ? (j-1) : (N-j) ;
      int slot = wingorbits[k].a[4*f + oz_wings] ;
      if ((j <= N/2)) return colors_of_wings[slot][0] ;
      else return colors_of_wings[ wingorbits[k].a[ PARTNER[4*f + oz_wings] ] ][1] ;
   }
   void show_wings() const {
      for (int i=0; i<(int)wingorbits.size(); i++) {
         std::cout << "edge orbit index " << i+1 << ": [" ;
         for (int j=0; j<24; j++) {
            if (j) std::cout << "," ;
            std::cout << (int)wingorbits[i].a[j] ;
         }
         std::cout << "]" << std::endl ;
      }
   }
   std::vector<std::vector<int>> get_wings(int i) {
      std::vector<std::vector<int>> r ;
      std::vector<int> orbit ;
      for (int j=0; j<24; j++)
         orbit.push_back(wingorbits[i].a[j]) ;
      r.push_back(orbit) ;
      return r ;
   }
   cubesizetype Mx, My ;
   std::vector<std::vector<xorbit>> rows ;
   std::vector<xorbit> wingorbits ;
   threecube tc ;
   void show() const ;
   int facerots[6] ;
   uchar ccs[6] ;
   char src ;
} ;

mv parse_move(const char *&s) ;
std::vector<mv> parse_moves(const std::string &alg) ;
void readmoveset() ;
int trymerge(const mv &a, const mv &b, mv &res) ;
void trycancellation(mv &m) ;
void showmove(const mv &mv) ;
void random_move_scramble(cubesizetype N, ll length) ;
