#pragma once
#include <algorithm>
#include "util.h"

using uchar = unsigned char ;

struct xorbit {
   uchar a[24] ;
   void rot4(int b, int c, int d, int e) {
      uchar t = a[b] ;
      a[b] = a[e] ;
      a[e] = a[d] ;
      a[d] = a[c] ;
      a[c] = t ;
   }
   void rot2(int b, int c, int d, int e) {
      std::swap(a[b], a[d]) ;
      std::swap(a[c], a[e]) ;
   }
   void rot4i(int b, int c, int d, int e) {
      rot4(e, d, c, b) ;
   }
} ;

inline void oclear(xorbit& x) {
   for (int i=0; i<24; i++)
      x.a[i] = i >> 2 ;
}
inline void uclear(xorbit& x) {
   for (int i=0; i<24; i++)
      x.a[i] = i ;
}
inline void scramble(xorbit& x) {
   oclear(x) ;
   std::shuffle(x.a, x.a+24, rng) ;
}
inline void permimul(xorbit &a, xorbit &b) {
   xorbit tmp ;
   for (int i=0; i<24; i++)
      tmp.a[i] = a.a[b.a[i]] ;
   a = tmp ;
}
