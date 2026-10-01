#include "facelets.h"
#include <cstring>
using namespace std ;

/*
 *   Sticker geometry is that of the usual net, face f of the cube drawn as an
 *   S x S grid (x = column, y = row, 0 .. S-1): corners at the grid corners,
 *   wings and middle edges along the border, centre orbits inside, as in
 *   brobdicube's frame writer.  Centre centres come from xcube::ccs, which
 *   matches physical cubes.
 */
// xcube face (U L F R B D) -> its position in a Kociemba string (U R F D L B)
static const int kocipos[6] = { 0, 4, 2, 1, 5, 3 } ;
// as xcube::wing_color: the colours of wing piece p (sticker 0 on the face
// of the slot it sits in, sticker 1 on the partner slot's face), and each
// wing slot's partner slot
static const int wingcols[24][2] = {
   {0,4}, {0,3}, {0,2}, {0,1}, {1,0}, {1,2}, {1,5}, {1,4},
   {2,0}, {2,3}, {2,5}, {2,1}, {3,0}, {3,4}, {3,5}, {3,2},
   {4,0}, {4,1}, {4,5}, {4,3}, {5,2}, {5,3}, {5,4}, {5,1}
} ;
static const int wingpartner[24] = {
   16,12,8,4,  3,11,23,17, 2,15,20,5,  1,19,21,9, 0,7,22,13, 10,14,18,6
} ;
static int parity(const int *p, int n) {
   int r = 0 ;
   for (int i=0; i<n; i++)
      for (int j=i+1; j<n; j++)
         r ^= p[i] > p[j] ;
   return r ;
}
// border stickers of a face: slot oz (0 top, 1 right, 2 bottom, 3 left) at
// position i along the border, clockwise; corners at the four grid corners
static void borderxy(int N, int oz, int i, int &x, int &y) {
   switch (oz) {
case 0: x = i ; y = 0 ; break ;
case 1: x = N+1 ; y = i ; break ;
case 2: x = N+1-i ; y = N+1 ; break ;
default: x = 0 ; y = N+1-i ; break ;
   }
}
static void cornerxy(int N, int oz, int &x, int &y) {
   x = (oz == 1 || oz == 2) ? N+1 : 0 ;
   y = (oz >= 2) ? N+1 : 0 ;
}
// the four stickers of centre orbit (i, j) on a face, quadrant q
static void centrexy(int N, int i, int j, int q, int &x, int &y) {
   switch (q) {
case 0: x = 1+j ; y = 1+i ; break ;
case 1: x = N-i ; y = 1+j ; break ;
case 2: x = N-j ; y = N-i ; break ;
default: x = 1+i ; y = N-j ; break ;
   }
}
string xcube_to_facelets(const xcube &cube) {
   const int N = (int)(cube.Mx + cube.My), S = N + 2 ;
   string r(6LL*S*S, '?') ;
   auto put = [&](int f, int x, int y, int c) {
      r[(size_t)kocipos[f]*S*S + (size_t)y*S + x] = faceorder[c] ;
   } ;
   int x, y ;
   for (int f=0; f<6; f++) {
      for (int i=0; i<cube.My; i++) {
         if (cube.rows[i].empty())
            continue ;
         for (int j=0; j<cube.Mx; j++)
            for (int q=0; q<4; q++) {
               centrexy(N, i, j, q, x, y) ;
               put(f, x, y, cube.rows[i][j].a[4*f+((q+cube.facerots[f])&3)]) ;
            }
      }
      for (int i=1; i<=N; i++) {
         if (i+i == N+1)
            continue ;
         for (int oz=0; oz<4; oz++) {
            int k = (i <= N/2) ? i-1 : N-i ;
            int slot = 4*f + oz ;
            int c = (i <= N/2) ? wingcols[cube.wingorbits[k].a[slot]][0]
                               : wingcols[cube.wingorbits[k].a[wingpartner[slot]]][1] ;
            borderxy(N, oz, i, x, y) ;
            put(f, x, y, c) ;
         }
      }
      for (int oz=0; oz<4; oz++) {
         cornerxy(N, oz, x, y) ;
         put(f, x, y, cube.tc.corner.a[4*f+oz]) ;
      }
      if (N & 1) {
         int half = (N+1)/2 ;
         for (int oz=0; oz<4; oz++) {
            borderxy(N, oz, half, x, y) ;
            put(f, x, y, cube.tc.edge.a[4*f+oz]) ;
         }
         put(f, half, half, cube.ccs[f]) ;
      }
   }
   return r ;
}
int xcube_from_facelets(xcube &cube, const string &in, string &err) {
   const int N = (int)(cube.Mx + cube.My), S = N + 2 ;
   string s ;
   for (char c: in)
      if (c > ' ')
         s.push_back(c) ;
   if ((ll)s.size() != 6LL*S*S) {
      err = "a " + to_string(S) + "x" + to_string(S) + "x" + to_string(S) + " cube needs " +
            to_string(6LL*S*S) + " facelets, not " + to_string(s.size()) ;
      return 0 ;
   }
   vector<uchar> col(s.size()) ;
   for (size_t i=0; i<s.size(); i++) {
      const char *p = strchr(faceorder, s[i]) ;
      if (p == 0) {
         err = string("facelet ") + to_string(i+1) + " is '" + s[i] + "', not one of URFDLB" ;
         return 0 ;
      }
      col[i] = (uchar)(p - faceorder) ;
   }
   auto get = [&](int f, int x, int y) {
      return col[(size_t)kocipos[f]*S*S + (size_t)y*S + x] ;
   } ;
   int x, y ;
   // centre orbits: four of each colour
   for (int f=0; f<6; f++)
      cube.facerots[f] = 0 ;
   for (int i=0; i<cube.My; i++) {
      cube.rows[i].resize(cube.Mx) ;
      for (int j=0; j<cube.Mx; j++) {
         int cnt[6] = {0} ;
         for (int f=0; f<6; f++)
            for (int q=0; q<4; q++) {
               centrexy(N, i, j, q, x, y) ;
               uchar c = get(f, x, y) ;
               cube.rows[i][j].a[4*f+q] = c ;
               cnt[c]++ ;
            }
         for (int c=0; c<6; c++)
            if (cnt[c] != 4) {
               err = "centre orbit (" + to_string(i+1) + "," + to_string(j+1) + ") has " +
                     to_string(cnt[c]) + " of " + faceorder[c] + ", not 4" ;
               return 0 ;
            }
      }
   }
   // wings: each orbit has each of the 24 pieces once
   for (int k=0; k<(int)cube.wingorbits.size(); k++) {
      int seen[24] = {0} ;
      for (int slot=0; slot<24; slot++) {
         int ps = wingpartner[slot] ;
         borderxy(N, slot & 3, k+1, x, y) ;
         int c0 = get(slot >> 2, x, y) ;
         borderxy(N, ps & 3, N-k, x, y) ;
         int c1 = get(ps >> 2, x, y) ;
         int p = 0 ;
         while (p < 24 && (wingcols[p][0] != c0 || wingcols[p][1] != c1))
            p++ ;
         if (p == 24) {
            err = string("wing orbit ") + to_string(k+1) + " has a " + faceorder[c0] + faceorder[c1] +
                  " wing, which no cube has" ;
            return 0 ;
         }
         cube.wingorbits[k].a[slot] = (uchar)p ;
         seen[p]++ ;
      }
      for (int p=0; p<24; p++)
         if (seen[p] != 1) {
            err = string("wing orbit ") + to_string(k+1) + " has " + to_string(seen[p]) + " " +
                  faceorder[wingcols[p][0]] + faceorder[wingcols[p][1]] + " wings, not 1" ;
            return 0 ;
         }
   }
   // corners: each once; the twists sum to 0 mod 3
   int cperm[8] ;
   {
      int seen[8] = {0}, twist = 0 ;
      for (int i=0; i<8; i++) {
         int c[3] ;
         for (int t=0; t<3; t++) {
            int v = threecube::cor[i][t] ;
            cornerxy(N, v & 3, x, y) ;
            c[t] = get(v >> 2, x, y) ;
            cube.tc.corner.a[v] = (uchar)c[t] ;
         }
         int found = -1, rot = 0 ;
         for (int j=0; j<8 && found < 0; j++)
            for (int r=0; r<3 && found < 0; r++)
               if ((threecube::cor[j][r] >> 2) == c[0] && (threecube::cor[j][(r+1)%3] >> 2) == c[1] &&
                   (threecube::cor[j][(r+2)%3] >> 2) == c[2]) {
                  found = j ;
                  rot = r ;
               }
         if (found < 0) {
            err = string("corner ") + faceorder[c[0]] + faceorder[c[1]] + faceorder[c[2]] + " is not a cube corner" ;
            return 0 ;
         }
         seen[found]++ ;
         twist += rot ;
         cperm[i] = found ;
      }
      for (int j=0; j<8; j++)
         if (seen[j] != 1) {
            err = "the corners are not all different" ;
            return 0 ;
         }
      if (twist % 3) {
         err = "the corners are twisted (no sequence of moves can do that)" ;
         return 0 ;
      }
   }
   // middle edges (odd cubes): each once; the flips sum to 0 mod 2
   for (int i=0; i<24; i++)
      cube.tc.edge.a[i] = 0 ;
   if (N & 1) {
      const int half = (N+1)/2 ;
      int seen[12] = {0}, flip = 0, eperm[12] ;
      for (int i=0; i<12; i++) {
         int c[2] ;
         for (int t=0; t<2; t++) {
            int v = threecube::edg[i][t] ;
            borderxy(N, v & 3, half, x, y) ;
            c[t] = get(v >> 2, x, y) ;
            cube.tc.edge.a[v] = (uchar)c[t] ;
         }
         int found = -1, rot = 0 ;
         for (int j=0; j<12 && found < 0; j++)
            for (int r=0; r<2 && found < 0; r++)
               if ((threecube::edg[j][r] >> 2) == c[0] && (threecube::edg[j][1-r] >> 2) == c[1]) {
                  found = j ;
                  rot = r ;
               }
         if (found < 0) {
            err = string("middle edge ") + faceorder[c[0]] + faceorder[c[1]] + " is not a cube edge" ;
            return 0 ;
         }
         seen[found]++ ;
         flip += rot ;
         eperm[i] = found ;
      }
      for (int j=0; j<12; j++)
         if (seen[j] != 1) {
            err = "the middle edges are not all different" ;
            return 0 ;
         }
      if (flip & 1) {
         err = "a middle edge is flipped (no sequence of moves can do that)" ;
         return 0 ;
      }
      // centre centres: an orientation of the cube (a corner's three faces
      // keep their cyclic order, and opposite faces stay opposite)
      for (int f=0; f<6; f++) {
         cube.ccs[f] = get(f, half, half) ;
         cube.tc.centercenters[f] = cube.ccs[f] ;
      }
      bool ok = true ;
      for (int f=0; f<6; f++)
         if (cube.ccs[oppface[f]] != oppface[cube.ccs[f]])
            ok = false ;
      if (ok) {  // U L B (corner 0's faces) must map to some corner's faces in cyclic order
         int a = cube.ccs[threecube::cor[0][0]>>2], b = cube.ccs[threecube::cor[0][1]>>2],
             c = cube.ccs[threecube::cor[0][2]>>2] ;
         ok = false ;
         for (int j=0; j<8; j++)
            for (int r=0; r<3; r++)
               if ((threecube::cor[j][r] >> 2) == a && (threecube::cor[j][(r+1)%3] >> 2) == b &&
                   (threecube::cor[j][(r+2)%3] >> 2) == c)
                  ok = true ;
      }
      if (!ok) {
         err = "the centre centres are not an orientation of the cube" ;
         return 0 ;
      }
      // a face turn is a 4-cycle of corners and one of middle edges, a middle
      // slice turn one of middle edges and one of centre centres
      int cc[6] ;
      for (int f=0; f<6; f++)
         cc[f] = cube.ccs[f] ;
      if (parity(cperm, 8) ^ parity(eperm, 12) ^ parity(cc, 6)) {
         err = "two pieces are swapped (no sequence of moves can do that)" ;
         return 0 ;
      }
   } else
      for (int f=0; f<6; f++) {
         cube.ccs[f] = f ;
         cube.tc.centercenters[f] = f ;
      }
   return 1 ;
}
