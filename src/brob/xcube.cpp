#include "xcube.h"
#include <iostream>
#include <string>
using namespace std ;

int rots[24][6] = {
   {0, 1, 2, 3, 4, 5}, {0, 2, 3, 4, 1, 5}, {2, 1, 5, 3, 0, 4},
} ;
const char *faceorder = "ULFRBD" ;
const int oppface[] = { 5, 3, 4, 1, 2, 0 } ;
vector<mv> moveset ;
vector<mv> prevmvs ;
vector<mv> *g_prep_recorder = nullptr ;

void coordtospeffz(int Mx, int My, int x, int y, int &ox, int &oy, int &oz) {
   assert(x >= 1) ;
   assert(x <= Mx + My) ;
   assert(y >= 1) ;
   assert(y <= Mx + My) ;
   if (y <= My && x <= Mx) {
      ox = x ;
      oy = y ;
      oz = 0 ;
   } else if (y <= Mx && x > Mx) {
      oy = Mx + My + 1 - x ;
      ox = y ;
      oz = 1 ;
   } else if (y > Mx && x > My) {
      ox = Mx + My + 1 - x ;
      oy = Mx + My + 1 - y ;
      oz = 2 ;
   } else if (y > My && x <= My) {
      ox = Mx + My + 1 - y ;
      oy = x ;
      oz = 3 ;
   } else {
      ox = 1 ;
      oy = 1 ;
      oz = 4 ;
   }
   assert(ox <= Mx) ;
   assert(oy <= My) ;
   assert(ox >= 1) ;
   assert(oy >= 1) ;
}

void border() {
   cout << "  " ;
}

void xcube::show() const {
   if (Mx+My > 1000) {
      cout << "Not showing a cube that big" << endl ;
      return ;
   }
   for (int y=1; y<=Mx+My; y++) {
      border() ;
      for (int x=1; x<=Mx+My; x++) {
         cout << " " ;
      }
      border() ;
      for (int x=1; x<=Mx+My; x++)
         cout << get(0, x, y) ;
      cout << endl ;
   }
   cout << endl ;
   for (int y=1; y<=Mx+My; y++) {
      for (int f=1; f<=4; f++) {
         border() ;
         for (int x=1; x<=Mx+My; x++)
            cout << get(f, x, y) ;
      }
      cout << endl ;
   }
   cout << endl ;
   for (int y=1; y<=Mx+My; y++) {
      border() ;
      for (int x=1; x<=Mx+My; x++) {
         cout << " " ;
      }
      border() ;
      for (int x=1; x<=Mx+My; x++)
         cout << get(5, x, y) ;
      cout << endl ;
   }
   cout << endl ;
}

mv parse_move(const char *& s) {
   mv r ;
   r.dep = 1 ;
   while (*s && *s <= ' ')
      s++ ;
   if ('0' <= *s && *s <= '9') {
      r.dep = 0 ;
      while ('0' <= *s && *s <= '9')
         r.dep = 10 * r.dep + *s++ - '0' ;
   }
   auto found = strchr(faceorder, *s) ;
   if (found == 0) {
      cerr << s << endl ;
      error("! could not parse move") ;
   }
   r.face = found - faceorder ;
   s++ ;
   r.twist = 1 ;
   if (*s && *s > ' ') {
      if ('1' <= *s && *s <= '3')
         r.twist = *s - '0' ;
      else if (*s == '\'')
         r.twist = 3 ;
      else {
         cerr << s << endl ;
         error("! could not parse move (2)") ;
      }
      s++ ;
   }
   while (*s && *s <= ' ')
      s++ ;
   return r ;
}

void readmoveset() {
   string s ;
   while (cin >> s) {
      const char *mstr = s.c_str() ;
      moveset.push_back(parse_move(mstr)) ;
   }
   cout << "Read move set; got " << moveset.size() << endl ;
}

vector<mv> parse_moves(const string &alg) {
   vector<mv> r ;
   const char *s = alg.c_str() ;
   while (*s) {
      while (*s && *s <= ' ')
         s++ ;
      if (*s == 0)
         break ;
      r.push_back(parse_move(s)) ;
   }
   return r ;
}

int trymerge(const mv &a, const mv &b, mv &res) {
   if (a.face == b.face && a.dep == b.dep) {
      int tw = (a.twist + b.twist) & 3 ;
      if (tw == 0)
         return 0 ;
      res = a ;
      res.twist = tw ;
      return 1 ;
   } else {
      return 2 ;
   }
}

void trycancellation(mv &m) {
   if (prevmvs.size() == 0) {
      prevmvs.push_back(m) ;
      return ;
   }
   if (m.dep > 1) {
      prevmvs.clear() ;
      return ;
   }
   mv merged ;
   int r = trymerge(prevmvs[prevmvs.size()-1], m, merged) ;
   if (r == 2) {
      prevmvs.push_back(m) ;
   } else if (r == 1) {
      cancellations++ ;
      prevmvs[prevmvs.size()-1] = merged ;
   } else {
      cancellations += 2 ;
      prevmvs.pop_back() ;
   }
}

void showmove(const struct mv &mv) {
   if (mv.dep > 1)
      cout << mv.dep ;
   cout << faceorder[mv.face] ;
   if (mv.twist == 2)
      cout << 2 ;
   else if (mv.twist == 3)
      cout << "'" ;
}

void random_move_scramble(cubesizetype N, ll length) {
   moveset.clear() ;
   int lastface = -1 ;
   cout << "SCRAMBLE:" ;
   for (ll i=0; i<length; i++) {
      mv m ;
      do { m.face = rng() % 6 ; } while (m.face == lastface) ;
      lastface = m.face ;
      m.dep   = 1 + rng() % ((N + 1) / 2) ;
      m.twist = 1 + rng() % 3 ;
      moveset.push_back(m) ;
      cout << " " ;
      showmove(m) ;
   }
   cout << endl ;
}

void domove(mv mv, xcube &cube, bool suppress_cancellations) {
   if (mv.dep > cube.Mx + cube.My + 3 - mv.dep) {
      mv.dep = cube.Mx + cube.My + 3 - mv.dep ;
      mv.face = oppface[mv.face] ;
      mv.twist = (4 - mv.twist) & 3 ;
   }
   if (g_prep_recorder)
      g_prep_recorder->push_back(mv) ;
   assert(mv.dep >= 1) ;
   assert(mv.dep <= cube.Mx + 1) ;
   if (!suppress_cancellations) trycancellation(mv) ;
   if (mv.dep == 1) {
      cube.tc.domove(mv) ;
      cube.facerots[mv.face] = (cube.facerots[mv.face] + 4 - mv.twist) & 3 ;
      for (int m=0; m<cube.My; m++){
         for (int i=0; i<mv.twist; i++) {
            if (mv.face == 0) {
               cube.wings_cycle(0,1,2,3, m);
               cube.wings_cycle(16,12,8,4, m);
            } else if (mv.face == 1) {
               cube.wings_cycle(4,5,6,7, m);
               cube.wings_cycle(3,11,23,17, m);
            } else if (mv.face == 2) {
               cube.wings_cycle(8,9,10,11, m);
               cube.wings_cycle(2,15,20,5, m);
            } else if (mv.face == 3) {
               cube.wings_cycle(12,13,14,15, m);
               cube.wings_cycle(1,19,21,9, m);
            } else if (mv.face == 4) {
               cube.wings_cycle(16,17,18,19, m);
               cube.wings_cycle(0,7,22,13, m);
            } else if (mv.face == 5) {
               cube.wings_cycle(20,21,22,23, m);
               cube.wings_cycle(10,14,18,6, m);
            }
         }
      }
      return ;
   }
   assert(mv.twist > 0) ;
   assert(mv.twist < 4) ;
   int dep = mv.dep - 2 ;
   if (dep >= cube.My) {
      cube.tc.domove({2, mv.face, mv.twist}) ;
      switch (mv.face * 3 + mv.twist) {
case 1:  cube.rot4ver2(1, 1, 2, 1, 3, 1, 4, 1, dep) ; break ;
case 2:  cube.rot2ver2(1, 1, 2, 1, 3, 1, 4, 1, dep) ; break ;
case 3:  cube.rot4ver2(4, 1, 3, 1, 2, 1, 1, 1, dep) ; break ;
case 4:  cube.rot4ver2(2, 0, 0, 0, 4, 2, 5, 0, dep) ; break ;
case 5:  cube.rot2ver2(2, 0, 5, 0, 4, 2, 0, 0, dep) ; break ;
case 6:  cube.rot4ver2(2, 0, 5, 0, 4, 2, 0, 0, dep) ; break ;
case 7:  cube.rot4ver2(0, 3, 1, 2, 5, 1, 3, 0, dep) ; break ;
case 8:  cube.rot2ver2(0, 3, 3, 0, 5, 1, 1, 2, dep) ; break ;
case 9:  cube.rot4ver2(0, 3, 3, 0, 5, 1, 1, 2, dep) ; break ;
case 10: cube.rot4ver2(2, 2, 5, 2, 4, 0, 0, 2, dep) ; break ;
case 11: cube.rot2ver2(2, 2, 5, 2, 4, 0, 0, 2, dep) ; break ;
case 12: cube.rot4ver2(2, 2, 0, 2, 4, 0, 5, 2, dep) ; break ;
case 13: cube.rot4ver2(0, 1, 3, 2, 5, 3, 1, 0, dep) ; break ;
case 14: cube.rot2ver2(0, 1, 3, 2, 5, 3, 1, 0, dep) ; break ;
case 15: cube.rot4ver2(0, 1, 1, 0, 5, 3, 3, 2, dep) ; break ;
case 16: cube.rot4ver2(4, 3, 3, 3, 2, 3, 1, 3, dep) ; break ;
case 17: cube.rot2ver2(1, 3, 2, 3, 3, 3, 4, 3, dep) ; break ;
case 18: cube.rot4ver2(4, 3, 1, 3, 2, 3, 3, 3, dep) ; break ;
      }
   } else {
      switch (mv.face * 3 + mv.twist) {
case 1:  cube.rot4hor(1, 0, 2, 0, 3, 0, 4, 0, dep) ; cube.rot4ver(1, 1, 2, 1, 3, 1, 4, 1, dep) ; cube.wings_cycle(5,17,13,9, dep) ; break ;
case 2:  cube.rot2hor(1, 0, 2, 0, 3, 0, 4, 0, dep) ; cube.rot2ver(1, 1, 2, 1, 3, 1, 4, 1, dep) ; cube.wings_cycle2(5,17,13,9, dep) ; break ;
case 3:  cube.rot4hor(4, 0, 3, 0, 2, 0, 1, 0, dep) ; cube.rot4ver(4, 1, 3, 1, 2, 1, 1, 1, dep) ; cube.wings_cycle(5,9,13,17, dep) ; break ;
case 4:  cube.rot4hor(2, 3, 0, 3, 4, 1, 5, 3, dep) ; cube.rot4ver(2, 0, 0, 0, 4, 2, 5, 0, dep) ; cube.wings_cycle(0,8,20,18, dep) ; break ;
case 5:  cube.rot2hor(2, 3, 5, 3, 4, 1, 0, 3, dep) ; cube.rot2ver(2, 0, 5, 0, 4, 2, 0, 0, dep) ; cube.wings_cycle2(0,8,20,18, dep) ; break ;
case 6:  cube.rot4hor(2, 3, 5, 3, 4, 1, 0, 3, dep) ; cube.rot4ver(2, 0, 5, 0, 4, 2, 0, 0, dep) ; cube.wings_cycle(0,18,20,8, dep) ; break ;
case 7:  cube.rot4hor(0, 2, 1, 1, 5, 0, 3, 3, dep) ; cube.rot4ver(0, 3, 1, 2, 5, 1, 3, 0, dep) ; cube.wings_cycle(3,12,21,6, dep) ; break ;
case 8:  cube.rot2hor(0, 2, 3, 3, 5, 0, 1, 1, dep) ; cube.rot2ver(0, 3, 3, 0, 5, 1, 1, 2, dep) ; cube.wings_cycle2(3,12,21,6, dep) ; break ;
case 9:  cube.rot4hor(0, 2, 3, 3, 5, 0, 1, 1, dep) ; cube.rot4ver(0, 3, 3, 0, 5, 1, 1, 2, dep) ; cube.wings_cycle(3,6,21,12, dep) ; break ;
case 10: cube.rot4hor(2, 1, 5, 1, 4, 3, 0, 1, dep) ; cube.rot4ver(2, 2, 5, 2, 4, 0, 0, 2, dep) ; cube.wings_cycle(2, 16, 22, 10, dep) ; break ;
case 11: cube.rot2hor(2, 1, 5, 1, 4, 3, 0, 1, dep) ; cube.rot2ver(2, 2, 5, 2, 4, 0, 0, 2, dep) ; cube.wings_cycle2(2, 16, 22, 10, dep) ; break ;
case 12: cube.rot4hor(2, 1, 0, 1, 4, 3, 5, 1, dep) ; cube.rot4ver(2, 2, 0, 2, 4, 0, 5, 2, dep) ; cube.wings_cycle(2, 10, 22, 16, dep) ; break ;
case 13: cube.rot4hor(0, 0, 3, 1, 5, 2, 1, 3, dep) ; cube.rot4ver(0, 1, 3, 2, 5, 3, 1, 0, dep) ; cube.wings_cycle(1, 4, 23, 14, dep) ; break ;
case 14: cube.rot2hor(0, 0, 3, 1, 5, 2, 1, 3, dep) ; cube.rot2ver(0, 1, 3, 2, 5, 3, 1, 0, dep) ; cube.wings_cycle2(1, 4, 23, 14, dep) ; break ;
case 15: cube.rot4hor(0, 0, 1, 3, 5, 2, 3, 1, dep) ; cube.rot4ver(0, 1, 1, 0, 5, 3, 3, 2, dep) ; cube.wings_cycle(1, 14, 23, 4, dep) ; break ;
case 16: cube.rot4hor(4, 2, 3, 2, 2, 2, 1, 2, dep) ; cube.rot4ver(4, 3, 3, 3, 2, 3, 1, 3, dep) ; cube.wings_cycle(7, 11, 15, 19, dep) ; break ;
case 17: cube.rot2hor(1, 2, 2, 2, 3, 2, 4, 2, dep) ; cube.rot2ver(1, 3, 2, 3, 3, 3, 4, 3, dep) ; cube.wings_cycle2(7, 19, 15, 11, dep) ; break ;
case 18: cube.rot4hor(4, 2, 1, 2, 2, 2, 3, 2, dep) ; cube.rot4ver(4, 3, 1, 3, 2, 3, 3, 3, dep) ; cube.wings_cycle(7, 19, 15, 11, dep) ; break ;
      }
   }
}

int encoderot(int *p) {
   return p[0] * 6 + p[1] ;
}

void expandrotations() {
   const int MAXENC = 6*6 ;
   char seen[MAXENC] ;
   if (rots[23][0] != 0 || rots[23][1] != 0)
      return ;
   for (int i=0; i<MAXENC; i++)
      seen[i] = 0 ;
   for (int i=0; i<3; i++)
      seen[encoderot(rots[i])] = 1 ;
   int qp = 3 ;
   int w[6] ;
   for (int qg=1; qg<24; qg++) {
      for (int g=1; g<3; g++) {
         for (int i=0; i<6; i++)
            w[i] = rots[qg][rots[g][i]] ;
         int ev = encoderot(w) ;
         if (seen[ev] == 0) {
            seen[ev] = 1 ;
            for (int i=0; i<6; i++)
               rots[qp][i] = w[i] ;
            qp++ ;
            if (qp > 24)
               error("! blew up expanding rotations") ;
         }
      }
   }
   if (qp != 24)
      error("! did not find the correct count of rotations") ;
}
