// Minimal NxNxN sticker simulator, just enough to turn a move sequence into its
// permutation on the (2,3)/(3,2) oblique centre pair and to check that it
// leaves everything else (except the diagonal orbits) alone.
//
// Move notation is twsearch / PuzzleGeometry: "F", "2F'", "3D2" -- a single
// layer at the given depth from the named face.  Metric is one move per token.
#pragma once
#include <array>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <functional>

// Face order U L F R B D (Speffz) = colour index.
static const char FACECH[] = "ULFRBD";
// outward normals
static const int NORMAL[6][3] = {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1},
                                 {1, 0, 0}, {0, 0, -1}, {0, -1, 0}};

inline int face_of_char(char c) {
  for (int f = 0; f < 6; f++)
    if (FACECH[f] == c) return f;
  return -1;
}
inline int face_of_normal(int x, int y, int z) {
  for (int f = 0; f < 6; f++)
    if (NORMAL[f][0] == x && NORMAL[f][1] == y && NORMAL[f][2] == z) return f;
  return -1;
}

struct Move {
  uint8_t face, depth, amt;  // amt 1,2,3 (quarter turns clockwise)
};

inline bool parse_move(const std::string &tok, Move &m) {
  size_t i = 0;
  int depth = 1;
  if (i < tok.size() && isdigit((unsigned char)tok[i])) {
    depth = 0;
    while (i < tok.size() && isdigit((unsigned char)tok[i])) depth = depth * 10 + (tok[i++] - '0');
  }
  if (i >= tok.size()) return false;
  int f = face_of_char(tok[i++]);
  if (f < 0) return false;
  int amt = 1;
  if (i < tok.size()) {
    if (tok[i] == '2') amt = 2;
    else if (tok[i] == '\'' || tok[i] == '3') amt = 3;
    else if (tok[i] == '1') amt = 1;
    else return false;
  }
  m.face = f;
  m.depth = depth;
  m.amt = amt;
  return true;
}

inline std::string move_str(const Move &m) {
  std::string s;
  if (m.depth > 1) s += std::to_string(m.depth);
  s += FACECH[m.face];
  if (m.amt == 2) s += "2";
  if (m.amt == 3) s += "'";
  return s;
}

inline std::vector<Move> parse_alg(const std::string &line) {
  std::vector<Move> r;
  size_t i = 0;
  while (i < line.size()) {
    while (i < line.size() && isspace((unsigned char)line[i])) i++;
    size_t j = i;
    while (j < line.size() && !isspace((unsigned char)line[j])) j++;
    if (j > i) {
      Move m;
      if (!parse_move(line.substr(i, j - i), m)) throw std::runtime_error("bad move " + line.substr(i, j - i));
      r.push_back(m);
    }
    i = j;
  }
  return r;
}

struct Cube {
  int n;
  // stickers: doubled coordinates of sticker centre (cubie centre + normal)
  std::vector<std::array<int, 3>> pos;
  std::vector<int> face;  // which face the sticker is on
  std::vector<int> row, col;
  // movePerm[face][depth-1] = where each sticker goes under one clockwise quarter turn
  std::vector<std::vector<std::vector<int>>> movePerm;
  std::vector<int> orbit;  // orbit id per sticker
  int norbits = 0;

  int find(const std::array<int, 3> &p) const {
    for (size_t i = 0; i < pos.size(); i++)
      if (pos[i] == p) return (int)i;
    return -1;
  }

  explicit Cube(int n_) : n(n_) {
    // cubie coordinates are in {-(n-1), ..., n-1} step 2; sticker = cubie + normal
    for (int f = 0; f < 6; f++) {
      for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++) {
          // build a face-local frame: u (row direction), v (col direction)
          int a = -(n - 1) + 2 * r, b = -(n - 1) + 2 * c;
          std::array<int, 3> p;
          const int *N = NORMAL[f];
          // pick two axes perpendicular to normal
          int ax = N[0] ? 0 : (N[1] ? 1 : 2);
          int o1 = (ax + 1) % 3, o2 = (ax + 2) % 3;
          p[ax] = N[ax] * n;
          p[o1] = a;
          p[o2] = b;
          pos.push_back(p);
          face.push_back(f);
          row.push_back(r);
          col.push_back(c);
        }
    }
    int S = pos.size();
    // lookup table
    int W = 2 * n + 1;
    std::vector<int> lut(W * W * W, -1);
    auto key = [&](const std::array<int, 3> &p) { return ((p[0] + n) * W + (p[1] + n)) * W + (p[2] + n); };
    for (int i = 0; i < S; i++) lut[key(pos[i])] = i;
    movePerm.assign(6, std::vector<std::vector<int>>(n));
    for (int f = 0; f < 6; f++) {
      const int *N = NORMAL[f];
      int ax = N[0] ? 0 : (N[1] ? 1 : 2);
      int sgn = N[ax];
      for (int d = 1; d <= n; d++) {
        // layer centre coordinate along ax (cubie coords), measured from face f
        int layer = sgn * (n - 1 - 2 * (d - 1));
        std::vector<int> perm(S);
        for (int i = 0; i < S; i++) {
          std::array<int, 3> p = pos[i];
          // cubie coordinate along ax: sticker minus normal if sticker is on a face perpendicular to ax
          int cub = p[ax];
          if (std::abs(p[ax]) == n) cub = p[ax] - (p[ax] > 0 ? 1 : -1);
          if (cub != layer) {
            perm[i] = i;
            continue;
          }
          // clockwise seen from face f = rotation by -90 deg about outward normal
          // rotation +90 about axis e_ax: (o1,o2) -> (-o2, o1) with o1=ax+1,o2=ax+2
          int o1 = (ax + 1) % 3, o2 = (ax + 2) % 3;
          std::array<int, 3> q = p;
          if (sgn > 0) {  // -90 about +e: (o1,o2)->(o2,-o1)
            q[o1] = p[o2];
            q[o2] = -p[o1];
          } else {  // -90 about -e = +90 about +e
            q[o1] = -p[o2];
            q[o2] = p[o1];
          }
          int j = lut[key(q)];
          if (j < 0) abort();
          perm[i] = j;
        }
        movePerm[f][d - 1] = perm;
      }
    }
    // orbits under all moves
    std::vector<int> uf(S);
    for (int i = 0; i < S; i++) uf[i] = i;
    std::function<int(int)> fd = [&](int x) { return uf[x] == x ? x : uf[x] = fd(uf[x]); };
    for (auto &fv : movePerm)
      for (auto &p : fv)
        for (int i = 0; i < S; i++) uf[fd(i)] = fd(p[i]);
    orbit.assign(S, -1);
    std::vector<int> id(S, -1);
    for (int i = 0; i < S; i++) {
      int r = fd(i);
      if (id[r] < 0) id[r] = norbits++;
      orbit[i] = id[r];
    }
  }
  // depth of sticker from the nearest edge in row / col direction (1-based)
  int drow(int i) const { return std::min(row[i], n - 1 - row[i]) + 1; }
  int dcol(int i) const { return std::min(col[i], n - 1 - col[i]) + 1; }

  // apply a sequence: returns loc[s] = where the sticker starting at s ends up
  std::vector<int> apply(const std::vector<Move> &alg) const {
    int S = pos.size();
    std::vector<int> loc(S);
    for (int i = 0; i < S; i++) loc[i] = i;
    for (auto &m : alg) {
      const auto &p = movePerm[m.face][m.depth - 1];
      for (int k = 0; k < m.amt; k++)
        for (int i = 0; i < S; i++) loc[i] = p[loc[i]];
    }
    return loc;
  }
};
