// Minimum over a beam table's slot scores, kept up to date in O(log size) per
// slot write: the beams' admission threshold (the weakest slot) without a
// scan of the whole table after every admission.  Empty slots hold INT32_MIN,
// so the minimum is INT32_MIN until every slot is filled.
#pragma once
#include <algorithm>
#include <climits>
#include <vector>

struct MinTree {
  std::vector<int> t;
  size_t m = 1;
  void reset(size_t n) {
    m = 1;
    while (m < n) m <<= 1;
    t.assign(2 * m, INT_MAX);
    for (size_t i = 0; i < n; i++) t[m + i] = INT_MIN;
    for (size_t i = m - 1; i >= 1; i--) t[i] = std::min(t[2 * i], t[2 * i + 1]);
  }
  void set(size_t i, int v) {
    i += m;
    t[i] = v;
    for (i >>= 1; i >= 1; i >>= 1) {
      const int x = std::min(t[2 * i], t[2 * i + 1]);
      if (t[i] == x) break;
      t[i] = x;
    }
  }
  int min() const { return t[1]; }
};
