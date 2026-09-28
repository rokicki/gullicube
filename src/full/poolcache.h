// Binary cache for pair-solver pools: the first load of FILE parses the text
// and writes FILE.bin; later loads read the binary (much faster for the
// multi-million-entry endgame pools).  Cache is rebuilt if the text is newer.
#pragma once
#include "../pair/pool.h"
#include <chrono>
#include <filesystem>

namespace poolcache {
template <class T> static void wr(FILE *f, const std::vector<T> &v) {
  uint64_t n = v.size();
  fwrite(&n, 8, 1, f);
  if (n) fwrite(v.data(), sizeof(T), n, f);
}
template <class T> static bool rd(FILE *f, std::vector<T> &v) {
  uint64_t n;
  if (fread(&n, 8, 1, f) != 1) return false;
  v.resize(n);
  return !n || fread(v.data(), sizeof(T), n, f) == n;
}
static double mtime(const std::string &fn) {
  std::error_code ec;
  auto t = std::filesystem::last_write_time(fn, ec);
  if (ec) return -1;
  return (double)std::chrono::duration_cast<std::chrono::seconds>(t.time_since_epoch()).count();
}
}  // namespace poolcache

// append FILE's entries to p (text or cached binary), then return p
inline void loadPoolCached(Pool &p, const std::string &fn) {
  using namespace poolcache;
  std::string bin = fn + ".bin";
  Pool q;
  bool ok = false;
  if (mtime(bin) >= mtime(fn) && mtime(fn) >= 0) {
    FILE *f = fopen(bin.c_str(), "rb");
    if (f) {
      ok = rd(f, q.S) && rd(f, q.src) && rd(f, q.cost) && rd(f, q.support) && rd(f, q.diag) && rd(f, q.goff) &&
           rd(f, q.gdat) && rd(f, q.head) && rd(f, q.second);
      fclose(f);
      q.text.assign(q.cost.size(), std::string());
    }
  }
  if (!ok) {
    q = Pool();
    q.load(fn.c_str(), 99, false);
    FILE *f = fopen(bin.c_str(), "wb");
    if (f) {
      wr(f, q.S); wr(f, q.src); wr(f, q.cost); wr(f, q.support); wr(f, q.diag); wr(f, q.goff); wr(f, q.gdat);
      wr(f, q.head); wr(f, q.second);
      fclose(f);
    }
  }
  // append q to p
  size_t off = p.gdat.size();
  p.S.insert(p.S.end(), q.S.begin(), q.S.end());
  p.src.insert(p.src.end(), q.src.begin(), q.src.end());
  p.cost.insert(p.cost.end(), q.cost.begin(), q.cost.end());
  p.support.insert(p.support.end(), q.support.begin(), q.support.end());
  p.diag.insert(p.diag.end(), q.diag.begin(), q.diag.end());
  p.text.insert(p.text.end(), q.text.begin(), q.text.end());
  p.head.insert(p.head.end(), q.head.begin(), q.head.end());
  p.second.insert(p.second.end(), q.second.begin(), q.second.end());
  p.gdat.insert(p.gdat.end(), q.gdat.begin(), q.gdat.end());
  for (size_t i = 1; i < q.goff.size(); i++) p.goff.push_back(q.goff[i] + off);
}
