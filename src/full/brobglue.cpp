// Globals that brobdicube's xcube.cpp / cornersolver.cpp / wingsolver.cpp
// expect from brobdicube.cpp, so those files can be linked unmodified.
#include "xcube.h"
#include "wingsolver.h"
#include <chrono>

int alwaysinmem = 1;
ll cancellations = 0;
ll runningmovecount = 0;
int writemoves = 0;
std::vector<mv> curblock;
ll wingmoves = 0;
ll blocks = 0;
void sendblock() { curblock.clear(); }

static double walltime() {
  using namespace std::chrono;
  return duration_cast<duration<double>>(steady_clock::now().time_since_epoch()).count();
}
double duration() {
  static double last = walltime();
  double now = walltime(), r = now - last;
  last = now;
  return r;
}
