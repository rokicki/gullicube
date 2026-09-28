#pragma once
#include "xcube.h"

extern std::vector<mv> curblock ;
extern ll wingmoves ;
extern ll blocks ;
extern ll runningmovecount ;
extern int writemoves ;

void sendblock() ;
void build_wing_cycle_index() ;
void loadwingalgs() ;
ll do_wing_cycles(xcube &cube) ;
