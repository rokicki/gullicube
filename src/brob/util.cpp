#include "util.h"
#include <iostream>
#include <cstdlib>

std::mt19937_64 rng(std::random_device{}()) ;

void error(const char *s) {
   std::cerr << s << std::endl ;
   std::exit(10) ;
}
