# gullicube: N x N x N cube solver.  `make` builds ./gullicube, which reads its
# *.cls data files from data/ beside the executable (as named on the command line).
CXX      ?= c++
CXXFLAGS ?= -O3 -std=c++20
# DEFS: extra defines, e.g. `make clean; make DEFS=-DGULLI_PHASESTAT=1` for the
# per-phase instruction / cycle counters (run with PHASESTAT=1)
DEFS     ?=
CPPFLAGS += -Isrc/brob -Isrc/full $(DEFS)
LDLIBS   += -lpthread

BROB = src/brob/xcube.o src/brob/cornersolver.o src/brob/wingsolver.o src/brob/util.o src/brob/facelets.o
OBJS = src/full/gullicube.o src/full/brobglue.o $(BROB)
HDRS = $(wildcard src/full/*.h src/pair/*.h src/brob/*.h)

gullicube: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDLIBS)

src/full/%.o: src/full/%.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -Wall -c -o $@ $<

# brobdicube's sources, unmodified (warnings off)
src/brob/%.o: src/brob/%.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -w -c -o $@ $<

clean:
	rm -f gullicube $(OBJS)

.PHONY: clean
