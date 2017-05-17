CXX ?= g++

# add -m32 if trying to build 32 bit binaries on 64 bit machines,
# on Linux will also need gcc-multilib installed
CXXFLAGS ?= -std=c++1y -O3 -s -fopenmp

# f.e. you can define suffix to "64g6" to distinguish 64-bit executables produced by GCC6
SUFFIX ?=

# flags to make vectorized SSE2/AVX2 builds
# (without -DSIMD the sources will use unvectorizable code path even on x64)
SSE2_FLAGS ?= -msse2 -DSIMD=SSE2
AVX2_FLAGS ?= -mavx2 -DSIMD=AVX2

SRCFILES = Makefile GF(p).cpp LargePages.cpp ntt.cpp SIMD.h wall_clock_timer.h
EXEFILES = ntt$(SUFFIX) ntt$(SUFFIX)-sse2 ntt$(SUFFIX)-avx2 rs$(SUFFIX) rs$(SUFFIX)-sse2 rs$(SUFFIX)-avx2 prime

all: $(EXEFILES)
.PHONY : all

ntt$(SUFFIX): main.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $<

ntt$(SUFFIX)-sse2: main.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $< $(SSE2_FLAGS)

ntt$(SUFFIX)-avx2: main.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $< $(AVX2_FLAGS)

rs$(SUFFIX): RS.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $<

rs$(SUFFIX)-sse2: RS.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $< $(SSE2_FLAGS)

rs$(SUFFIX)-avx2: RS.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $< $(AVX2_FLAGS)

prime: prime.cpp
	$(CXX) -O2 -s -o $@ $<

.PHONY: clean
clean:
	rm -f $(EXEFILES)
