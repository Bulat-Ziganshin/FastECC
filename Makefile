CXX=g++

#add -m32 if trying to build 32 bit binaries on 64 bit machines,
#on Linux will also need gcc-multilib installed
CXXFLAGS=-std=c++1y -O3 -s -fopenmp

#flags to make vectorized SSE2/AVX2 builds
#(without -DSIMD the sources will use unvectorized code path even on x64)
SSE2_FLAGS=-msse2 -DSIMD=SSE2
AVX2_FLAGS=-mavx2 -DSIMD=AVX2

SRCFILES=Makefile GF(p).cpp LargePages.cpp ntt.cpp SIMD.h wall_clock_timer.h
EXEFILES=nttg nttg-sse2 nttg-avx2 rsg rsg-sse2 rsg-avx2 prime

all: $(EXEFILES)

nttg: main.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $<

nttg-sse2: main.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $< $(SSE2_FLAGS)

nttg-avx2: main.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $< $(AVX2_FLAGS)

rsg: RS.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $<

rsg-sse2: RS.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $< $(SSE2_FLAGS)

rsg-avx2: RS.cpp $(SRCFILES)
	$(CXX) $(CXXFLAGS) -o $@ $< $(AVX2_FLAGS)

prime: prime.cpp
	$(CXX) -O2 -s -o $@ $<

clean:
	rm -f $(EXEFILES)
