CXX=g++
#add -m32 if trying to build 32 bit binaries on 64 bit machines, will
#also need gcc-multilib installed
CFLAGS=-std=c++1y -O3 -s -fopenmp

nttg: nttg-avx2 nttg-sse2 rsg-avx2 rsg-sse2 rsg 
	$(CXX) $(CFLAGS) -o nttg      main.cpp

nttg-avx2: 
	$(CXX) $(CFLAGS) -o nttg-avx2 main.cpp -mavx2 -DSIMD=AVX2 

nttg-sse2:
	$(CXX) $(CFLAGS) -o nttg-sse2 main.cpp -msse2 -DSIMD=SSE2

rsg: 
	$(CXX) $(CFLAGS) -o rsg      RS.cpp

rsg-avx2: 
	$(CXX) $(CFLAGS) -o rsg-avx2 RS.cpp -mavx2 -DSIMD=AVX2 

rsg-sse2:
	$(CXX) $(CFLAGS) -o rsg-sse2 RS.cpp -msse2 -DSIMD=SSE2

clean:
	rm -f nttg nttg-avx2 nttg-sse2 rsg rsg-avx2 rsg-sse2
