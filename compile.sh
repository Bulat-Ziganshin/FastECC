g++ -std=c++14 -m64 -O3 main.cpp -ontt64g-avx2  -static -s -fopenmp -mavx2   -DSIMD=AVX2
g++ -std=c++14 -m64 -O3 main.cpp -ontt64g-sse2  -static -s -fopenmp          -DSIMD=SSE2
g++ -std=c++14 -m64 -O3 main.cpp -ontt64g       -static -s -fopenmp

g++ -std=c++14 -m32 -O3 main.cpp -ontt32g-avx2  -static -s -fopenmp -mavx2   -DSIMD=AVX2
g++ -std=c++14 -m32 -O3 main.cpp -ontt32g-sse2  -static -s -fopenmp -msse2   -DSIMD=SSE2
g++ -std=c++14 -m32 -O3 main.cpp -ontt32g       -static -s -fopenmp -mmmx

g++ -std=c++14 -m64 -O3 RS.cpp -ors64g-avx2     -static -s -fopenmp -mavx2   -DSIMD=AVX2
g++ -std=c++14 -m64 -O3 RS.cpp -ors64g-sse2     -static -s -fopenmp          -DSIMD=SSE2
g++ -std=c++14 -m64 -O3 RS.cpp -ors64g          -static -s -fopenmp

g++ -std=c++14 -m32 -O3 RS.cpp -ors32g-avx2     -static -s -fopenmp -mavx2   -DSIMD=AVX2
g++ -std=c++14 -m32 -O3 RS.cpp -ors32g-sse2     -static -s -fopenmp -msse2   -DSIMD=SSE2
g++ -std=c++14 -m32 -O3 RS.cpp -ors32g          -static -s -fopenmp -mmmx
