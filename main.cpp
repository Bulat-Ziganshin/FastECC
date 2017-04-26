/// Implementation of the Number-Theoretical Transform in GF(P) Galois field
#include <iostream>
#include <algorithm>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <cassert>
#include <algorithm>
#include <utility>
#include <functional>
#include <vector>
#include <memory>

#include "wall_clock_timer.h"
#include "LargePages.cpp"
#include "GF(p).cpp"
#include "NTT.cpp"


/***********************************************************************************************************************
*** Testing/benchmarking routines **************************************************************************************
************************************************************************************************************************/

template <typename T, T P>
void Test_GF_Inv()
{
    int cnt = 0;
    for (T i=1; i<P; i++)
    {
        if (i%(1<<20)==0)  std::cout << std::hex << "\r0x" << i << "...";
        if (GF_Mul<T,P>(i, GF_Inv<T,P>(i)) != 1)
        {
            std::cout << i << "\n";
            if (++cnt==10) break;
        }
    }
}


// Find first few primary roots of 1 of power N
template <typename T, T P>
void FindRoot (T N)
{
    int cnt = 0;
    for (T i=2; i<P; i++)
    {
        if (i<256 || (i%(1024*1024))==0)
            std::cout << "\r" << i << "**" << std::hex << N << std::dec << "...";
        T q = GF_Pow<T,P> (i,N);
        if (q==1)
        {

            if (P==0x10001 || P==0xFFFFFFFF || P==0xFFFFFFFFFFFFFFFF) {
                if (1 == GF_Pow<T,P> (i,N/2))  goto next;
            } else {
                assert (P==0xFFF00001);  // other P aren' supported
                if (1 == GF_Pow<T,P> (i,N/2) ||
                    1 == GF_Pow<T,P> (i,N/3) ||
                    1 == GF_Pow<T,P> (i,N/5) ||
                    1 == GF_Pow<T,P> (i,N/7) ||
                    1 == GF_Pow<T,P> (i,N/13) ||
                    0)
                    goto next;
            }
/*


            for (size_t n=2; n<N; n++)
            {
                if (1 == GF_Pow<T,P> (i,n))
                    goto next;
            }

            T nn = N;
            for (size_t n=2; nn>1 && n<nn; n++)
            {
                while (nn>1 && nn%n==0)
                {
                    nn /= n;
                    if (1 == GF_Pow<T,P> (i,N/n))
                        goto next;
                }
            }
*/
            std::cout << i << "\n";
            if (++cnt==10) break;
        }
        next:;
    }
}


// Test the GF_Mul correctness
template <typename T, T P>
void Test_GF_Mul()
{
    int n = 0;
    for (T i=P-1; i>0; i--)
    {
        if (i%0x1000==0)  std::cout << std::hex << "\r0x" << i << "...";
        for (T j=P-1; j>=i; j--)
        {
            using DoubleT = typename Double<T>::T;
            auto a = (DoubleT(i)*j) % P;
            auto b = GF_Mul<T,P> (i,j);
            if (P==0xFFFFFFFF && b==P)  b=0;
            if (a != b)
            {
                std::cout << std::hex << "\r" << i << "*" << j << "=" << a << " != " << b << "\n" ;
                if (++n>10) return;
            }
        }
    }
}
template <> void Test_GF_Mul<uint64_t,0xFFFFFFFFFFFFFFFF>()
{
    printf("Test_GF_Mul<uint64_t>: unsupported\n");
}



// Print dividers count & density
template <typename T, T P>
void DividersDensity()
{
    std::vector<bool>  a(P, false);

    if (P == 0xFFF00001)
    {
        // fast algo manually optimized for 0xFFF00001
        for (T i=1; i<=9; i*=3)
        for (T k=1; k<=5; k*=5)
        for (T l=1; l<=7; l*=7)
        for (T m=1; m<=13; m*=13)
        for (T n=1; n<=1<<20; n*=2)
            a[i*k*l*m*n] = true;
    } else
    {
        // slow generic algo
        a[P-1] = true;
        for (T i=1; i<=P/2; i++)
            if ((P-1) % i  ==  0)
                a[i] = true;
    }

    // Compute dividers count & density
    T count = 0;  double prod = 1,  last = 0;
    T low = 1,  high = P;
    for (T i=low; i<high; i++)
        if (a[i]) {
            count++;
            if (last)  prod *= i / last;
            last = i+1;
        }

    // Density is the average "distance" to the next largest divider.
    // Density of 1.03 means that, at average, FFT order will be 3% higher than the block count.
    // Also you may want to check density for paricular range, f.e. 1000..10000 and/or without some P-1 factors.
    printf("%.0lf..%.0lf:  dividers count: %.0lf   density: %lf  ", low*1.0, high*1.0, count*1.0, pow(prod, 1/(2.0*(count-1))));
}


// Butterfly operation
template <typename T, T P>
void Butterfly (T* a, T* b, int TIMES, int SIZE, T root)
{
    for (int n=0; n<TIMES; n++) {
        for (int k=0; k<SIZE; k++) {                     // cycle over SIZE elements of the single block
            T temp = GF_Mul<T,P> (root, b[k]);
            b[k] = GF_Sub<T,P> (a[k], temp);
            a[k] = GF_Add<T,P> (a[k], temp);
        }
    }
}


// Benchmark speed of the Butterfly operation, processing 10 GB data == 500 MB processed with NTT<2**20>
template <typename T, T P>
int BenchButterfly()
{
    T x=0;
    #pragma omp parallel for
    for (int n=0; n<4096/sizeof(T); n++)
    {
        const int sz = 1280;
        T a[sz], b[sz];
        for (int i=0; i<sz; i++)
            a[i] = i*7+1, b[i] = i*15+8;
        Butterfly<T,P> (a, b, 1024, sz, 1557);
        x += a[0];
    }
    return x?1:0;
}


/***********************************************************************************************************************
*** The benchmarking driver ********************************************************************************************
************************************************************************************************************************/

// Return hash of the data
template <typename T>
uint32_t hash (T** data, size_t N, size_t SIZE)
{
    uint32_t hash = 314159253;
    for (size_t i=0; i<N; i++) {
        uint32_t* ptr = (uint32_t*) data[i];
        for (size_t k=0; k<SIZE*sizeof(T)/sizeof(uint32_t); k++)
            hash = (hash+ptr[k])*123456791 + (hash>>17);
    }
    return hash;
}


// Benchmark and verify two NTT implementations: Rec_NTT() & MFA_NTT(), compare results to definitive Slow_NTT()
template <typename T, T P>
void BenchNTT (bool RunOld, bool RunCanonical, size_t N, size_t SIZE, const char* P_str)
{
    bool RunNTT3  =  (N%3 == 0);
    bool RunNTT6  =  (N%6 == 0);
    bool RunNTT9  =  (N%9 == 0);

    T *data0 = VAlloc<T> (uint64_t(N)*SIZE);
    if (data0==0)  {printf("Can't alloc %.0lf MiB of memory!\n", (N/1048576.0)*SIZE*sizeof(T)); return;}

    for (size_t i=0; i<N*SIZE; i++)
        data0[i] = i%P;

    T **data = new T* [N];      // pointers to blocks
    for (size_t i=0; i<N; i++)
        data[i] = data0 + i*SIZE;

    uint32_t hash0 = hash(data, N, SIZE);    // hash of original data

    char title[99];
    int divider = RunOld? N : RunNTT9? 9 : RunNTT6? 6 : RunNTT3? 3 : N;
    sprintf (title, "NTT%d<%d*%.0lf,%.0lf,P=%s>", divider, divider, N*1.0/divider, SIZE*1.0*sizeof(T), P_str);
    for (int i=64; i--; )
        if (T(1)<<i == N)
            sprintf (title, "%s<2^%d,%.0lf,P=%s>", RunCanonical?"Slow_NTT":RunOld?"Rec_NTT":"MFA_NTT", i, SIZE*1.0*sizeof(T), P_str);

    double processed_size = (P==0x10001? 0.5:1.0) * N*SIZE*sizeof(T);   // In my GF(0x10001) implementation 4-byte value represents only 2 bytes of real data

         if (RunOld)       time_it (processed_size, title, [&]{Rec_NTT <T,P> (data, N, SIZE, false);});
    else if (RunNTT9)      time_it (processed_size, title, [&]{NTT9<T,P,false> (data, N/divider, SIZE);});
    else if (RunNTT6)      time_it (processed_size, title, [&]{NTT6<T,P,false> (data, N/divider, SIZE);});
    else if (RunNTT3)      time_it (processed_size, title, [&]{NTT3<T,P,false> (data, N/divider, SIZE);});
    else if (RunCanonical) time_it (processed_size, title, [&]{Slow_NTT<T,P> (data0,N, SIZE, false);});
    else                   time_it (processed_size, title, [&]{MFA_NTT <T,P> (data, N, SIZE, false);});

    // Pack results into 0..P-1 range
    for (size_t i=0; i<N*SIZE; i++)
        data0[i] = GF_Normalize<T,P> (data0[i]);

    uint32_t hash1 = hash(data, N, SIZE);    // hash after NTT

    // Inverse NTT
         if (RunOld)       Rec_NTT <T,P> (data, N, SIZE, true);
    else if (RunNTT9)      NTT9<T,P,true>(data, N/9, SIZE);
    else if (RunNTT6)      NTT6<T,P,true>(data, N/6, SIZE);
    else if (RunNTT3)      NTT3<T,P,true>(data, N/3, SIZE);
    else if (RunCanonical) Slow_NTT<T,P> (data0,N, SIZE, true);
    else                   MFA_NTT <T,P> (data, N, SIZE, true);

    // Normalize the result by dividing by N and pack results into 0..P-1 range
    T inv_N = GF_Inv<T,P>(divider);
    for (size_t i=0; i<N*SIZE; i++)
        data0[i] = GF_Normalize<T,P> (GF_Mul<T,P> (data0[i], inv_N));

    // Now we should have exactly the input data
    uint32_t hash2 = hash(data, N, SIZE);    // hash after NTT+iNTT
    if (hash2 == hash0) {
        if (verbose)  printf("Verified!  Original %.0lf,  after NTT: %.0lf\n", double(hash0), double(hash1));
    } else {
        printf("Checksum mismatch: original %.0lf,  after NTT: %.0lf,  after NTT+iNTT %.0lf\n", double(hash0), double(hash1), double(hash2));
    }
}


// Parse cmdline and invoke appropriate benchmark/test routine
template <typename T, T P>
void Code (int argc, char **argv, const char* P_str)
{
    if (argc>=2 && argv[1][0]=='.') {
        argv[1]++;
        verbose = false;
    }
    char opt  =  (argc>=2?  argv[1][0] : ' ');
    if (opt=='i')  {Test_GF_Inv<T,P>();  return;}
    if (opt=='m')  {Test_GF_Mul<T,P>();  return;}
    if (opt=='r')  {FindRoot<T,P>(P==0xFFFFFFFFFFFFFFFF?(uint64_t(65536)*2*5*17449):P==0xFFFFFFFF?65536:P-1);  printf ("GF_Root %s\n", GF_Root<T,P>(2)==P-1? "OK": "failed");  return;}
    if (opt=='d')  {DividersDensity<T,P>();  return;}
    if (opt=='b')  {time_it (10LL<<30, "Butterfly", [&]{BenchButterfly<T,P>();});  return;}

    size_t N = 1<<19;   // NTT order
    size_t SIZE = 2052; // Block size, in bytes
                        // 1 GB total

    if (argc>=3)  N = 1<<atoi(argv[2]);
    if (argc>=4)  SIZE = atoi(argv[3]);

    assert(N<P);  // Too long NTT for the such small P
    BenchNTT<T,P> (opt=='o', opt=='s', N, SIZE/sizeof(T), P_str);
}


// Deal with argv[1] prefix:
//   '.': quiet mode (on success, print only benchmark results)
//   '=': switch to P=0x10001
//   '-': switch to P=2^32-1 (not a primary number!)
//   '+': switch to P=2^64-1 (not a primary number!)
int main (int argc, char **argv)
{
    InitLargePages();
    if (argc>=2 && argv[1][0]=='.') {
        argv[1]++;
        verbose = false;
    }
    if (argc>=2 && argv[1][0]=='=') {
        argv[1]++;
        Code <uint32_t,0x10001> (argc, argv, "65537");
    } else if (argc>=2 && argv[1][0]=='-') {
        argv[1]++;
        Code <uint32_t,0xFFFFFFFF> (argc, argv, "2^32-1");
    } else if (argc>=2 && argv[1][0]=='+') {
#ifdef MY_CPU_64BIT
        argv[1]++;
        Code <uint64_t,0xFFFFFFFFFFFFFFFF> (argc, argv, "2^64-1");
#else         
        printf("Computations modulo 2^64-1 are supported only in 64-bit program versions\n");
#endif
    } else {
        Code <uint32_t,0xFFF00001> (argc, argv, "0xFFF00001");
    }
    return 0;
}

/* to do:
"b N SIZE" in cmdline
replace "(res>X)*P" in GF_Sub with bit arithmetics (compiler can do it itself?)
MS GF_Mul64 should became faster with the same algo as GCC one
IterativeNTT_Steps: optional extra twiddle factors in the last cycle so we can avoid them in MFA_NTT
template<typename GF> {operators +-*^/; static constexpr const GF root3, root3_2...;}
try to use "double" for GF(p)&Mod(p) operations, it may be faster
std::async: http://www.numberworld.org/y-cruncher/guides/multithreading.html

Order Multiplications
2     0
3     2
4     1
5     6
6     4
7     14
8     5
9     16
12    11
13    34
16    17
32    49 = 2*17 + 32/2-1 (MFA)
36    73 = 4*16 + 9*1
64    129
65    248 = 5*34+13*6
128   321
130   496
256   769
260   1057
512   1693
520   2309 = 8*248+65*5
585   3272 = 9*248+65*16
1024  3897
1040  5073 = 16*248+65*17
*/
