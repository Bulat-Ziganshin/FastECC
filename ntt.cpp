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

#if defined(_M_X64) || defined(_M_AMD64) || defined(__x86_64__)
#define MY_CPU_AMD64
#endif

#if defined(MY_CPU_AMD64) || defined(_M_IA64)
#define MY_CPU_64BIT
#endif

#ifndef __GNUC__
#define __restrict__
#endif


void time_it (int64_t size, const char* name, std::function<void()> Code)
{
    static int _ = (StartTimer(),0);
    double start = GetTimer(), KernelTime[2], UserTime[2];
    GetProcessKernelUserTimes (KernelTime, UserTime);
    Code();
    GetProcessKernelUserTimes (KernelTime+1, UserTime+1);
    double wall_time = GetTimer()-start;
    double cpu_time  = (UserTime[1] - UserTime[0]) * 1000;
    printf("%s: %.0lf ms = %.0lf MiB/s,  cpu %.0lf ms = %.0lf%%,  os %.0lf ms\n",
        name, wall_time, (size / wall_time)*1000 / (1<<20),
        cpu_time, cpu_time/wall_time*100,
        (KernelTime[1]-KernelTime[0])*1000);
}


/***********************************************************************************************************************
*** GF(P) **************************************************************************************************************
************************************************************************************************************************/

// Twice wider type to hold intermediate MUL results
template<typename Type> struct Double           {};
template<>              struct Double<uint32_t> {typedef uint64_t T;};


template <typename T, T P>
T GF_Sub (T X, T Y)
{
    T res = X - Y;
    return res + (res>X)*P;   // res<=X? res : res+P
}

template <typename T, T P>
T GF_Add (T X, T Y)
{
    return GF_Sub<T,P> (X, P-Y);
}


#if __GNUC__ && defined(MY_CPU_64BIT)
// Alternative GF_Mul64 implementation for 64-bit GCC

#include <inttypes.h>
typedef unsigned __int128 uint128_t;
template<>              struct Double<uint64_t>    {typedef uint128_t T;};

// 4x wider type to hold intermediate MUL results
template<typename Type> struct Quadruple           {};
template<>              struct Quadruple<uint32_t> {typedef uint128_t T;};

// The binary logarithm, rounded down
template <typename T>
static constexpr int trunc_log2 (T x)
{return x<=1? 0 : 1+trunc_log2(x/2);}

template <typename T, T P>
T GF_Mul64 (T X, T Y)
{
    using DoubleT = typename Double<T>::T;
    using QuadT   = typename Quadruple<T>::T;

    // See chapter "16.9 Division" in the http://www.agner.org/optimize/optimizing_assembly.pdf
    constexpr int     BITS  = trunc_log2(P) + 8*sizeof(DoubleT);
    constexpr QuadT   invP2 = (QuadT(2) << BITS) / P;  // double the invP value
    constexpr DoubleT invP  = (invP2+1) / 2;           // rounded invP value
    constexpr DoubleT extra = 1 - (invP2 & 1);         // 1 if invP was rounded down, 0 otherwise

    DoubleT res = DoubleT(X)*Y;
    res -= DoubleT(((res+extra)*QuadT(invP)) >> BITS) * P;
    return T(res);
}

#elif _MSC_VER && defined(MY_CPU_64BIT)
// Alternative GF_Mul64 implementation made with MSVC intrinsics
#include <intrin.h>

template <typename T, T P>
T GF_Mul64 (T X, T Y)
{
    using DoubleT = typename Double<T>::T;
    DoubleT estInvP = ((DoubleT(1)<<63) / P) << 1;
    DoubleT invP    = (estInvP*P > (estInvP+1)*P? estInvP : estInvP+1);

    DoubleT res = DoubleT(X)*Y;
    res -= __umulh(res,invP) * P;
    return T(res>=P? res-P : res);
}

#else

// GF_Mul64 is optimized for 64-bit CPUs
template <typename T, T P>
T GF_Mul64 (T X, T Y)
{
    using DoubleT = typename Double<T>::T;
    return T((DoubleT(X)*Y) % P);
}

#endif

// GF_Mul32 is optimized for 32-bit CPUs, SIMD and GPUs
template <typename T, T P>
T GF_Mul32 (T X, T Y)
{
    using DoubleT = typename Double<T>::T;
    // invP32 := (2**64)/P - 2**32  :  if 2**31<P<2**32, then 2**32 < (2**64)/P < 2**33, and invP32 is a 32-bit value
    const DoubleT estInvP = ((DoubleT(1)<<63) / P) << 1;                        // == invP & (~1)
    const T       invP32  = T(estInvP*P > (estInvP+1)*P? estInvP : estInvP+1);  // we can't use 1<<64 for exact invP computation so we add the posible 1 in other way

    DoubleT res = DoubleT(X)*Y;
    res  -=  ((res + (res>>32)*invP32) >> 32) * P;    // The same as res -= ((res*invP) >> 64) * P, where invP = (2**64)/P, but optimized for 32-bit computations
    return T(res>=P? res-P : res);
}

#ifdef MY_CPU_64BIT
#define GF_Mul GF_Mul64
#else
#define GF_Mul GF_Mul32
#endif


template <typename T, T P>
T GF_Pow (T X, T N)
{
    T res = 1;
    for ( ; N; N/=2)
    {
        if (N&1)  res = GF_Mul<T,P> (res,X);
        X = GF_Mul<T,P> (X,X);
    }
    return res;
}


// Some primary root of 1 of power N
template <typename T, T P>
T GF_Root (T N)
{
    assert (P==0xFFF00001 || P==0x10001);
    T main_root  =  (P==0x10001? 3 : 19);  // root of power P-1 in the GF(P)

    assert ((P-1) % N  ==  0);
    return GF_Pow<T,P> (main_root, (P-1) / N);
}


template <typename T, T P>
T GF_Inv (T X)
{
    return GF_Pow<T,P> (X, P-2);
}



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
        std::cout << "\r" << i << "**" << std::hex << N << std::dec << "...";
        T q = GF_Pow<T,P> (i,N);
        if (q==1)
        {
            assert (P==0xFFF00001 || P==0x10001);

            if (P==0x10001) {
                if (1 == GF_Pow<T,P> (i,N/2))  goto next;
            } else {
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
            if (a != b)
            {
                std::cout << std::hex << "\r" << i << "*" << j << "=" << a << " != " << b << "\n" ;
                if (++n>10) return;
            }
        }
    }
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
    for (int n=0; n<1024; n++)
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
*** Number-Theoretical Transform in GF(P) ******************************************************************************
************************************************************************************************************************/

/* re-order data */
template <typename T, T P>
void revbin_permute (T** data, size_t n)
{
    if (n<=2)  return;
    size_t mr = 0; // the reversed 0
    for (size_t m=1; m<n; ++m) {
        // revbin_upd(r,n)
        size_t l = n;
        do {
        	l >>= 1;
        } while (mr+l >= n);
        mr = (mr & (l-1)) + l;

        if (mr > m) {
            std::swap (data[m], data[mr]);
        }
    }
}


// Recursive NTT implementation
template <typename T, T P>
void RecursiveNTT_Steps (T** data, size_t FirstN, size_t N, size_t SIZE, T* roots)
{
    N /= 2;
    if (N >= FirstN) {
#if _OPENMP>=200805
        #pragma omp task if (N>16384)
#endif
        RecursiveNTT_Steps<T,P> (data,   FirstN, N, SIZE, roots+1);
#if _OPENMP>=200805
        #pragma omp task if (N>16384)
#endif
        RecursiveNTT_Steps<T,P> (data+N, FirstN, N, SIZE, roots+1);
#if _OPENMP>=200805
        #pragma omp taskwait
#endif
    }

    T root = *roots,   root_i = 1;                      // zeroth root of power 2N of 1
    for (size_t i=0; i<N; i++) {
        T* __restrict__ block1 = data[i];
        T* __restrict__ block2 = data[i+N];
        for (size_t k=0; k<SIZE; k++) {                 // cycle over SIZE elements of the single block
            T temp    = GF_Mul<T,P> (block2[k], root_i);
            block2[k] = GF_Sub<T,P> (block1[k], temp);
            block1[k] = GF_Add<T,P> (block1[k], temp);
        }
        root_i = GF_Mul<T,P> (root_i, root);            // next root of power 2N of 1
    }
}


// Iterative NTT implementation
template <typename T, T P>
void IterativeNTT_Steps (T** data, size_t FirstN, size_t LastN, size_t SIZE, T* root_ptr)
{
    for (size_t N=FirstN; N<LastN; N*=2)
    {
        T root = *--root_ptr;
        for (size_t x=0; x<LastN; x+=2*N)
        {
            // first cycle optimized for root_i==1
            T* __restrict__ block1 = data[x];
            T* __restrict__ block2 = data[x+N];
            for (size_t k=0; k<SIZE; k++) {                     // cycle over SIZE elements of the single block
                T temp    = block2[k];                          // optimized for root_i==1
                block2[k] = GF_Sub<T,P> (block1[k], temp);
                block1[k] = GF_Add<T,P> (block1[k], temp);
            }

            // remaining cycles with root_i!=1
            T root_i = root;                                    // first root of power 2N of 1
            for (size_t i=1; i<N; i++) {
                T* __restrict__ block1 = data[x+i];
                T* __restrict__ block2 = data[x+i+N];
                for (size_t k=0; k<SIZE; k++) {                 // cycle over SIZE elements of the single block
                    T temp    = GF_Mul<T,P> (block2[k], root_i);
                    block2[k] = GF_Sub<T,P> (block1[k], temp);
                    block1[k] = GF_Add<T,P> (block1[k], temp);
                }
                root_i = GF_Mul<T,P> (root_i, root);            // next root of power 2N of 1
            }
        }
    }
}


// GF(P) NTT of N==2**X points of type T. Each point represented by SIZE elements (sequential in memory), so we perform SIZE transforms simultaneously
template <typename T, T P>
void Rec_NTT (size_t N, size_t SIZE, T** data, bool InvNTT)
{
    // Fill roots[] with roots of 1 of powers N, N/2, ... 2;  root_ptr points after the last entry
    T root = GF_Root<T,P>(N),  roots[66],  *root_ptr = roots;
    if (InvNTT)  root = GF_Inv<T,P>(root);
    while (root != 1) {
        *root_ptr++ = root;
        root = GF_Mul<T,P> (root, root);
    }

    revbin_permute<T,P> (data, N);

    #pragma omp parallel
    {
        // Smaller N values up to S are processed iteratively
#if defined(_OPENMP) && (_OPENMP < 200805)
        const size_t S = N/16;  // optimized for OpenMP 2.0 - do as much work as possible in the parallelized for loop
#else
        const size_t S = 1 << int(logb (99000/(SIZE*sizeof(T)) ));    // otherwise stay in L2 cache
#endif
        #pragma omp for
        for (int64_t i=0; i<N; i+=S)
            IterativeNTT_Steps<T,P> (data+i, 1, S, SIZE, root_ptr);

        // Larger N values are processed recursively
        #pragma omp master
        RecursiveNTT_Steps<T,P> (data, 2*S, N, SIZE, roots);
    }
}


// Iterative NTT implementation
template <typename T, T P>
void IterativeNTT (T** data, size_t N, size_t SIZE, T* root_ptr)
{
    revbin_permute<T,P> (data, N);
    IterativeNTT_Steps<T,P> (data, 1, N, SIZE, root_ptr);
}


// Transpose matrix R*C (rows*columns) into matrix C*R
template <typename T>
void TransposeMatrix (T* data, size_t R, size_t C)
{
    #pragma omp single
    if (R==C) {
        for (int r=0; r<R; r++) {
            for (size_t c=0; c<r; c++) {
                std::swap (data[r*C+c], data[c*R+r]);
            }
        }
    } else {
        T* tmp = new T [R*C];  std::unique_ptr<T> _tmp{tmp};
        for (int r=0; r<R; r++) {
            for (size_t c=0; c<C; c++) {
                tmp[c*R+r] = data[r*C+c];
            }
        }
        memcpy (data, tmp, R*C*sizeof(T));
    }
}


// The matrix Fourier algorithm (MFA)
template <typename T, T P>
void MFA_NTT (size_t N, size_t SIZE, T** data, bool InvNTT)
{
    // Split N-size problem into R rows * C columns
    size_t R = 1;   while (R*R < N)  R*=2;
    size_t C = N/R;

    // Fill roots[] with roots of 1 of powers N, N/2, ... 2;  root_ptr points after the last entry
    T root = GF_Root<T,P>(N),  roots[66],  *root_ptr = roots;
    if (InvNTT)  root = GF_Inv<T,P>(root);
    while (root != 1) {
        *root_ptr++ = root;
        root = GF_Mul<T,P> (root, root);
    }

    // Fill root_arr[i] with roots[0] ** i, where roots[0] ** N == 1
    T root_r = roots[0],  root_arr[1024];
    assert(1024 >= R);  // root_arr[] is too small
    for (size_t r=1; r<R; r++) {
        root_arr[r] = root_r;
        root_r = GF_Mul<T,P> (root_r, roots[0]);    // next root of power N
    }

    // MFA is impossible or inefficient
    if (N < 4  ||  N*SIZE < 256*1024)
    {
        IterativeNTT<T,P> (data, N, SIZE, root_ptr);
        return;
    }


    #pragma omp parallel
    {
        // 1. Apply a (length R) FFT on each column
        TransposeMatrix (data, R, C);
        #pragma omp for
        for (int64_t i=0; i<N; i+=R)
            IterativeNTT<T,P> (data+i, R, SIZE, root_ptr);
        TransposeMatrix (data, C, R);

        // 2. Multiply each matrix element (index r,c) by roots[0] ** (r*c)
        #pragma omp for
        for (int r=1; r<R; r++) {
            T root_c = root_arr[r];                             // roots[0] ** r
            for (size_t c=1; c<C; c++) {
                T* __restrict__ block = data[r*C+c];
                for (size_t k=0; k<SIZE; k++) {                 // cycle over SIZE elements of the single block
                    block[k] = GF_Mul<T,P> (block[k], root_c);
                }
                root_c = GF_Mul<T,P> (root_c, root_arr[r]);     // roots[0] ** r*c for the next c
            }
        }

        // 3. Apply a (length C) FFT on each row
        #pragma omp for
        for (int64_t i=0; i<N; i+=C)
            IterativeNTT<T,P> (data+i, C, SIZE, root_ptr);

        // 4. Transpose the matrix by transposing block pointers in the data[]
        TransposeMatrix (data, R, C);
    }
}


// Number theoretic transform by definition (slow!)
template <typename T, T P>
void Slow_NTT (size_t N, size_t SIZE, T* data, bool InvNTT)
{
    T *outdata = new T[N*SIZE];

    T root = GF_Root<T,P>(N);
    if (InvNTT)  root = GF_Inv<T,P>(root);

    T dw = 1;
    for (T i=0; i<N; ++i)
    {
        #pragma omp parallel for
        for (int k=0; k<SIZE; k++)      // cycle over SIZE elements of the single block
        {
            T t = 0;
            T w = 1;

            for (T x=0; x<N; ++x)
            {
                T tmp = GF_Mul<T,P> (w, data[x*SIZE+k]);
                t = GF_Add<T,P> (t, tmp);
                w = GF_Mul<T,P> (w, dw);
            }

            outdata[i*SIZE+k] = t;
        }

        dw = GF_Mul<T,P> (dw, root);      // next root of power N
    }

    memcpy (data, outdata, N*SIZE*sizeof(T));
    delete[] outdata;
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
void BenchNTT (bool RunMFA, bool RunCanonical, size_t N, size_t SIZE)
{
    T *data0 = new T[N*SIZE];
    for (size_t i=0; i<N*SIZE; i++)
        data0[i] = i%P;

    T **data = new T* [N];      // pointers to blocks
    for (size_t i=0; i<N; i++)
        data[i] = data0 + i*SIZE;

    uint32_t hash0 = hash(data, N, SIZE);    // hash of original data

    char title[99];
    for (int i=64; i--; )
        if (T(1)<<i == N)
            sprintf (title, "%s<2^%d,%.0lf>", RunCanonical?"Slow_NTT":RunMFA?"MFA_NTT":"Rec_NTT", i, SIZE*1.0*sizeof(T));

    if (RunCanonical) time_it (N*SIZE*sizeof(T), title, [&]{Slow_NTT<T,P> (N, SIZE, data0,false);});
    else if (RunMFA)  time_it (N*SIZE*sizeof(T), title, [&]{MFA_NTT <T,P> (N, SIZE, data, false);});
    else              time_it (N*SIZE*sizeof(T), title, [&]{Rec_NTT <T,P> (N, SIZE, data, false);});

    uint32_t hash1 = hash(data, N, SIZE);    // hash after NTT

    // Inverse NTT
    if (RunCanonical) Slow_NTT<T,P> (N, SIZE, data0,true);
    else if (RunMFA)  MFA_NTT <T,P> (N, SIZE, data, true);
    else              Rec_NTT <T,P> (N, SIZE, data, true);

    // Normalize the result by dividing by N
    T inv_N = GF_Inv<T,P>(N);
    for (size_t i=0; i<N*SIZE; i++)
        data0[i] = GF_Mul<T,P> (data0[i], inv_N);

    // Now we should have exactly the input data
    uint32_t hash2 = hash(data, N, SIZE);    // hash after NTT+iNTT
    if (hash2 == hash0)
        printf("Verified!  Original %.0lf,  after NTT: %.0lf  ", double(hash0), double(hash1));
    else
        printf("Checksum mismatch: original %.0lf,  after NTT: %.0lf,  after NTT+iNTT %.0lf", double(hash0), double(hash1), double(hash2));
}


// Parse cmdline and invoke appropriate benchmark/test routine
template <typename T, T P>
void Code (int argc, char **argv)
{
    char opt  =  (argc>=2?  argv[1][0] : ' ');
    if (opt=='i')  {Test_GF_Inv<T,P>();  return;}
    if (opt=='m')  {Test_GF_Mul<T,P>();  return;}
    if (opt=='r')  {FindRoot<T,P>(P-1);  printf ("GF_Root %s\n", GF_Root<T,P>(2)==P-1? "OK": "failed");  return;}
    if (opt=='d')  {DividersDensity<T,P>();  return;}
    if (opt=='b')  {time_it (10LL<<30, "Butterfly", [&]{BenchButterfly<T,P>();});  return;}

    size_t N = 1<<20;   // NTT order
    size_t SIZE = 512;  // Block size, in bytes
                        // 512 MB total

    if (argc>=3)  N = 1<<atoi(argv[2]);
    if (argc>=4)  SIZE = atoi(argv[3]);

    assert(N<P);  // Too long NTT for the such small P
    BenchNTT<T,P> (opt=='n', opt=='s', N, SIZE/sizeof(T));
}

// Deal with '=' argv[1] prefix: switch to P=0x10001
int main (int argc, char **argv)
{
    if (argc>=2 && argv[1][0]=='=')
    {
        argv[1]++;
        Code <uint32_t,0x10001> (argc, argv);
    } else {
        Code <uint32_t,0xFFF00001> (argc, argv);
    }
    return 0;
}

// MS: both P have the same speed, GCC: 65537 is 10% faster
// to do:
// "b N SIZE" in cmdline
// ntt32*.exe/GF_Mul32 doesn't work with 65537, probably due to hardcoded assumptions about P in GF_Mul32
// replace "(res>X)*P" in GF_Sub with bit arithmetics
// MS GF_Mul64 should became faster with the same algo as GCC one
