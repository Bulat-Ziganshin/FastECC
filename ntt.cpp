/// Implementation of the Number-Theoretical Transform in GF(P) Galois field
#include <iostream>
#include <algorithm>
#include <stdint.h>

#if defined(_M_X64) || defined(_M_AMD64) || defined(__x86_64__)
#define MY_CPU_AMD64
#endif

#if defined(MY_CPU_AMD64) || defined(_M_IA64)
#define MY_CPU_64BIT
#endif


/***********************************************************************************************************************
*** GF(P) **************************************************************************************************************
************************************************************************************************************************/

typedef uint32_t ElementT;        // data items type, 32-bit unsigned integer for GF(P) computations with P>65536
typedef uint64_t DoubleElementT;  // twice wider type to hold intermediate results


template <typename T, T P>
ElementT GF_Add (ElementT X, ElementT Y)
{
    ElementT res = X + Y;
    return res - ((res>=P)+(res<X))*P;   // (res>=P || res<X)? res-P : res;
}

template <typename T, T P>
ElementT GF_Sub (ElementT X, ElementT Y)
{
    ElementT res = X - Y;
    return res + (res>X)*P;   // res<=X? res : res+P
}

#if 0
// Alternative GF_Mul64 implementation for GCC - unfortunately, GCC 4.9 generates over-smart code for it

#include <inttypes.h>
typedef unsigned __int128 uint128_t;
typedef uint128_t FourElement;    // 4x wider type to hold intermediate MUL results

template <ElementT P>
ElementT GF_Mul64 (ElementT X, ElementT Y)
{
    DoubleElementT res = DoubleElementT(X)*Y;
    DoubleElementT invP = DoubleElementT((FourElement(1)<<64) / P);
    res -= DoubleElementT(((res)*FourElement(invP)) >> 64) * P;
    return ElementT(res>=P? res-P : res);
}

#elif _MSC_VER
// Alternative GF_Mul64 implementation made with MSVC intrinsics

template <ElementT P>
ElementT GF_Mul64 (ElementT X, ElementT Y)
{
    const DoubleElementT estInvP = ((DoubleElementT(1)<<63) / P) << 1;
    const DoubleElementT invP    = (estInvP*P > (estInvP+1)*P? estInvP : estInvP+1);
    DoubleElementT res = DoubleElementT(X)*Y;
    res -= __umulh(res,invP) * P;
    return ElementT(res>=P? res-P : res);
}

#else

// GF_Mul64 is optimized for 64-bit CPUs
template <ElementT P>
ElementT GF_Mul64 (ElementT X, ElementT Y)
{
    return ElementT( (DoubleElementT(X)*Y) % P);
}

#endif

// GF_Mul32 is optimized for 32-bit CPUs, SIMD and GPUs
template <ElementT P>
ElementT GF_Mul32 (ElementT X, ElementT Y)
{
    // invP32 := (2**64)/P - 2**32  :  if 2**31<P<2**32, then 2**32 < (2**64)/P < 2**33, and invP32 is a 32-bit value
    const DoubleElementT estInvP = ((DoubleElementT(1)<<63) / P) << 1;                          // == invP & (~1)
    const ElementT       invP32  = ElementT(estInvP*P > (estInvP+1)*P? estInvP : estInvP+1);    // we can't use 1<<64 for exact invP computation so we add the posible 1 in other way
    DoubleElementT res = DoubleElementT(X)*Y;
    res  -=  ((res + (res>>32)*invP32) >> 32) * P;    // The same as res -= ((res*invP) >> 64) * P, where invP = (2**64)/P, but optimized for 32-bit computations
    return ElementT(res>=P? res-P : res);
}

#ifdef MY_CPU_64BIT
#define GF_Mul GF_Mul64
#else
#define GF_Mul GF_Mul32
#endif


template <typename T, T P>
ElementT GF_Pow (T X, size_t N)
{
    T res = 1;
    for ( ; N; N/=2)
    {
        if (N&1)  res = GF_Mul<P> (res,X);
        X = GF_Mul<P> (X,X);
    }
    return res;
}



/***********************************************************************************************************************
*** Testing/benchmarking routines **************************************************************************************
************************************************************************************************************************/

// Find first root of 1 of power N
template <typename T, T P>
void FindRoot (size_t N)
{
    int cnt = 0;
    for (T i=2; i<P; i++)
    {
        std::cout << "\r" << i << "**" << std::hex << N << std::dec << "...";
        T q = GF_Pow<T,P> (i,N);
        if (q==1)
        {
            if (1 == GF_Pow<T,P> (i,N/2) ||
                1 == GF_Pow<T,P> (i,N/3) ||
                1 == GF_Pow<T,P> (i,N/5) ||
                1 == GF_Pow<T,P> (i,N/7) ||
                1 == GF_Pow<T,P> (i,N/13) ||
                0)
                goto next;
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


// Test the GF_Mul32 correctness
template <ElementT P>
void Test_GF_Mul32()
{
    int n = 0;
    for (ElementT i=0; i<P; i++)
    {
        if (i%4096==0)  std::cout << std::hex << "\r0x" << i << "...";
        for (ElementT j=0; j<i; j++)
            if (GF_Mul64<P> (i,j)  !=  GF_Mul32<P> (i,j))
            {
                std::cout << std::hex << "\r" << i << "*" << j << "=" << GF_Mul64<P> (i,j) << " != " << GF_Mul32<P> (i,j) << "\n" ;
                if (++n>10) return;
            }
    }
}


// Butterfly operation
template <typename T, T P>
void Butterfly (T* a, T* b, int TIMES, int SIZE, T root)
{
    for (int n=0; n<TIMES; n++) {
        for (int k=0; k<SIZE; k++) {                     // cycle over SIZE elements of the single block
            T temp = GF_Mul<P> (root, b[k]);
            b[k] = GF_Sub<T,P> (a[k], temp);
            a[k] = GF_Add<T,P> (a[k], temp);
        }
    }
}


// Benchmark speed of the Butterfly operation, processing 10 GB data == 500 MB processed with NTT<2**20>
template <typename T, T P>
int BenchButterfly()
{
    ElementT x=0;
    #pragma omp parallel for
    for (int n=0; n<8; n++)
    {
        const int sz = 1280;
        ElementT a[sz], b[sz];
        for (int i=0; i<sz; i++)
            a[i] = i*7+1, b[i] = i*15+8;
        Butterfly<ElementT,P> (a, b, 128*1024, sz, 1557);
        x += a[0];
    }
    return x?1:0;
}


/***********************************************************************************************************************
*** Number-Theoretical Transform in GF(P) ******************************************************************************
************************************************************************************************************************/

// Recursive NTT implementation
template <typename T, T P>
void RecursiveNTT (T* data, size_t FirstN, size_t N, size_t TOTAL, size_t SIZE, T* roots)
{
    N /= 2;  TOTAL /= 2;
    if (N >= FirstN) {
#if _OPENMP>=200805
        #pragma omp task if (N>16384)
#endif
        RecursiveNTT<T,P> (data,       FirstN, N, TOTAL, SIZE, roots+1);
#if _OPENMP>=200805
        #pragma omp task if (N>16384)
#endif
        RecursiveNTT<T,P> (data+TOTAL, FirstN, N, TOTAL, SIZE, roots+1);
#if _OPENMP>=200805
        #pragma omp taskwait
#endif
    }

    T root = *roots,   root_i = root;                   // first root of power 2N of 1
    for (size_t i=0; i<N*SIZE; i+=SIZE) {
        for (size_t k=0; k<SIZE; k++) {                 // cycle over SIZE elements of the single block
            size_t i1 = i+k, i2 = i+k+N*SIZE;
            T temp   = GF_Mul<P> (root_i, data[i2]);
            data[i2] = GF_Sub<T,P> (data[i1], temp);
            data[i1] = GF_Add<T,P> (data[i1], temp);
        }
        root_i = GF_Mul<P> (root_i, root);              // next root of power 2N of 1
    }
}


// Iterative NTT implementation
template <typename T, T P>
void IterativeNTT (T* data, size_t FirstN, size_t LastN, size_t TOTAL, size_t SIZE, T* root_ptr)
{
    if (FirstN == 1) {
        FirstN = 2;  --root_ptr;
        for (size_t x=0; x<TOTAL; x+=2*SIZE) {
            for (size_t k=0; k<SIZE; k++) {                     // cycle over SIZE elements of the single block
                size_t i1 = x+k, i2 = x+k+SIZE;
                T temp   = data[i2];   // root of power 2 is always -1, so we just exchange GF_Add and GF_Sub
                data[i2] = GF_Add<T,P> (data[i1], temp);
                data[i1] = GF_Sub<T,P> (data[i1], temp);
            }
        }
    }

    for (size_t N=FirstN; N<LastN; N*=2) {
        T root = *--root_ptr;
        for (size_t x=0; x<TOTAL; x+=2*N*SIZE) {
            T root_i = root;                                    // first root of power 2N of 1
            for (size_t i=0; i<N*SIZE; i+=SIZE) {
                for (size_t k=0; k<SIZE; k++) {                 // cycle over SIZE elements of the single block
                    size_t i1 = x+i+k, i2 = x+i+k+N*SIZE;
                    T temp   = GF_Mul<P> (root_i, data[i2]);
                    data[i2] = GF_Sub<T,P> (data[i1], temp);
                    data[i1] = GF_Add<T,P> (data[i1], temp);
                }
                root_i = GF_Mul<P> (root_i, root);              // next root of power 2N of 1
            }
        }
    }
}


// GF(P) NTT of N==2**X points of type T. Each point represented by SIZE elements (sequential in memory), so we perform SIZE transforms simultaneously
template <typename T, T P>
void NTT (size_t N, size_t SIZE, T* data)
{
    T root = 1557;                      // init 'root' with root of 1 of power 2**20 in GF(0xFFF00001)
    for (size_t i=1<<20; i>N; i/=2)
        root = GF_Mul<P> (root, root);  // find root of 1 of power N
    T roots[33], *root_ptr = roots;
    while (root!=1)
        *root_ptr++ = root,
        root = GF_Mul<P> (root, root);

    #pragma omp parallel
    {
        // Smaller N values up to S are processed iteratively
#if defined(_OPENMP) && (_OPENMP < 200805)
        const size_t S = N/16;  // optimized for OpenMP 2.0 - do as much work as possible in the paralleled for loop
#else
        const size_t S = 65536/(SIZE*sizeof(T));    // otherwise stay in L2 cache
#endif
        #pragma omp for
        for (int64_t i=0; i<N*SIZE; i+=S*SIZE)
            IterativeNTT<T,P> (data+i, 1, S, S*SIZE, SIZE, root_ptr);

        // Larger N values are processed recursively
        #pragma omp master
        RecursiveNTT<T,P> (data, 2*S, N, N*SIZE, SIZE, roots);
    }
}



int main()
{
    const ElementT P = 0xFFF00001;
    // Test_GF_Mul32<P>();
    // FindRoot<ElementT,P>(P-1);  // prints 19
    // if (BenchButterfly<ElementT,P>())  return 0;

    const size_t N = 1<<20;   // NTT order
    const size_t SIZE = 128;  // Block size, in 32-bit elements
    ElementT *data = new ElementT[N*SIZE];  // 512 MB
    for (int i=0; i<N*SIZE; i++)
        data[i] = i;

    NTT<ElementT,P> (N, SIZE, data);

    uint32_t sum = 314159253;
    for (int i=0; i<N*SIZE; i++)
        sum = (sum+data[i])*123456791;
    if (sum != 3267607014UL)
        printf("checksum failed: %.0lf", double(sum));

    return 0;
}
