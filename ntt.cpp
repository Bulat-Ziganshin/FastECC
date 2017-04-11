/// Implementation of the Number-Theoretical Transform in GF(P) Galois field
#include <iostream>
#include <algorithm>
#include <stdint.h>
#include <string.h>
#include <cassert>

#if defined(_M_X64) || defined(_M_AMD64) || defined(__x86_64__)
#define MY_CPU_AMD64
#endif

#if defined(MY_CPU_AMD64) || defined(_M_IA64)
#define MY_CPU_64BIT
#endif


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


// Find first root of 1 of power N
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
    for (int n=0; n<8; n++)
    {
        const int sz = 1280;
        T a[sz], b[sz];
        for (int i=0; i<sz; i++)
            a[i] = i*7+1, b[i] = i*15+8;
        Butterfly<T,P> (a, b, 128*1024, sz, 1557);
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
            T temp   = GF_Mul<T,P> (root_i, data[i2]);
            data[i2] = GF_Sub<T,P> (data[i1], temp);
            data[i1] = GF_Add<T,P> (data[i1], temp);
        }
        root_i = GF_Mul<T,P> (root_i, root);              // next root of power 2N of 1
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
                T temp   = data[i2];                            // here we use only root**0==1
                data[i2] = GF_Sub<T,P> (data[i1], temp);
                data[i1] = GF_Add<T,P> (data[i1], temp);
            }
        }
    }

    for (size_t N=FirstN; N<LastN; N*=2) {
        T root = *--root_ptr;
        for (size_t x=0; x<TOTAL; x+=2*N*SIZE) {
            T root_i = 1;                                       // zeroth power of the root of power 2N of 1
            for (size_t i=0; i<N*SIZE; i+=SIZE) {
                for (size_t k=0; k<SIZE; k++) {                 // cycle over SIZE elements of the single block
                    size_t i1 = x+i+k, i2 = x+i+k+N*SIZE;
                    T temp   = GF_Mul<T,P> (root_i, data[i2]);
                    data[i2] = GF_Sub<T,P> (data[i1], temp);
                    data[i1] = GF_Add<T,P> (data[i1], temp);
                }
                root_i = GF_Mul<T,P> (root_i, root);            // next root of power 2N of 1
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
        root = GF_Mul<T,P> (root, root);  // find root of 1 of power N
    T roots[33], *root_ptr = roots;
    while (root!=1)
        *root_ptr++ = root,
        root = GF_Mul<T,P> (root, root);

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

template <typename T, T P>
void MFA_NTT (size_t N, size_t SIZE, T* data)
{
    T root = 1557;                        // init 'root' with root of 1 of power 2**20 in GF(0xFFF00001)
    for (size_t i=1<<20; i>N; i/=2)
        root = GF_Mul<T,P> (root, root);  // find root of 1 of power N
    T roots[33], *root_ptr = roots;
    while (root!=1)
        *root_ptr++ = root,
        root = GF_Mul<T,P> (root, root);

    #pragma omp parallel
    {
        const size_t R = 1024, C = N/R;
        
        // 1. Apply a (length R) FFT on each column
        #pragma omp for
        for (int64_t i=0; i<N*SIZE; i+=R*SIZE)
            IterativeNTT<T,P> (data+i, 1, R, R*SIZE, SIZE, root_ptr);   // PROCESS COLUMN, NOT ROW!!!

        // 2. Multiply each matrix element (index r,c) by *roots ** (r*c)
        T root_r = *roots,  root_arr[R];        
        for (size_t r=1; r<R; r++) {
            root_arr[r] = root_r;
            root_r = GF_Mul<T,P> (root_r, *roots);              // next root of power N
        }

        #pragma omp for
        for (int r=1; r<R; r++) {
            T root_r = root_arr[r];
            for (size_t c=1; c<C; c++) {
                T root_c = root_r;
                T* sector = data + (r*C+c)*SIZE;
                for (size_t k=0; k<SIZE; k++) {                 // cycle over SIZE elements of the single block
                    sector[k] = GF_Mul<T,P> (sector[k], root_c);
                }
                root_c = GF_Mul<T,P> (root_c, root_r);          // next root of power N/R
            }
        }
        
        // 3. Apply a (length C) FFT on each row
        #pragma omp for
        for (int64_t i=0; i<N*SIZE; i+=C*SIZE)
            IterativeNTT<T,P> (data+i, 1, C, C*SIZE, SIZE, root_ptr);

        // 4. Transpose the matrix
        assert(R==C);  // traspose algo doesn't support R!=C
        #pragma omp for
        for (int r=0; r<R; r++) {
            for (size_t c=0; c<r; c++) {
                T* sector1 = data + (r*C+c)*SIZE;
                T* sector2 = data + (c*R+r)*SIZE;
                for (size_t k=0; k<SIZE; k++) {                 // cycle over SIZE elements of the single block
                    T tmp = sector1[k];
                    sector1[k] = sector2[k];
                    sector2[k] = tmp;
                }
            }
        }
    }
}


int main (int argc, char **argv)
{
    typedef uint32_t T;        // data items type, 32-bit unsigned integer for GF(P) computations with P>65536
    const T P = 0xFFF00001;
    char opt  =  (argc==2 && strlen(argv[1])==1?  argv[1][0] : ' ');
    if (opt=='i')  {Test_GF_Inv<T,P>(); return 0;}
    if (opt=='m')  {Test_GF_Mul<T,P>(); return 0;}
    if (opt=='r')  {FindRoot<T,P>(P-1); return 0;} // prints 19
    if (opt=='b')  {if (BenchButterfly<T,P>())  return 0;}

    const size_t N = 1<<20;   // NTT order
    const size_t SIZE = 128;  // Block size, in 32-bit elements
    T *data = new T[N*SIZE];  // 512 MB
    for (int i=0; i<N*SIZE; i++)
        data[i] = i;

//    NTT<T,P> (N, SIZE, data);
    MFA_NTT<T,P> (N, SIZE, data);

    uint32_t sum = 314159253;
    for (int i=0; i<N*SIZE; i++)
        sum = (sum+data[i])*123456791;
    if (sum != 562770418UL)
        printf("checksum failed: %.0lf", double(sum));

    return 0;
}
