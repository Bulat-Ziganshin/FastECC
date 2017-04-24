/// Implementation of the Number-Theoretical Transform in GF(P) Galois field
#ifndef __GNUC__
#define __restrict__
#endif


/***********************************************************************************************************************
*** Small-order NTT codelets *******************************************************************************************
************************************************************************************************************************/

// Perform a single order-2 NTT
template <typename T, T P>
__forceinline void NTT2 (T& __restrict__ f0, T& __restrict__ f1)
{
    T u = f0, v = f1;
    f0 = GF_Add<T,P> (u, v);
    f1 = GF_Sub<T,P> (u, v);
}


// Perform a single order-3 NTT
template <typename T, T P, bool InvNTT>
__forceinline void NTT3 (T& __restrict__ f0, T& __restrict__ f1, T& __restrict__ f2)
{
    static constexpr T root    = GF_Root<T,P> (3);
    static constexpr T root1   = InvNTT? GF_Inv<T,P>(root) : root;
    static constexpr T root2   = GF_Mul <T,P> (root1, root1);
    static constexpr T const_1 = GF_Div <T,P> (GF_Add <T,P> (root1, root2), 2);
    static constexpr T const_2 = GF_Div <T,P> (GF_Sub <T,P> (root1, root2), 2);

    T u = GF_Add<T,P> (f1, f2);     // u = f1+f2
    T v = GF_Sub<T,P> (f1, f2);     // v = f1-f2
    T f1_f2 = u;                    // f1+f2

    u = GF_Mul<T,P> (u, const_1);   // u*(X+Y)/2,  X**3==1, Y==X**2
    v = GF_Mul<T,P> (v, const_2);   // v*(X-Y)/2
    u = GF_Add<T,P> (f0, u);        // f0 + u*(X+Y)/2

    f0 = GF_Add<T,P> (f0, f1_f2);   // f0 := f0+f1+f2
    f1 = GF_Add<T,P> (u, v);        // f1 := f0 + u*(X+Y)/2 + v*(X-Y)/2
    f2 = GF_Sub<T,P> (u, v);        // f2 := f0 + u*(X+Y)/2 - v*(X-Y)/2
}


// Perform a single order-4 NTT
template <typename T, T P, bool InvNTT>
__forceinline void NTT4 (T& __restrict__ f0, T& __restrict__ f1, T& __restrict__ f2, T& __restrict__ f3)
{
    static constexpr T root  = GF_Root<T,P> (4);
    static constexpr T root1 = InvNTT? GF_Inv<T,P>(root) : root;

    NTT2<T,P> (f0, f2);                // classic MFA algo
    NTT2<T,P> (f1, f3);
    f3 = GF_Mul<T,P> (f3, root1);
    NTT2<T,P> (f0, f1);
    NTT2<T,P> (f2, f3);
    std::swap(f1,f2);
}


#if NEVER  // require NTT(roots) which can't be efficiently implemented with current GF type model
// Perform a single order-5 NTT
template <typename T, T P, bool InvNTT>
__forceinline void NTT5 (T& __restrict__ f0, T& __restrict__ f1, T& __restrict__ f2, T& __restrict__ f3, T& __restrict__ f4)
{
    static constexpr T root  = GF_Root<T,P> (5);
    static constexpr T root1 = InvNTT? GF_Inv<T,P>(root) : root;
    static constexpr T root2 = GF_Mul <T,P> (root1, root1);
    static constexpr T root3 = GF_Mul <T,P> (root2, root1);
    static constexpr T root4 = GF_Mul <T,P> (root3, root1);

    NTT4<T,P, InvNTT> (f1, f3, f4, f2);
    // to do: NTT(roots)
    f1 = GF_Mul <T,P> (f1, root1);
    f2 = GF_Mul <T,P> (f1, root2);
    f3 = GF_Mul <T,P> (f1, root3);
    f4 = GF_Mul <T,P> (f1, root4);
    NTT4<T,P,!InvNTT> (f1, f3, f4, f2);
    std::swap(f2,f3);

    T u = f1;
    u = GF_Add <T,P> (u, f2);
    u = GF_Add <T,P> (u, f3);
    u = GF_Add <T,P> (u, f4);

    f1 = GF_Add <T,P> (f0, f1);
    f2 = GF_Add <T,P> (f0, f2);
    f3 = GF_Add <T,P> (f0, f3);
    f4 = GF_Add <T,P> (f0, f4);
    f0 = GF_Add <T,P> (f0, u);
}
#endif


// Perform a single order-6 NTT
template <typename T, T P, bool InvNTT>
__forceinline void NTT6 (T& __restrict__ f0, T& __restrict__ f1, T& __restrict__ f2,
                         T& __restrict__ f3, T& __restrict__ f4, T& __restrict__ f5)
{
    NTT3<T,P,InvNTT> (f0, f2, f4);
    NTT3<T,P,InvNTT> (f3, f5, f1);

    NTT2<T,P> (f0, f3);
    NTT2<T,P> (f2, f5);
    NTT2<T,P> (f4, f1);
}


// Perform a single order-9 NTT
template <typename T, T P, bool InvNTT>
__forceinline void NTT9 (T& __restrict__ f0, T& __restrict__ f1, T& __restrict__ f2,
                         T& __restrict__ f3, T& __restrict__ f4, T& __restrict__ f5,
                         T& __restrict__ f6, T& __restrict__ f7, T& __restrict__ f8)
{
    static constexpr T root  = GF_Root<T,P> (9);
    static constexpr T root1 = InvNTT? GF_Inv<T,P>(root) : root;
    static constexpr T root2 = GF_Mul <T,P> (root1, root1);
    static constexpr T root4 = GF_Mul <T,P> (root2, root2);

    // 4-step MFA on 3x3 matrix:

    // 1. NTT(3) on rows
    NTT3<T,P,InvNTT> (f0, f3, f6);
    NTT3<T,P,InvNTT> (f1, f4, f7);
    NTT3<T,P,InvNTT> (f2, f5, f8);

    // 2. Multiply by twiddle factors
    f4 = GF_Mul<T,P> (f4, root1);
    f5 = GF_Mul<T,P> (f5, root2);
    f7 = GF_Mul<T,P> (f7, root2);
    f8 = GF_Mul<T,P> (f8, root4);

    // 3. NTT(3) on columns
    NTT3<T,P,InvNTT> (f0, f1, f2);
    NTT3<T,P,InvNTT> (f3, f4, f5);
    NTT3<T,P,InvNTT> (f6, f7, f8);

    // 4. Transpose the matrix
    std::swap(f1,f3);
    std::swap(f2,f6);
    std::swap(f5,f7);
}


// Perform N order-2 NTTs
template <typename T, T P>
void NTT2 (T** data, size_t N, size_t SIZE)
{
    for (size_t i=0; i<N; i++)
        for (size_t k=0; k<SIZE; k++)         // cycle over SIZE elements of the single block
            NTT2<T,P> (data[i][k], data[i+N][k]);
}


// Perform N order-3 NTTs
template <typename T, T P, bool InvNTT>
void NTT3 (T** data, size_t N, size_t SIZE)
{
    for (size_t i=0; i<N; i++)
        for (size_t k=0; k<SIZE; k++)         // cycle over SIZE elements of the single block
            NTT3<T,P,InvNTT> (data[i][k], data[i+N][k], data[i+2*N][k]);
}


// Perform N order-4 NTTs
template <typename T, T P, bool InvNTT>
void NTT4 (T** data, size_t N, size_t SIZE)
{
    for (size_t i=0; i<N; i++)
        for (size_t k=0; k<SIZE; k++)         // cycle over SIZE elements of the single block
            NTT4<T,P,InvNTT> (data[i][k], data[i+N][k], data[i+2*N][k], data[i+3*N][k]);
}


// Perform N order-5 NTTs
template <typename T, T P, bool InvNTT>
void NTT5 (T** data, size_t N, size_t SIZE)
{
    for (size_t i=0; i<N; i++)
        for (size_t k=0; k<SIZE; k++)         // cycle over SIZE elements of the single block
            NTT5<T,P,InvNTT> (data[i][k], data[i+N][k], data[i+2*N][k], data[i+3*N][k], data[i+4*N][k]);
}


// Perform N order-6 NTTs
template <typename T, T P, bool InvNTT>
void NTT6 (T** data, size_t N, size_t SIZE)
{
    for (size_t i=0; i<N; i++)
        for (size_t k=0; k<SIZE; k++)         // cycle over SIZE elements of the single block
            NTT6<T,P,InvNTT> (data[i    ][k], data[i+  N][k], data[i+2*N][k],
                              data[i+3*N][k], data[i+4*N][k], data[i+5*N][k]);
}


// Perform N order-9 NTTs
template <typename T, T P, bool InvNTT>
void NTT9 (T** data, size_t N, size_t SIZE)
{
    for (size_t i=0; i<N; i++)
        for (size_t k=0; k<SIZE; k++)         // cycle over SIZE elements of the single block
            NTT9<T,P,InvNTT> (data[i    ][k], data[i+  N][k], data[i+2*N][k],
                              data[i+3*N][k], data[i+4*N][k], data[i+5*N][k],
                              data[i+6*N][k], data[i+7*N][k], data[i+8*N][k]);
}


/***********************************************************************************************************************
*** NTT steps **********************************************************************************************************
************************************************************************************************************************/

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


/***********************************************************************************************************************
*** Auxiliary NTT procedures *******************************************************************************************
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


/***********************************************************************************************************************
*** Three NTT implementations ******************************************************************************************
************************************************************************************************************************/

// GF(P) NTT of N==2**X points of type T. Each point represented by SIZE elements (sequential in memory), so we perform SIZE transforms simultaneously
template <typename T, T P>
void Rec_NTT (T** data, size_t N, size_t SIZE, bool InvNTT)
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
        const size_t S = 1 << int(logb (99000/(SIZE*sizeof(T)) ));    // otherwise stay in L2 cache (usually at least 256 KB / 2 threads minus memory lost due to only 4/8-associative hashing)
#endif
        #pragma omp for
        for (ptrdiff_t i=0; i<N; i+=S)
            IterativeNTT_Steps<T,P> (data+i, 1, S, SIZE, root_ptr);

        // Larger N values are processed recursively
        #pragma omp master
        RecursiveNTT_Steps<T,P> (data, 2*S, N, SIZE, roots);
    }
}


// The matrix Fourier algorithm (MFA)
template <typename T, T P>
void MFA_NTT (T** data, size_t N, size_t SIZE, bool InvNTT)
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
        for (ptrdiff_t i=0; i<N; i+=R)
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
        for (ptrdiff_t i=0; i<N; i+=C)
            IterativeNTT<T,P> (data+i, C, SIZE, root_ptr);

        // 4. Transpose the matrix by transposing block pointers in the data[]
        TransposeMatrix (data, R, C);
    }
}


// Number theoretic transform by definition (slow - O(N^2)!)
template <typename T, T P>
void Slow_NTT (T* data, size_t N, size_t SIZE, bool InvNTT)
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

        dw = GF_Mul<T,P> (dw, root);    // next root of power N
    }

    memcpy (data, outdata, N*SIZE*sizeof(T));
    delete[] outdata;
}


#if NEVER
// Generic driver handling any orders
template <typename T, T P>
void Generic_NTT (T* data, size_t N, size_t SIZE, bool InvNTT)
{
    R = max. 2^n power of N
    C = N/R
    Employ PFA with R*C elements
    0. fill roots[] with unity roots starting with root**R==1
    1. MFA_NTT handles R-sized NTTs
    2. PFA_NTT handles C-sized NTTs
}
#endif
