/// Implementation of the Reed-Solomon algo in O(N*log(N)) using Number-Theoretical Transform in GF(p)
#include <iostream>
#include <algorithm>
#include <stdint.h>
#include <cstddef>
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
#include "ntt.cpp"


// Benchmark encoding using the Reed-Solomon algo
template <typename T, T P>
void EncodeReedSolomon (size_t N, size_t SIZE)
{
    T *data0 = VAlloc<T> (uint64_t(N)*SIZE);
    if (data0==0)  {printf("Can't alloc %.0lf MiB of memory!\n", (N/1048576.0)*SIZE*sizeof(T)); return;}

    for (size_t i=0; i<N*SIZE; i++)
        data0[i] = i%P;

    T **data = new T* [N];      // pointers to blocks
    for (size_t i=0; i<N; i++)
        data[i] = data0 + i*SIZE;

    char title[999];
    sprintf (title, "Reed-Solomon encoding (2^%.0lf source blocks => 2^%.0lf ECC blocks, %.0lf bytes each)", logb(N), logb(N), SIZE*1.0*sizeof(T));

    time_it (2.0*N*SIZE*sizeof(T), title, [&]
    {
        // 1. iNTT: polynomial interpolation. We find coefficients of order-N polynomial describing the source data
        MFA_NTT<T,P> (data, N, SIZE, true);
        // Now we should divide results by N in order to get coefficients, but we combined this operation with the multiplication below

        // Now we can evaluate the polynomial at 2*N points.
        // Points with even index will contain the source data,
        // while points with odd indexes may be used as ECC data.
        // But more efficient approach is to compute only odd-indexed points.
        // This is accomplished by the following steps:

        // 2. Multiply the polynomial coefficients by root(2*N)**i
        T root_2N = GF_Root<T,P>(2*N),  inv_N = GF_Inv<T,P>(N);
        #pragma omp parallel for
        for (ptrdiff_t i=0; i<N; i++) {
            T root_i = GF_Mul<T,P> (inv_N, GF_Pow<T,P>(root_2N,i));    // root_2N**i / N (combine division by N with multiplication by powers of the root)
            T* __restrict__ block = data[i];
            for (size_t k=0; k<SIZE; k++) {         // cycle over SIZE elements of the single block
                block[k] = GF_Mul<T,P> (block[k], root_i);
            }
        }

        // 3. NTT: polynomial evaluation. This evaluates the modified polynomial at root(N)**i points,
        // that is equivalent to evaluation of the original polynomial at root(2*N)**(2*i+1) points.
        MFA_NTT<T,P> (data, N, SIZE, false);

        // Further optimization: in order to compute only even-indexed points,
        // it's enough to compute order-N/2 NTT of data[i]+data[i+N/2]. And so on...
    });
}


int main (int argc, char **argv)
{
    size_t N = 1<<19;   // NTT order
    size_t SIZE = 2052; // Block size, in bytes
                        // 1 GB total

    if (argc>=2 && argv[1][0]=='.') {
        argv[1]++;
        verbose = false;
        if (argv[1][0]==0)  argv++, argc--;
    }
    if (argc>=2)  N = 1<<atoi(argv[1]);
    if (argc>=3)  SIZE = atoi(argv[2]);

    // InitLargePages();
    EncodeReedSolomon<uint32_t,0xFFF00001> (N,SIZE/sizeof(uint32_t));
}
