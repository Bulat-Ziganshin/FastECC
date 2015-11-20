#include <stdint.h>

typedef uint32_t ElementT;  // data items type, 32-bit unsigned integer for GF(P) computations with P>65536


template <typename T=ElementT>
void scramble (T* data, size_t nn)
{
    size_t n, mmax, m, j, istep, i;
    
    // reverse-binary reindexing
    n = nn<<1;
    j=1;
    for (i=1; i<n; i+=2) {
        if (j>i) {
            swap(data[j-1], data[i-1]);
            swap(data[j], data[i]);
        }
        m = nn;
        while (m>=2 && j>m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    };
}


template <size_t N, typename T=ElementT>
class DanielsonLanczos {
   DanielsonLanczos<N/2,T> next;
public:
   void apply(T* data) {
      next.apply(data);
      next.apply(data+N);
 
      for (size_t i=0; i<N; i++) {
        T root = 1; /// to do: i'th root of power 2N of 1
        T temp = root * data[i+N];
        data[i+N] = data[i] - temp;
        data[i] += temp;
      }
   }
};
 
template<typename T>
class DanielsonLanczos<0,T> {
public:
   void apply(T* data) { }
};



template <size_t P, typename T=ElementT>
class GFFT {
   enum { N = 1<<P };
   DanielsonLanczos<N,T> recursion;
public:
   void fft(T* data) {
      scramble(data,N);
      recursion.apply(data);
   }
};
