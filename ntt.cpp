/// Implementation of Number-Theoretical transform in the GF(P) Galois field
#include <stdint.h>

typedef uint32_t ElementT;        // data items type, 32-bit unsigned integer for GF(P) computations with P>65536
typedef uint64_t DoubleElementT;  // twice wider type to hold intermediate results


// Reverse-binary reindexing
template <typename T>
void scramble (T* data, size_t nn)
{
    size_t n, mmax, m, j, istep, i;
    
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



template <ElementT P>
ElementT ADD (ElementT X, ElementT Y)
{
    ElementT res = X + Y;
    return (res>=P || res<X)? res-P : res;
}

template <ElementT P>
ElementT SUB (ElementT X, ElementT Y)
{
    ElementT res = X - Y;
    return res<=X? res : res+P;
}

template <ElementT P>
ElementT MUL (ElementT X, ElementT Y)
{
    return ElementT( (DoubleElementT(X)*Y) % P);
}



template <size_t N, typename T, T P>
class DanielsonLanczos {
   DanielsonLanczos<N/2,T,P> next;
public:
   void apply(T* data, T root) {
      T root_sqr = MUL<P> (root, root);  // first root of power N of 1
      next.apply (data,   root_sqr);
      next.apply (data+N, root_sqr);

      T root_i = root;   // first root of power 2N of 1 
      for (size_t i=0; i<N; i++) {
        T temp    = MUL<P> (root_i, data[i+N]);
        data[i+N] = SUB<P> (data[i], temp);
        data[i]   = ADD<P> (data[i], temp);
        root_i    = MUL<P> (root_i, root);  // next root of power 2N of 1
      }
   }
};
 
template<typename T, T P>
class DanielsonLanczos<0,T,P> {
public:
   void apply(T* data, T root) { }
};



template <size_t Exp, typename T, T P>
class GFFT {
   enum { N = 1<<Exp };
   DanielsonLanczos<N,T,P> recursion;
public:
   void fft(T* data) {
      scramble(data,N);
      recursion.apply(data,2781828);  /// to do: replace 2781828 with root of power 2N of 1
   }
};
