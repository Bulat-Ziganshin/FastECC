/// Implementation of the Number-Theoretical Transform in GF(P) Galois field
#include <iostream>
#include <algorithm> 
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
            std::swap(data[j-1], data[i-1]);
            std::swap(data[j], data[i]);
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
ElementT GF_Add (ElementT X, ElementT Y)
{
    ElementT res = X + Y;
    return (res>=P || res<X)? res-P : res;
}

template <ElementT P>
ElementT GF_Sub (ElementT X, ElementT Y)
{
    ElementT res = X - Y;
    return res<=X? res : res+P;
}

template <ElementT P>
ElementT GF_Mul (ElementT X, ElementT Y)
{
    return ElementT( (DoubleElementT(X)*Y) % P);
}



template <size_t N, typename T, T P>
class DanielsonLanczos {
   DanielsonLanczos<N/2,T,P> next;
public:
   void apply(T* data, T root) {
      T root_sqr = GF_Mul<P> (root, root);  // first root of power N of 1
      if (N>8192) {
          #pragma omp task
          next.apply (data,   root_sqr);
          #pragma omp task
          next.apply (data+N, root_sqr);
          #pragma omp taskwait
      } else {
          next.apply (data,   root_sqr);
          next.apply (data+N, root_sqr);
      }

      T root_i = root;   // first root of power 2N of 1 
      for (size_t i=0; i<N; i++) {
        T temp    = GF_Mul<P> (root_i, data[i+N]);
        data[i+N] = GF_Sub<P> (data[i], temp);
        data[i]   = GF_Add<P> (data[i], temp);
        root_i    = GF_Mul<P> (root_i, root);  // next root of power 2N of 1
      }
   }
};
 
template<typename T, T P>
class DanielsonLanczos<0,T,P> {
public:
   void apply(T* data, T root) { }
};



template <size_t Exp, typename T, T P>
class NTT {
   enum { N = 1<<Exp };
   DanielsonLanczos<N,T,P> recursion;
public:
   void ntt(T* data) {
      T root = 1557;                      // init 'root' with root of 1 of power 2**20 in GF(0xFFF00001)
      for (int i=20; --i>Exp; )
          root = GF_Mul<P> (root, root);  // find root of 1 of power 2N
      scramble(data,N);
      recursion.apply(data,root);
   }
};



// Find first root of 1 of power 2**N
template <typename T, T P>
void FindRoot (int N)
{
    for (T i=2; i<P; i++)
    {
        T q = i;
        for (int j=0; j<N; j++)
        {
            if (q==1)  goto next;
            q = GF_Mul<P> (q,q);
        }
        if (q==1)
        { 
            std::cout << i << "\n";
            break;
        }
next:;        
    }
}



int main()
{
    const ElementT P = 0xFFF00001;
    // FindRoot<ElementT,P>(20);  // prints 1557

    ElementT *data = new ElementT[1<<20];
    for (int i=0; i<(1<<20); i++)
        data[i] = i;
    NTT<19,ElementT,P> transformer;
    #pragma omp parallel num_threads(16)
    #pragma omp single
    for (int i=0; i<100; i++)
        transformer.ntt(data);
    return 0;
}
