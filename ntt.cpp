#define _USE_MATH_DEFINES
#include <cmath> 


template<typename T=double>
void scramble(T* data, unsigned nn)
{
    unsigned n, mmax, m, j, istep, i;
    
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


template<unsigned N, typename T=double>
class DanielsonLanczos {
   DanielsonLanczos<N/2,T> next;
public:
   void apply(T* data) {
      next.apply(data);
      next.apply(data+N);
 
      T tempr,tempi,c,s;
 
      for (unsigned i=0; i<N; i+=2) {
        c = cos(i*M_PI/N);
        s = -sin(i*M_PI/N);
        tempr = data[i+N]*c - data[i+N+1]*s;
        tempi = data[i+N]*s + data[i+N+1]*c;
        data[i+N] = data[i]-tempr;
        data[i+N+1] = data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
   }
};
 
template<typename T>
class DanielsonLanczos<1,T> {
public:
   void apply(T* data) { }
};





template<unsigned P,
         typename T=double>
class GFFT {
   enum { N = 1<<P };
   DanielsonLanczos<N,T> recursion;
public:
   void fft(T* data) {
      scramble(data,N);
      recursion.apply(data);
   }
};
