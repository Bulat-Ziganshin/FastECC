// ************************************************
// ************************************************

// Sample program to use reed_solomon.c
// (c) 2009 Frederic Didier.

#include "stdio.h"
#include "stdlib.h"
#include "time.h"

typedef unsigned char byte;

/*****************************************************************/
/** Random number generator -> 32bits                           **/
/** https://en.wikipedia.org/wiki/Linear_congruential_generator **/
/*****************************************************************/

static unsigned long last_rnd = 0;

void sgenrand (unsigned long seed)
{
    last_rnd = seed;
}

unsigned int genrand()
{
    last_rnd = 29943829*last_rnd + 1013904223;
    return last_rnd;
}

double double_genrand() {
    return genrand() * (1.0/4294967295.0);
}

// *******************************************************

// The Art of Computer programming - Knuth
// vol 2 - section 3.4.2 page 137
// Algorithm S (Selection sampling technique)

void generate_positions(int N, int K, int *pos)
{
    int size=N;
    int w=K;
    do {
        if (double_genrand()*size <= w) {
            pos[K-w] = N-size;
            w--;
        }
        size--;
    } while (size);
}

// *******************************************************

double get_sec(clock_t diff)
{
    return (double)diff / (double)CLOCKS_PER_SEC;
}

double get_rate(clock_t diff, int byte)
{
    return (double)(byte)/((double)(1<<20) * get_sec(diff));
}

double get_KB(int byte)
{
    return (double)byte/(double)(1<<10);
}

/***********************************************************************************************************************
*** GF(P) **************************************************************************************************************
************************************************************************************************************************/

#include <stdint.h>
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

template <typename T, T P>
ElementT GF_Mul (T X, T Y)
{
    return T( (DoubleElementT(X)*Y) % P);
}

template <typename T, T P>
ElementT GF_Pow (T X, T N)
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
ElementT GF_Inv (T X)
{
    return GF_Pow<T,P> (X,P-2);
}

// *******************************************************

#define GF ((1<<16)+1)

int inline field_mult(unsigned int a, unsigned int b)
{
    return GF_Mul<ElementT,GF>(a,b);
}

int inline field_mult_no(unsigned int a, unsigned int b)
{
    return GF_Mul<ElementT,GF>(a,b);
}

int inline field_diff(unsigned int a, unsigned int b)
{
    return GF_Sub<ElementT,GF>(a,b);
}

int inline field_sum(unsigned int a, unsigned int b)
{
    return GF_Add<ElementT,GF>(a,b);
}


int *exp, *inv;

template <typename T, T P>
void init_field()
{
    exp = (int *) malloc(sizeof(int) * P);
    inv = (int *) malloc(sizeof(int) * P);

    T p = 1;
    for (T i=0; i+1<P; i++) {
        exp[i]=p;
        p = GF_Mul<T,P>(3,p);
    }
    exp[P-1]=1;

    for (T i=0; i<P; i++) {
        inv[i] = GF_Inv<T,P>(i);
    }
    inv[0]=0;
    inv[1]=1;

    //test
    for (T i=0; i<P; i++) {
        if (field_mult(i, inv[i]) != 1) printf("%d ",i);
    }
    printf("\n");
}

//*********************************************************************
//*********************************************************************

// return log_2(b);
int get_log(int n)
{
    int i=0;
    while (n>>(i+1))
        i++;
    return i;
}

// decimation in frequency
void fft_inc(int *vect, int n)
{
    int i,j;

    int number=1;
    int step=n/2;
    int mult = 16 - get_log(n);

    while (step>0)
    {
        int *p=vect;
        int *q=vect+step;
        for (i=0; i<number; i++) {
            for (j=0; j<step; j++) {
                int a = *p;
                int b = *q;
                *(p++) = field_sum(a,b);
                *(q++) = field_mult_no(field_diff(a,b), exp[j<<mult]);
            }
            p+=step;
            q+=step;
        }
        step>>=1;
        number<<=1;
        mult++;
    }

//    reverse(vect,n);
}

// decimation in time
void ifft_inc(int *vect, int n)
{
//    reverse(vect,n);

    int i,j;
    int number=n/2;

    int step=1;
    int mult=15;

    int *root=exp + (1<<16);
    while (number>0)
    {
        int *p=vect;
        int *q=vect+step;
        for (i=0; i<number; i++) {
            for (j=0; j<step; j++) {
                int a = *p;
                int b = field_mult_no(*q, *(root - (j<<mult)));
                *(p++) = field_sum(a,b);
                *(q++) = field_diff(a,b);
            }
            p+=step;
            q+=step;
        }
        step<<=1;
        number>>=1;
        mult--;
    }
}

void fft_rec_special(int *vect, int n)
{
    int i;
    int mult = 16 - get_log(n);

    n/=2;
    int *l = vect;
    int *h = vect + n;
    for (i=0; i<n; i++) {
        int a = *(l++);
        *(h++) = field_mult_no(a, exp[i<<mult]);
    }

    fft_inc(vect, n);
    fft_inc(vect+n, n);
}

void ifft_rec_special(int *vect, int n)
{
    int i;
    int mult = 16 - get_log(n);

    n/=2;
    ifft_inc(vect, n);
    ifft_inc(vect+n, n);

    int *l = vect;
    int *h = vect + n;
    int *root = exp + (1<<16);
    for (i=0; i<n; i++) {
        int a = *l;
        int b = field_mult_no(*h, *(root - (i<<mult)));
        *(l++) = field_sum(a,b);
        h++;
    }
}

void (*fft)(int *, int) = fft_inc;
void (*ifft)(int *, int) = ifft_inc;

//*********************************************************************
//*********************************************************************

template <typename T, T P>
void compute_prod(int *prod, int *pos, int k, int n)
{
    int x,i;
    for (x=0; x<n; x++) {
        T t=1;
        for (i=0; i<k; i++) {
            if (x!=pos[i])
                t = GF_Mul<T,P> (t, GF_Sub<T,P>(x,pos[i]));
        }
        prod[x]=t;
    }
}

//*********************************************************************
//*********************************************************************

int *temp, *prod;
int *inv_fft;
int *enc_fft;
int *prod_enc;

template <typename T, T P>
void init_code(int n)
{
    temp = (int *) malloc(sizeof(int) * n * 2);
    prod = (int *) malloc(sizeof(int) * n);
    inv_fft = (int *) malloc(sizeof(int) * n * 2);
    enc_fft = (int *) malloc(sizeof(int) * n);

    int x;
    for (x=0; x<n; x++) {
        inv_fft[x]=GF_Inv<T,P>(x);
        enc_fft[x]=GF_Inv<T,P>(x);
        if (x>0)
            inv_fft[2*n-x] = GF_Inv<T,P> (GF_Sub<T,P> (0,x));
    }
    fft(inv_fft,2*n);
    fft(enc_fft,n);
    for (x=0; x<n; x++) {
        inv_fft[x]   = GF_Mul<T,P> (inv_fft[x],   GF_Inv<T,P>(2*n));
        inv_fft[n+x] = GF_Mul<T,P> (inv_fft[n+x], GF_Inv<T,P>(2*n));
        enc_fft[x]   = GF_Mul<T,P> (enc_fft[x],   GF_Inv<T,P>(n));
    }
    prod_enc = (int *) malloc(sizeof(int) * n);
}

//*********************************************************************
//*********************************************************************

// field_mult_no work here ?!
// means not both are -1 in all multiplication...
template <typename T, T P>
void encode(int *dst, int *src, int k, int n)
{
    int x,i;

    // put received packet in place
    // divide by prod[pos[i]]
    for (x=k; x<n; x++) temp[x] = 0;
    for (i=0; i<k; i++) temp[i] = field_mult_no(src[i], inv[prod_enc[i]]);

    // convolve with inverse function
    fft(temp,n);
    for (x=0; x<n; x++)
        temp[x] = GF_Mul<T,P>(temp[x], enc_fft[x]);
    ifft(temp,n);

    // multiply by prod[x] the parity positions
    for (x=k; x<n; x++) {
        dst[x] = GF_Mul<T,P>(temp[x], prod_enc[x]);
    }

    // put systematic symbol in place
    for (i=0; i<k; i++) dst[i]=src[i];
}

template <typename T, T P>
void decode(int *dst, int *src, int *pos, int k, int n)
{
    int i;
//    compute_prod(prod, pos, k, n);

    // put received packet in place
    // divide by prod[pos[i]]
    for (i=0; i<2*n; i++) temp[i]=0;
    for (i=0; i<k; i++) temp[pos[i]] = GF_Mul<T,P>(src[i], inv[prod[pos[i]]]);

    // convolve with inverse function
    fft_rec_special(temp,2*n);
    for (i=0; i<2*n; i++)
        temp[i] = GF_Mul<T,P>(temp[i], inv_fft[i]);
    ifft_rec_special(temp,2*n);

    // multiply by prod[i]
    for (i=0; i<k; i++)
        dst[i] = GF_Mul<T,P>(temp[i], prod[i]);

    // replace received position by the correct symbol
    for (i=0; i<k; i++)
        if (pos[i]<k)
            dst[pos[i]]=src[i];
}


template <typename T, T P>
int test(int argc, char *argv[])
{
    int i,j,n;
    clock_t tick;

    // get parameters
    int N,K,S, nb_bloc;

    // help message
    if (argc<=3) {
        printf("usage: %s K N S\n",argv[0]);
        return 0;
    }

    K = atoi(argv[1]);
    N = atoi(argv[2]);
    nb_bloc = atoi(argv[3]);
    S = 2;

    // power of two just greater than N
    int n_walsh=1;
    while ((1<<n_walsh) < N) n_walsh++;
    N = 1<<n_walsh;

    printf("GF(%d)\n", GF);
    printf("K=%d\n",K);
    printf("N=%d\n",N);
    printf("nb_bloc=%d\n", nb_bloc);
    printf("message size = %f KB\n", get_KB(S*K*nb_bloc));

    // ****************************************
    // ****************************************
    printf("[initialisation (memory + randomness)]\n");
    tick  = clock();

    // code init
    init_field<T,P>();
    init_code <T,P>(N);

    // memory for the full message
    int *positions = (int *)malloc(sizeof(int)*K*nb_bloc);
    int *message   = (int *)malloc(sizeof(int)*K*nb_bloc);
    int *received  = (int *)malloc(sizeof(int)*K*nb_bloc);


    // init random number generator
//     sgenrand(time(NULL));
    sgenrand(123);
//    sgenrand(321);

    // Generate the random message
    for (i=0; i<K*nb_bloc; i++) {
        message[i] = genrand() % GF;
    }

    // Generate the random positions
//    for (n=0; n<nb_bloc; n++) {
//        generate_positions(N, K, positions + n*K);
//    }
    generate_positions(N, K, positions);

    // memory for encoding/decoding one bloc
    int *codeword = (int *)malloc(sizeof(int)*N);

    // for encoding
    int *kpos=(int *)malloc(sizeof(int)*K);
    for(i=0; i<K; i++) kpos[i] = i;
    compute_prod<T,P> (prod_enc, kpos, K, N);

    // end of initialisation
    tick = clock() - tick;
    printf("%f s\n", get_sec(tick));
    printf("%f MB/s\n", get_rate(tick, S*K*nb_bloc));
    printf("\n");

    // ****************************************
    // ****************************************
    printf("[encoding message]\n");
    tick = clock();

    int *pos = positions;
    for (n=0; n<nb_bloc; n++)
    {
//        int *pos = positions + n*K;
        int *sys = message + n*K;
        int *rec = received + n*K;

        encode<T,P> (codeword, sys, K, N);

        // simulate errors...
        for (i=0; i<K; i++) rec[i] = codeword[pos[i]];
    }

    tick = clock() - tick;
    printf("%f s\n", get_sec(tick));
    printf("%f MB/s\n", get_rate(tick, S*K*nb_bloc));
    printf("%f MB/s\n", get_rate(tick, S*N*nb_bloc));
    printf("\n");

    // ****************************************
    // ****************************************
    printf("[decoding]\n");
    tick = clock();

    double syst=0;

    compute_prod<T,P> (prod, pos, K, N);
    for (n=0; n<nb_bloc; n++)
    {
        // current bloc
        int *rec = received + n*K;

        decode<T,P> (codeword, rec, pos, K, N);

        // put result back into received
        for (i=0; i<K; i++) {
            if (pos[i]<K) syst++;
            rec[i]=codeword[i];
        }
    }

    tick = clock() - tick;
    printf("%f s\n", get_sec(tick));
    printf("%f MB/s\n", get_rate(tick, S*K*nb_bloc));
    printf("%f percent of systematic packets received\n", syst / (double)(K*nb_bloc));
    printf("\n");

    // ****************************************
    // ****************************************

    // verify that we recovered the message correctly
    printf("[errors]\n");
    int count=0;
    for (i=0; i<K*nb_bloc; i++) {
        if (message[i]!=received[i]) count++;
    }
    printf("%i\n",count);

    // ****************************************
    // ****************************************

    // end;
    return 0;
}

int main(int argc, char *argv[])
{
    return test<ElementT,GF>(argc,argv);
}
