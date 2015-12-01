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

void generate_message(byte *data, int size)
{
    int *p = (int *)data;
    size >>=2;
    int i;
    for (i=0; i<size; i++) {
        *p++ = genrand();
    }
}

int compare_message(byte *p, byte *q, int size)
{
    int *pi = (int *)p;
    int *qi = (int *)q;
    size >>=2;
    int res=0;
    int i;
    for (i=0; i<size; i++) {
        if (pi[i]!=qi[i]) {
            res++;
            if (res<100) {
                printf("%d ",i*4);
            }
        }
    }
    printf("\n");
    return res;
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

// *******************************************************

#define GF ((1<<16)+1)

int *log, *exp, *inv;

int inline reduce(unsigned int a)
{
    a = (a& ((1<<16)-1)) - (a>>16);
    a += (((int)a)>>31)&GF;
    return a;
}

int inline field_mult(unsigned int a, unsigned int b)
{
    if (a==(1<<16)) return -(int)b + (((-(int)b)>>31) & GF);
    return reduce(a*b);
}

int inline field_mult_no(unsigned int a, unsigned int b)
{
    return reduce(a*b);
}

int inline field_diff(unsigned int a, unsigned int b)
{
    a -= b;
    return a + ((((int)a)>>31)&GF);
}

int inline field_sum(unsigned int a, unsigned int b)
{
    a -= GF-b;
    return a + ((((int)a)>>31)&GF);
}


void init_field()
{
    log = (int *) malloc(sizeof(int) * GF);
    exp = (int *) malloc(sizeof(int) * GF);
    inv = (int *) malloc(sizeof(int) * GF);

    int p = 1;
    int i;
    for (i=0; i+1<GF; i++) {
        exp[i]=p;
        log[p]=i;
        p = reduce(3*p);
    }
    exp[GF-1]=1;
    log[0]=0;

    for (i=0; i<GF; i++) {
        inv[i]=exp[GF-1-log[i]];
    }
    inv[0]=0;
    inv[1]=1;

    //test
    for (i=0; i<GF; i++) {
        if (field_mult(i, inv[i]) != 1) printf("%d ",i);
        if (exp[log[i]]!=i) printf("%d ", i);
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

void reverse(int *vect, int n)
{
    int i,j;

    j=n >> 1;
    for (i=1; i<n; i++)
    {
        if (j > i) {
            int temp=vect[i];
            vect[i]=vect[j];
            vect[j]=temp;
        }

        int m = n >> 1;
        while (m >= 1 && (j & m)) {
            j ^= m;
            m >>= 1;
        }
        j += m;
    }
}

void fft_dit(int *vect, int n)
{
    reverse(vect, n);

    int i,j;

    int step=1;
    int number=n/2;
    int mult=15;
    while (number>0)
    {
        int *p=vect;
        int *q=vect+step;
        for (i=0; i<number; i++) {
            for (j=0; j<step; j++) {
                int a = *p;
                int b = field_mult_no(*q, exp[j<<mult]);
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

void ifft_dit(int *vect, int n)
{
    reverse(vect, n);

    int i,j;

    int step=1;
    int number=n/2;
    int mult=15;
    while (number>0)
    {
        int *p=vect;
        int *q=vect+step;
        for (i=0; i<number; i++) {
            for (j=2*step; j>step; j--) {
                int a = *p;
                int b = field_mult_no(*q, exp[j<<mult]);
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


void fft2(int *vect, int n)
{
    int i=n/2;
    while (i--) {
        int a = *vect;
        int b = *(vect+1);
        *(vect++) = field_sum(a,b);
        *(vect++) = field_diff(a,b);
    }
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

void fft_rec(int *vect, int n)
{
//    if (n==1024) return fft_inc(vect, n);
    if (n==2) {
        int a = vect[0];
        int b = vect[1];
        vect[0] = field_sum(a,b);
        vect[1] = field_diff(a,b);
        return;
    }

    int i;
    int mult = 16 - get_log(n);

    n/=2;
    int *l = vect;
    int *h = vect + n;
    for (i=0; i<n; i++) {
        int a = *l;
        int b = *h;
        *(l++) = field_sum(a,b);
        *(h++) = field_mult_no(field_diff(a,b), exp[i<<mult]);
    }

    fft_rec(vect, n);
    fft_rec(vect+n, n);
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

void compute_prod_old(int *prod, int *pos, int k, int n)
{
    int x,i;
    for (x=0; x<n; x++) {
        long long t=1;
        for (i=0; i<k; i++) {
            if (x!=pos[i])
                t = field_mult(t, field_diff(x,pos[i]));
        }
        prod[x]=t;
    }
}

void compute_prod(int *prod, int *pos, int k, int n)
{
    int x,i;
    for (x=0; x<n; x++) {
        unsigned int t=0;
        for (i=0; i<k; i++) {
            t += log[field_diff(x,pos[i])];
        }
        prod[x]=exp[t % (1<<16)];
    }
}

//*********************************************************************
//*********************************************************************

int *temp, *prod;
int *inv_fft;
int *enc_fft;
int *prod_enc;

void init_code(int n)
{
    temp = (int *) malloc(sizeof(int) * n * 2);
    prod = (int *) malloc(sizeof(int) * n);
    inv_fft = (int *) malloc(sizeof(int) * n * 2);
    enc_fft = (int *) malloc(sizeof(int) * n);

    int x;
    for (x=0; x<n; x++) {
        inv_fft[x]=inv[x];
        enc_fft[x]=inv[x];
        if (x>0) inv_fft[2*n-x]=inv[GF-x];
    }
    fft(inv_fft,2*n);
    fft(enc_fft,n);
    for (x=0; x<n; x++) {
        inv_fft[x] = field_mult(inv_fft[x],inv[2*n]);
        inv_fft[n+x] = field_mult(inv_fft[n+x],inv[2*n]);
        enc_fft[x] = field_mult(enc_fft[x],inv[n]);
    }
    prod_enc = (int *) malloc(sizeof(int) * n);
}

//*********************************************************************
//*********************************************************************

//void fast_convolution(int *dst, int *src, int n)
//{
//    int x,y;
//
//    for (x=0; x<n; x++) {
//        ta[x]=src[x];
//        ta[n+x]=0;
//    }
//
//    fft(ta,2*n);
//    for (x=0; x<2*n; x++) {
//        ta[x] = field_mult(ta[x],inv_fft[x]);
//    }
//    ifft(ta,2*n);
//
//    for (x=0; x<n; x++) {
//        dst[x] = ta[x];
//    }
//}

void convolution(int *dst, int *src, int n)
{
    int x,y;
    for (x=0; x<n; x++) {
        int t=0;
        for (y=0; y<n; y++) {
            t = (t + field_mult(src[y], inv[(x + GF - y)%GF])) % GF;
        }
        dst[x]= t;
    }
}

//*********************************************************************
//*********************************************************************

// field_mult_no work here ?!
// means not both are -1 in all multiplication...
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
        temp[x] = field_mult(temp[x], enc_fft[x]);
    ifft(temp,n);

    // multiply by prod[x] the parity positions
    for (x=k; x<n; x++) {
        dst[x] = field_mult(temp[x], prod_enc[x]);
    }

    // put systematic symbol in place
    for (i=0; i<k; i++) dst[i]=src[i];
}

void decode(int *dst, int *src, int *pos, int k, int n)
{
    int i;
//    compute_prod(prod, pos, k, n);

    // put received packet in place
    // divide by prod[pos[i]]
    for (i=0; i<2*n; i++) temp[i]=0;
    for (i=0; i<k; i++) temp[pos[i]] = field_mult(src[i], inv[prod[pos[i]]]);

    // convolve with inverse function
    fft_rec_special(temp,2*n);
    for (i=0; i<2*n; i++)
        temp[i] = field_mult(temp[i], inv_fft[i]);
    ifft_rec_special(temp,2*n);

    // multiply by prod[i]
    for (i=0; i<k; i++)
        dst[i] = field_mult(temp[i], prod[i]);

    // replace received position by the correct symbol
    for (i=0; i<k; i++)
        if (pos[i]<k)
            dst[pos[i]]=src[i];
}


int main(int argc, char *argv[])
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
    init_field();
    init_code(N);

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
    compute_prod(prod_enc, kpos, K, N);

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

        encode(codeword, sys, K, N);

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

    compute_prod(prod, pos, K, N);
    for (n=0; n<nb_bloc; n++)
    {
        // current bloc
        int *rec = received + n*K;

        decode(codeword, rec, pos, K, N);

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
