
### Program usage

`RS [N=19 [SIZE=2052]]` - benchmark NTT-based Reed-Solomon encoding using 2^N input (data) blocks and 2^N output (parity) blocks, each block SIZE bytes long


### Prior art

The encoding and decoding algorithms implemented by FastECC were described in the paper
[An Efficient (n,k) Information Dispersal Algorithm based on Fermat Number Transforms](https://pdfs.semanticscholar.org/141d/c4ee4cca45b4ed1c07f890f758e427597db8.pdf)
published in 2013 by Sian-Jheng Lin and Wei-Ho Chung.

The following is my own investigations written prior to reading this great paper :)


### Encoding in O(N*log(N))

In this scheme, we consider N input words as values at N fixed points of some polynomial p(x) with order<N,
and for parity words compute values of p(x) at another M fixed points ("fixed" here means that these points depend only on N and M).
It requires intermediate step of computing N polynomial coefficients.
We can compute polynomial coefficients from values with iNTT(N) (i.e. order-N inverse NTT), and then compute polynomial values with NTT(M1) where M1>=max(N,M).

But again, since we can compute (multiplicative) NTT in GF(p) only with orders dividing p-1, the actual algorithm is:
1. Find N1>=N such that N divides p-1
2. Extend input vector of size N with zeroes (considered as polynomial values at extra points) to the size N1
3. Perform order-N1 iNTT
4. Find M1>=max(N1,M) such that M1 divides p-1
5. Extend intermediate vector of size N1 with zeroes (considered as higher polynomial coefficients) to the size M1
6. Perform order-M1 NTT
7. Output some M values from NTT result as parity data

In order to allow fast decoding, the M1 points employed in the second NTT should inlcude the N1 points from the first iNTT, that requires M1=k*N1.
So, we need to correct fourth step of the algorithm:

4. Find M1>=N1+M such that M1=k*N1 and M1 divides p-1. This means that M1 points will include N1 points from the first iNTT and at least M other points.



### [Speeding up RS-encoding (external link)](https://www.livebusinesschat.com/smf/index.php?topic=5952.msg44167#msg44167)


### Decoding in O(N*log(N))

On decoding stage, we have N remaining values out of N+M encoded, and we know that they represent polynomial f(x) with only N coefficients:
```
f(x) = a[0]+...+a[N-1]*x^(N-1))
```
So, we have to solve classic "polynomial interpolation" problem - restore f(x) from N equations:
```
f(x[i]) = y[i], i=1..N
```
Googling by "fast polynomial interpolation", we can find that it can be solved in generic way in O(N*log(N)^2) time.
But we can do it faster by mentioning that those N remaining points is a subset of N+M points that can be used for NTT/iNTT evaluation.
This allows us to build even more efficient algorithm:

1. Compute special polynomial l(x) of order==M
2. Multiply them: p(x) = f(x)*l(x)
3. Compute f(x) = p(x) / l(x)

On the first sight, it looks meaningless. The catch is that l(x) is so special that f(x)*l(x) product can be fully defined although we know f(x) values only at N points.

Let's denote points where f(x) values was lost as e[i], i=1..M. So, we have:
```
f(x[1..N],e[1..M]) = (y[1..N], ?[1..M])
```
We will use the following polynomial `l(x) = (x-e[1])*....*(x-e[M])`. It will have the following values at the same points:
```
l(x[1..N],e[1..M]) = (l[1..N], 0,0..0)
```
And the p(x) = f(x)*l(x) product has the following values:
```
p(x[1..N],e[1..M]) = (y[1]*l[1], ... y[N]*l[N], 0,0..0)
```

As you see, we know p(x) values at N+M points. We also know that it has order<N+M, because l(x) has order M and f(x) has order<N.
So, now we can use NTT of order N+M to convert values at those N+M points into coefficients of p(x) polynomial.

Once we have computed p(x), we can employ usual "polynomial division" algorithm to compute f(x) = p(x) / l(x) in the O(N*log(N)) time.
As result, we will get all coefficients of f(x) and then call NTT again to compute f(x) values corresponding to the source data
