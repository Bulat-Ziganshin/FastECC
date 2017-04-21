
### Program usage

`RS [N=19 [SIZE=1024]]` - benchmark NTT-based Reed-Solomon encoding using 2^N input (source) blocks and 2^N output (ECC) blocks, each block SIZE bytes long


### Lacan scheme - encoding in O(N*logN)

Actually, it seems that this scheme was not proposed by Plank, but anyway it was known at least since 2004.

In the new scheme, we consider N input words as values at N fixed points of some polynom p(x) with order<N, and compute for ECC data values of p(x) at another M fixed points ("fixed" here means that these points depend only on N and M). It requires intermediate step of computing N polynom coefficients. We can compute polynom coefficients from values with IFFT(N) (i.e. order-N Inverse FFT), and then compute polynom values with FFT(M1) where M1>=max(N,M).

But again, since we can compute (multiplicative) FFT in GF(p) only with orders dividing p-1, the actual algorithm is:
1. Find N1>=N such that N divides p-1
2. Extend input vector of size N with zeroes (considered as polynom values at extra points) to the size N1
3. Perform order-N1 IFFT
4. Find M1>=max(N1,M) such that M1 divides p-1
5. Extend intermediate vector of size N1 with zeroes (considered as higher polynom coefficients) to the size M1
6. Perform order-M1 FFT
7. Output some M values from FFT result as ECC data

In order to allow fast decoding, the M1 points employed in the second FFT should inlcude the N1 points from the first IFFT, that requires M1=k*N1. So, we need to correct fourth step of the algorithm:

4. Find M1>=N1+M such that M1=k*N1 and M1 divides p-1. This means that M1 points will include N1 points from the first IFFT and at least M other points.



### [Speeding up RS-encoding (external link)](https://www.livebusinesschat.com/smf/index.php?topic=5952.msg44167#msg44167)


### Lacan scheme - decoding in O(N*logN)

On decoding stage, we have N remaining values out of N+M encoded, and we know that they represent polynom p(x) with only N coefficients:
```
p(x) = a[0]+...+a[N-1]*x^(N-1))
```
So, we have to solve classical "polynom interpolation" problem - restore p(x) from N equations:
```
p(x[j]) = y[j], j=1..N
```
Googling by "fast polynom interpolation", we can find that it can be solved in generic way in O(N*log(N)^2) time. But we can do it faster by mentioning that those N remaining points is a subset of N+M points that can be used for FFT/iFFT evaluation. This allows us to build even more efficient algorithm:

1. Compute special polynom g(x) of order==M
2. Multiply them: h(x) = p(x)*g(x)
3. Compute p(x) = h(x) / g(x)

On the first sight, it looks meaningless. The catch is that g(x) is so special that p(x)*g(x) product can be fully although we know p(x) values only at N points.

Let's denote points where p(x) values was lost as z[j], j=1..M. So, we have:
```
p(x[1..N],z[1..M]) = (y[1..N], ?[1..M])
```
We will use the following polynom g(x) = (x-z[1])*....*(x-z[M]). It will have the following values at the same points:
```
g(x[1..N],z[1..M]) = (g[1..N], 0,0..0)
```
And the h(x) = p(x)*g(x) product has the following values:
```
h(x[1..N],z[1..M]) = (y[1]*g[1], ... y[N]*g[N], 0,0..0)
```

As you see, we know h(x) values at N+M points. We also know that it has order<N+M, because g(x) has order M and p(x) has order<N. So, now we can use FFT of order N+M to convert values at those N+M points into coefficients of h(x) polynom.

Once we have computed h(x), we can employ usual "polynom division" algorithm to compute p(x) = h(x) / g(x) in the O(N*logN) time. As result, we will get all coefficients of p(x) and then call FFT again to compute p(x) values corresponding to the source data
