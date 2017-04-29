FastECC provides `O(N*log(N))` [Reed-Solomon coder], running at 1 GB/s on i7-4770 with 2^20 blocks.
Version 0.1 implements only encoding, so it isn't yet ready for real use.


<a name="what"/>

## What

Almost all existing Reed-Solomon ECC implementations have `O(N^2)` speed behavior,
i.e. they can produce N ECC blocks in `O(N^2)` time, thus spending `O(N)` time per block.
F.e. the fastest implementation I know, [MultiPar2], can compute 1000 ECC blocks at the speed ~50MB/s,
but only at ~2 MB/s in its maximum configuration, 32000 ECC blocks.
And 32-bit GF, implemented in the same way, will compute one million ECC blocks at 50 KB/s.

The only exception is [RSC32 by persicum] with `O(N*log(N))` speed,
i.e. it spends `O(log(N))` time per ECC block.
Its speed with million ECC blocks is 100 MB/s,
i.e. it computes one million of 4 KB ECC blocks from one million of source blocks
(processing 8 GB overall) in just 80 seconds.
Note that all speeds mentioned here are measured on i7-4770, employing all features available in a particular program -
including multi-threading, SIMD and x64 support.

FastECC is open-source library implementing `O(N*log(N))` encoding algorithm.
Depending on SIMD extension used, it's 4-10 times faster than RSC32, computing million ECC blocks at 400-1000 MB/s.
Future versions will implement decoding that's also `O(N*log(N))`, although 4-7 times slower than encoding.
Current implementation is limited to 2^20 blocks, removing this limit is the main priority for future work
aside of decoder implementation.


<a name="how"/>

## How

All `O(N*log(N))` Reed-Solomon implementations I'm aware of, use fast transforms like FFT or FWT.
FastECC employs Number-Theoretic Transform that is just an FFT over integer field or ring.
Let's see how it works. Note that below by `order-N polynomial` I mean any polynomial with order < N.

For any given set of N points, only one order-N polynomial may go through all these points.
Let's consider N input words as y-values of some order-N polynomial at N fixed x-points,
only one such polynomial may exist.

Typical Reed-Solomon encoding computes coefficients of this unique polynomial (the so-called [polynomial interpolation]),
evaluates the polynomial at M another fixed points (the `polynomial evaluation`)
and outputs these M words as the resulting ECC data.

At the decoding stage, we may receive any subset of N values out of N source words and M computed ECC words.
But since they all belong to the original order-N polynomial, we may recover this polynomial from N known points
and then compute its values in other points, in particular those N points assigned to original data, thus restoring them.


<a name="fast"/>

## Fast

Usually, Reed-Solomon libraries implement encoding by multiplication with [Vandermonde matrix] (`O(N^2)` algo)
and decoding by multiplication with the matrix inverse.

But with special choice of fixed x-points we can perform polynomial interpolation and evaluation at these points
in `O(N*log(N))` time, using NTT for evaluation and inverse NTT for interpolation. So, the fast encoding is as simple as:
- consider N input words as y-values of order-N polynomial in the special set of N x-points
- compute the polynomial coefficients in `O(N*log(N))` time using inverse NTT
- evaluate the polynomial at another M special x-points in `O(M*log(M))` time using NTT

Decoding is more involved. We have N words representing values of order-N polynomial at **some** N points.
Since we can't choose these points, we can't just use iNTT to compute the polynomial coefficients.
So it is a generic [polynomial interpolation] problem that can be solved in `O(N*log(N)^2)` time.
Or we can use customized interpolation algorithm (described below)
that can solve our specific problem in `O((N+M)*log(N+M))` time.
In the usual case of M<=N, the algorithm require three NTT(2*N) computations plus one NTT(N) computation,
making it 3.5-7 times slower than NTT(N)+NTT(M) performed on encoding.


<a name="more"/>

## More

- [NTT: Number-theoretic transform](Overview.md): what one needs to know in order to implement O(N*log(N)) Reed-Solomon error-correcting codes
- [GF(p).cpp: fast computations in integer rings and fields](GF.md)
- [NTT.cpp: NTT implementation plus benchmarks](NTT.md)
- [RS.cpp: Reed-Solomon coder](RS.md)


[Reed-Solomon coder]: https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction
[MultiPar2]: https://www.livebusinesschat.com/smf/index.php?board=396.0
[RSC32 by persicum]: https://www.livebusinesschat.com/smf/index.php?board=399.0
[polynomial interpolation]: https://en.wikipedia.org/wiki/Polynomial_interpolation
[Vandermonde matrix]: https://en.wikipedia.org/wiki/Vandermonde_matrix
