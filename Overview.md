What one needs to know in order to implement O(N*log(N)) Reed-Solomon error-correcting codes


# NTT implementation

Prerequsites:
* [Complex arithmetic](https://en.wikipedia.org/wiki/Complex_number), in particular add/sub/mul, powers, roots and polar form
* [Finite fields](https://en.wikipedia.org/wiki/Finite_field), with the same topics in mind

Topics:
* [Discrete Fourier transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform)
* [Fast Fourier transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform) and in particular
[Cooley–Tukey FFT algorithm](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm) as O(N*log(N)) algorithm implementing DFT
* [Fast Number-Theoretic Transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform_(general)) as modified FFT
employing the same add/sub/mul operations and unity roots, but in Galois Field

Once you grasped all these topics, you can grab some FFT implementation and convert it to NTT.


# NTT meaning

If source vector `A = (A0, A1...A[n-1])` represents polynomial `A(x) = A0 + A1*x + ... + A[n-1]*x**(n-1)`,
then NTT(A,a) computes polynomial values at the points `1, a, a**2 ... a**(n-1)`,
representing `n` roots of power `n` of 1 (if `a` is a primitive root of power `n` of 1).
It's equivalent to multiplying vector `A` to [Vandermonde Matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix) defined by these points.
So, NTT converts polynomial coefficients to polynomial values at these specific points, and iNTT (inverse NTT) reverses this transformation.
iNTT is almost the same as NTT, i.e. `iNTT(A,a) = NTT(A,1/a)/n` (with an element-wise division by n) reverses the effect of `NTT(A,a)`.

This means that `NTT(A,a)` provides the fast way to convert polynomial coefficients `A` into values at `a**i` points and vice versa.


# RS meaning

The idea behind Reed-Solomon `n+k` encoding (i.e. with `n` source words and `k` ecc words) is simple:
multiply source vector `A` to some fixed matrix `M[n*k]` and get resulting vector `C`.
Decoding needs a mutiplication of combined vector (A,C) to inverse of combined matrix `(I,M)` where `I[n*n]` is identity matrix.
Of course, the decoding required only when we lost some values in `A`, so real decoding combines alive values in `A` and `C` to vector of size `n`
(if more values are alive, they are just ignored), then deletes from `(I,M)` combined matrix the `k` lines corresponding to deleted values
and then mutiplies partial `(A,C)` vector of size `n` to inverse of partial `(I,M)` matrix of size `n*n`.
This operation computes the original `A` vector. And of course, all computations are performed in some Galois Field to ensure
that all the intermediate results still fit to some 8-32 bit word.

Therefore, all that we need are guarantees that any `n` lines of the `(I,M)` combined matrix form a invertible matrix.
If that's true, then we have an [MDS code](https://en.wikipedia.org/wiki/Erasure_code#Optimal_erasure_codes) - we can lose any `k` words
and restore the original data from `n` remaining ones.

There are so many ways to ensure that matrices are invertible:
* We can use [non-systematic code](https://en.wikipedia.org/wiki/Systematic_code), i.e. code where all `n+k` codes are computed
rather than include `n` original and `k` computed codes. Non-systematic codes are ok, for example, for data transfer over noisy channels.
So we just multiply source `A` vector by Vandermonde `(n+k)*n` matrix generated from `n+k` different `a[i]` numbers.
It's guaranteed that any `n` different `a[i]` numbers form an invertible Vandermonde matrix, so we can restore from any `n` remaining words after a loss.
* [Plank proposed](http://web.eecs.utk.edu/~plank/plank/papers/SPE-04.html) to start with Vandermonde `(n+k)*n` matrix
and then apply the [Gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination) in order to convert it to some `(I,M)` matrix.
As far as we perform this operation only once per a lot of parity computations, we can ignore the time required by this operation.
* PAR2 format employs `(I,V)` encoding matrix, i.e. it employs Vandermonde `k*n` matrix to compute `k` ecc words while employing the systematic code.
Despite of special form of `a[i]` used in their Vandermonde matrix, the restoration matrix is sometimes non-invertible.
But it seems to be a good compromise between the speed/complexity of computations and recovery strength.
* Persicum proposed us to use the [Cauchy matrix](https://en.wikipedia.org/wiki/Cauchy_matrix) of the form `M[i,j] = 1/(n+i+j)` (and it probably was empoyed by RAR5).
He said that it guarantees invertability of any `n*n` submatrix of `(I,M)`, but i have no idea whether it's true.


# NTT for RS encoding and decoding in O(N*log(N)) time

The first and most obvious idea to perform RS encoding in O(N*log(N)) time is to compute non-systematic code:
extend `n` source words with zeroes to `n+k`, NTT them and send `n+k` words produced to the channel. We are done!

This can be written as `C = V(a)*ext(A)` where `ext(A)` is a source vector `A` extended with zeroes to `n+k` elements,
`C` is encoded vector to be sent to the channel, `V(a)` is Vandermonde `(n+k)*(n+k)` matix produced by powers of `a`,
i.e. `(1, a, a**2 ... a**(n+k-1))` and `a` is a primitive root of 1 of power `n+k`.

Converting this to systematic code means the equation `(A,C) = V(a)*X`, so `A` becomes a part of encoded word
and we need to compute only remaining `C` part, and `X` is some unknown vector we have to compute.
