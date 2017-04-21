FastECC will provide O(N*log(N)) Reed-Solomon coder, running at the speed around 1 GB/s with 2^20 blocks.
Version 0.1 implements only encoding, so it isn't yet ready for real use.

Further reading:
- [NTT: Number-theoretic transform](Overview.md): what one needs to know in order to implement O(N*log(N)) Reed-Solomon error-correction codes
- [GF(p).cpp: fast computations in integer rings and fields](GF.md)
- [NTT.cpp: NTT implementation](NTT.md)
- [RS.cpp: Reed-Solomon coder](RS.md)
- [Benchmarks](bench.txt)


### Program usage

`prime N` - check whether N is prime, print divisors of N, and search, starting at N+1, for prime numbers as well as numbers only with small divisors
