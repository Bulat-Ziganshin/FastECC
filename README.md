FastECC will provide O(N*log(N)) Reed-Solomon coder, running at the speed around 1 GB/s with 2^20 blocks.
Version 0.1 implements only encoding, so it isn't yet ready for real use.

- [NTT: Number-theoretic transform](Overview.md): what one needs to know in order to implement O(N*log(N)) Reed-Solomon error-correction codes
- [GF(p).cpp: fast computations in integer rings and fields](GF.md)
- [NTT.cpp: NTT implementation plus benchmarks](NTT.md)
- [RS.cpp: Reed-Solomon coder](RS.md)
