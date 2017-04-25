FastECC provides O(N*log(N)) Reed-Solomon coder, running at 600 MB/s on i7-4770 with 2^20 blocks.
Version 0.1 implements only encoding, so it isn't yet ready for real use.

Almost all existing Reed-Solomon ECC implementations have O(N^2) speed behavior,
i.e. they can produce N ECC blocks in O(N^2) time, i.e. spend O(N) time per block.
F.e. the fastest implementation I know, MultiPar2, can compute 1000 ECC blocks at the speed ~50MB/s,
but only at ~2 MB/s in its maximum configuration, 32000 ECC blocks.
And 32-bit GF, implemented with the same approach, will compute one million ECC blocks at 50 KB/s.

The only exception is RSC32 by persicum, that has O(N*log(N)) speed,
i.e. it can produce ECC data in O(log(N)) time per ECC block.
Its speed with million ECC blocks is 100 MB/s,
i.e. it computes one million of 4 KB ECC blocks from one million of source blocks
(processing 8 GB overall) in just 80 seconds.
Note that all speeds mentioned here are measured on i7-4770, employing all features available in any particular program -
i.e. multi-threading, SIMD and x64 version.

- [NTT: Number-theoretic transform](Overview.md): what one needs to know in order to implement O(N*log(N)) Reed-Solomon error-correction codes
- [GF(p).cpp: fast computations in integer rings and fields](GF.md)
- [NTT.cpp: NTT implementation plus benchmarks](NTT.md)
- [RS.cpp: Reed-Solomon coder](RS.md)
