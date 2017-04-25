FastECC provides `O(N*log(N))` Reed-Solomon coder, running at 600 MB/s on i7-4770 with 2^20 blocks.
Version 0.1 implements only encoding, so it isn't yet ready for real use.


## What

Almost all existing Reed-Solomon ECC implementations have `O(N^2)` speed behavior,
i.e. they can produce N ECC blocks in `O(N^2)` time, thus spending `O(N)` time per block.
F.e. the fastest implementation I know, MultiPar2, can compute 1000 ECC blocks at the speed ~50MB/s,
but only at ~2 MB/s in its maximum configuration, 32000 ECC blocks.
And 32-bit GF, implemented with the same approach, will compute one million ECC blocks at 50 KB/s.

The only exception is RSC32 by persicum with `O(N*log(N))` speed, i.e. it spends `O(log(N))` time per ECC block.
Its speed with million ECC blocks is 100 MB/s,
i.e. it computes one million of 4 KB ECC blocks from one million of source blocks
(processing 8 GB overall) in just 80 seconds.
Note that all speeds mentioned here are measured on i7-4770, employing all features available in any particular program -
i.e. multi-threading, SIMD and x64 version.

FastECC is open-source library implementing `O(N*log(N))` encoding algorithm.
It's 6 times faster than RSC32, i.e. it computes million ECC blocks at 600 MB/s.
Future versions will implement decoding that's also `O(N*log(N))`, although 4-7 times slower than encoding.
Further optimization is possible, although is not my top priority.
Current implementation is limited to 2^20 blocks, removing this limit is main priority for future work,
aside of decoder implementation.


## How

All `O(N*log(N))` Reed-Solomon implementations I'm aware of, use fast transforms like FFT or FWT.
FastECC employs Number-Theoretic Transform that is just an FFT over integer field or ring.
Let's see how it works. Note that by `order-N polynoms` I mean polynoms with any order < N.

For any given set of N points, only one order-N polynom may go through all these points.
Let's consider N input words as values of some order-N polynom at N fixed points,
only one such polynom may exist.

Typical Reed-Solomon encoding computes coefficients of this unique polynom (the so-called `polynom interpolation`),
evaluates the polynom at M another fixed points (the `polynom evaluation`)
and outputs these M words as the resulting ECC data.

At the decoding stage, we may receive any subset of N values out of N source words and M computed ECC words.
But since they all belong to the original order-N polynom, we may recover this polynom from N known points
and then compute its values in other points, in particular those N points assigned to original data, thus restoring them.


## More

- [NTT: Number-theoretic transform](Overview.md): what one needs to know in order to implement O(N*log(N)) Reed-Solomon error-correcting codes
- [GF(p).cpp: fast computations in integer rings and fields](GF.md)
- [NTT.cpp: NTT implementation plus benchmarks](NTT.md)
- [RS.cpp: Reed-Solomon coder](RS.md)
