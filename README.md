FastECC will provide O(N*log(N)) Reed-Solomon coder, running at the speed around 1 GB/s with 2^20 blocks.
Version 0.1 implements only encoding, so it isn't yet ready for real use.

Additional info:
- [NTT: Number-theoretic transform](NTT.md): what one need to know in order to implement N*logN Reed-Solomon error-correction codes
- [Benchmarks](bench.txt)

### Program usage

prime N - check whether N is prime, print divisors of N, and search, starting at N+1, for prime numbers as well as numbers only with small divisors

RS [N=19 [SIZE=1024]] - benchmark NTT-based Reed-Solomon encoding using 2^N input (source) blocks and 2^N output (ECC) blocks, each block SIZE bytes long

---

NTT [=|-][i|r|m|d|b|s|o|n] [N=20 [SIZE=512]] - test/benchmark GF(p) and NTT implementation

First argument is one of chars "irmdbson", optionally prefixed with "=" or "-" (character "n" may be omitted). Remaining arguments are used only for options "son".

By default, all computations are performed in GF(0xFFF00001). Prefix "=" switches to GF(0x10001), while prefix "-" switches to computations modulo 0xFFFFFFFF.
Note that 0xFFFFFFFF (2^32-1) isn't a prime number, nevertheless it supports NTT up to order 65536, and more than 2x faster than computations in GF(0xFFF00001).
Computations modulo 0xFFFFFFFF require normalisation (GF_Normalize call) after all computations.

The remainder of the first option is interpreted as following:
- i: test GF(p) implementation: check that each number in GF(p) has proper inverse (this check will fail for computations modulo 0xFFFFFFFF)
- m: test GF(p) implementation: check multiplication correctness (this check will also fail for computations modulo 0xFFFFFFFF since GF_Normalize isn't called here)
- r: find primary root of maximum order (P-1 for primary P, 65536 for P=0xFFFFFFFF)
- d: check divisors count and density, i.e. average "distance" to the next largest divider of the field order
- b: benchmark Butterfly operation (i.e. `a+b*K`) on 10 GB of input data (cosidered as 1.25G of (a,b) records). This is roughly equivalent to computing NTT(2^20) over 512 MB of data
- s: benchmark slow NTT (i.e. O(N^2) algo)
- o: benchmark old, recursive radix-2 NTT implementation
- n: benchmark new, faster MFA-based NTT implementation

NTT algorithms are performed using 2^N blocks SIZE bytes each. By default, N=20 and SIZE=512, these values can be overwritten in cmdline.
For every NTT, inverse operation is also performed and program verifies that NTT+iNTT results are equivalent to original data.


### Performance

Now the best possible performance of Reed-Solomon encoding is 600 MB/s on i7-4770 using ALL cores.
It can be reached with 2^16 source blocks and 2^16 ECC blocks, each 4 KB large,
i.e. encoding 256 MB of source data into 256 MB of ECC data, which will be processed in 0.8 seconds.

- Block size >= 4 KB. Smaller block sizes means more cache misses, it can be avoided only by careful prefetching.
- Overall data size limited to ~500 MB, processing larger datasizes are up to 1.5x slower. Efficient processing of larger data sizes will require
[Large Pages](https://msdn.microsoft.com/en-us/library/windows/desktop/aa366720(v=vs.85).aspx) support.
- Number of source blocks is 2^N. Current implementation supports only NTT of orders 2^N, so number of blocks is rounded up to the next 2^N value, thus making real performance up to 2x lower.
We need to implement PFA NTT as well as NTT kernels of orders 3,5,7,9,13 (since `0xFFF00000 = 2^20*3*3*5*7*13`) in order to get efficient support of arbitrary block counts.
0xFFF00000 has 504 dividers, so for random N the next divider of 0xFFF00000 is only a few percents larger than N itself.
- Number of ECC blocks (M) is equal to the number of source blocks (N). Current implementation performs backward transform (from N polynom coefficients to ECC block values)
using order-N NTT, even if M is much smaller than N, so backward transform always takes the same time as the forward one, making its effective speed N/M times lower.
But when M<=N/2, it's possible to compute only even points of the second transform, thus halving the execution time - we just need to perform a[i]+=a[i+N/2] and then run NTT on the first N/2 points.
When M<=N/4, we can compute only 1/4 of all points and so on, effectively running at the same effective speed as the first transform.
- Current MFA implementation is slightly inefficent. Ideally, it should split data into blocks of 512 KB or smaller (fitting into L3 cache even with m/t processing),
so the entire processing will require only 2-3 passes over memory.

If the final program version will implement all the features mentioned, it will run at 10/logb(N) GB/s with SSE2, and twice as fast with AVX2 - for ANY source+ECC configuration.

The speed can be further doubled by using computations modulo 2^32-1, but this ring supports only NTT of orders 2,4...65536.
Since this base is already supported by underlying GF(p).cpp library, required changes in RS.cpp are trivial - replace 0xFFF00001 with 0xFFFFFFFF
and post-process ECC blocks with GF_Normalize prior to writing.


