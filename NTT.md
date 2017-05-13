
### Program usage

`NTT [.][=-+][irmdbqson] [N=19 [SIZE=2052]]` - test/benchmark GF(p) and NTT implementations

First argument is one of chars "irmdbqson", optionally prefixed with "." for quiet mode and "=", "-" or "+" for GF(p) choice (character "n" may be omitted).
Remaining arguments are used only for options "qson".

By default, all computations are performed in GF(0xFFF00001). Prefix "=" switches to GF(0x10001),
while prefixes "-" and "+" switches to computations modulo 2^32-1 and 2^64-1, correspondingly.
Note that 2^32-1 and 2^64-1 aren't prime numbers, nevertheless they support NTT up to order 65536,
and with proper implementation more than 2x faster than computations in GF(0xFFF00001).
Computations modulo 2^32-1 and 2^64-1 require normalisation (GF_Normalize call) after all computations.

The remainder of the first option is interpreted as following:
- i: test GF(p) implementation: check that each number in GF(p) has proper inverse (this check will fail for computations modulo 2^n-1)
- m: test GF(p) implementation: check multiplication correctness (this check will also fail for computations modulo 2^n-1 since GF_Normalize isn't called here)
- r: find primary root of maximum order (P-1 for primary P, 65536 for P=2^32-1, `65536*2*5*17449` for P=2^64-1)
- d: check divisors count and density, i.e. average "distance" to the next largest divider of the field order
- b: benchmark Butterfly operation (i.e. `a+b*K`) on 20 GiB of input data (considered as 2.5Gi of (a,b) pairs). This is roughly equivalent to computing NTT(2^21) over 1 GiB of data,
but without overheads of NTT management - i.e. shows maximum NTT performance possible.
- q: benchmark slow quadratic NTT (i.e. O(N^2) algo)
- s: benchmark small NTT orders (run multiple times in single thread)
- o: benchmark old, recursive radix-2 NTT implementation
- n: benchmark new, faster MFA-based NTT implementation

NTT algorithms are performed using 2^N blocks SIZE bytes each. By default, N=19 and SIZE=2052 (N=5 for small NTT), other values can be specified as the second and third program options.
For every but small NTT, inverse operation is also performed and program verifies that NTT+iNTT results are equivalent to original data.
Incorrect results are reported like that:
```
Checksum mismatch: original 1690540224,  after NTT: 3386487444,  after NTT+iNTT 141226615
```


### Performance

Now the best possible performance of large-order Reed-Solomon encoding is 1.2 GB/s on i7-4770 using AVX2 and all cores.
It can be reached with 2^19 data blocks and 2^19 parity blocks, each 2 KiB large,
i.e. encoding 1 GiB of data into 1 GiB of parity, that is finished in 1.8 seconds.

With current implementation, maximum performance is reached only when all of the following conditions are met:
- Block size >= 2 KB. Smaller block sizes means more cache misses, it can be avoided only by careful prefetching.
- Number of source data blocks is 2^N. Current implementation supports only NTT of orders 2^N, so number of blocks is rounded up to the next 2^N value, thus making real performance up to 2x lower.
We need to implement PFA NTT as well as NTT kernels of orders 3,5,7,9,13 (since `0xFFF00000 = 2^20*3*3*5*7*13`) in order to get efficient support of arbitrary block counts.
0xFFF00000 has 504 dividers, so for random N the next divider of 0xFFF00000 is only a few percents larger than N itself.
- Number of parity blocks (M) is equal to the number of data blocks (N). Current implementation performs backward transform (from N polynomial coefficients to parity block values)
using order-N NTT, even if M is much smaller than N, so backward transform always takes the same time as the forward one, making its effective speed N/M times lower.
But when M<=N/2, it's possible to compute only even points of the second transform, thus halving the execution time - we just need to perform a[i]+=a[i+N/2] and then run NTT on the first N/2 points.
When M<=N/4, we can compute only 1/4 of all points and so on, effectively running at the same effective speed as the first transform.

If the final program version will implement all the features mentioned, it will run at ~20/logb(N) GB/s with SSE2, and twice as fast with AVX2 - for ANY data+parity configuration.

The speed can be further doubled by using computations modulo 2^32-1 or 2^64-1, but these rings support only NTT of orders 2,4...65536.
Since this base is already supported by underlying GF(p).cpp library, required changes in RS.cpp are trivial - replace 0xFFF00001 with 0xFFFFFFFF
and post-process parity blocks with GF_Normalize prior to writing.
