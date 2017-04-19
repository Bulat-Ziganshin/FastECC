More info:
- [NTT: Number-theoretic transform](NTT.md): what one need to know in order to implement N*logN Reed-Solomon error-correction codes
- [Benchmarks](bench.txt)

### Program usage

prime N - check whether N is prime, print divisors of N, and search, starting at N+1, for prime numbers as well as numbers only with small divisors

RS [N=19 [SIZE=1024]] - benchmark NTT-based Reed-Solomon encoding using 2^N input (source) blocks and 2^N output (ECC) blocks, each block SIZE bytes long

---

main [=|-][i|r|m|d|b|s|o|n] [N=20 [SIZE=512]] - test/benchmark GF(p) and NTT implementation

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
For every NTT, inverse operation is also performed and program verifies whether NTT+iNTT results are equivalent to original data.

