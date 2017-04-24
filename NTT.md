
### Program usage

`NTT [=|-][i|r|m|d|b|s|o|n] [N=20 [SIZE=512]]` - test/benchmark GF(p) and NTT implementation

First argument is one of chars "irmdbson", optionally prefixed with "=", "-" or "+" (character "n" may be omitted). Remaining arguments are used only for options "son".

By default, all computations are performed in GF(0xFFF00001). Prefix "=" switches to GF(0x10001), 
while prefixes "-" and "+" switches to computations modulo 2^32-1 and 2^64-1, correspondingly.
Note that 2^32-1 and 2^64-1 aren't prime numbers, nevertheless they support NTT up to order 65536, and more than 2x faster than computations in GF(0xFFF00001).
Computations modulo 2^32-1 and 2^64-1 require normalisation (GF_Normalize call) after all computations.

The remainder of the first option is interpreted as following:
- i: test GF(p) implementation: check that each number in GF(p) has proper inverse (this check will fail for computations modulo 2^n-1)
- m: test GF(p) implementation: check multiplication correctness (this check will also fail for computations modulo 2^n-1 since GF_Normalize isn't called here)
- r: find primary root of maximum order (P-1 for primary P, 65536 for P=2^32-1)
- d: check divisors count and density, i.e. average "distance" to the next largest divider of the field order
- b: benchmark Butterfly operation (i.e. `a+b*K`) on 10 GB of input data (cosidered as 1.25G of (a,b) records). This is roughly equivalent to computing NTT(2^20) over 512 MB of data
- s: benchmark slow NTT (i.e. O(N^2) algo)
- o: benchmark old, recursive radix-2 NTT implementation
- n: benchmark new, faster MFA-based NTT implementation

NTT algorithms are performed using 2^N blocks SIZE bytes each. By default, N=20 and SIZE=512, these values can be overwritten in cmdline.
For every NTT, inverse operation is also performed and program verifies that NTT+iNTT results are equivalent to original data.


### Performance

Now the best possible performance of Reed-Solomon encoding is 600 MB/s on i7-4770 using ALL cores.
It can be reached with 2^20 source blocks and 2^20 ECC blocks, each 4 KB large,
i.e. encoding 2 GB of source data into 2 GB of ECC data, that is processed in 7 seconds.

With current implementation, maximum performance is reached only when all of the following conditions are met:
- Block size >= 4 KB. Smaller block sizes means more cache misses, it can be avoided only by careful prefetching.
- Number of source blocks is 2^N. Current implementation supports only NTT of orders 2^N, so number of blocks is rounded up to the next 2^N value, thus making real performance up to 2x lower.
We need to implement PFA NTT as well as NTT kernels of orders 3,5,7,9,13 (since `0xFFF00000 = 2^20*3*3*5*7*13`) in order to get efficient support of arbitrary block counts.
0xFFF00000 has 504 dividers, so for random N the next divider of 0xFFF00000 is only a few percents larger than N itself.
- Number of ECC blocks (M) is equal to the number of source blocks (N). Current implementation performs backward transform (from N polynom coefficients to ECC block values)
using order-N NTT, even if M is much smaller than N, so backward transform always takes the same time as the forward one, making its effective speed N/M times lower.
But when M<=N/2, it's possible to compute only even points of the second transform, thus halving the execution time - we just need to perform a[i]+=a[i+N/2] and then run NTT on the first N/2 points.
When M<=N/4, we can compute only 1/4 of all points and so on, effectively running at the same effective speed as the first transform.

If the final program version will implement all the features mentioned, it will run at 12/logb(N) GB/s with SSE2, and twice as fast with AVX2 - for ANY source+ECC configuration.

The speed can be further doubled by using computations modulo 2^32-1 or 2^64-1, but these rings support only NTT of orders 2,4...65536.
Since this base is already supported by underlying GF(p).cpp library, required changes in RS.cpp are trivial - replace 0xFFF00001 with 0xFFFFFFFF
and post-process ECC blocks with GF_Normalize prior to writing.


### Benchmarks

```
C:\!FreeArc\public\FastECC>timer ntt64g.exe s 20 32
Slow_NTT<2^20,32>: 9123729 ms = 0 MiB/s,  cpu 60703889 ms = 665%,  os 51808 ms
Verified! Original 2679569933,  after NTT: 1187104119

Kernel Time  =   101.041 = 00:01:41.041 =   0%
User Time    = 121412.286 = 33:43:32.286 = 665%
Process Time = 121513.328 = 33:45:13.328 = 665%
Global Time  = 18257.312 = 05:04:17.312 = 100%



C:\!FreeArc\public\FastECC>timer ntt64g.exe o 20 32
Rec_NTT<2^20,32>: 365 ms = 88 MiB/s,  cpu 1388 ms = 380%,  os 0 ms
Verified! Original 2679569933,  after NTT: 1187104119

Kernel Time  =     0.000 = 00:00:00.000 =   0%
User Time    =     1.794 = 00:00:01.794 = 302%
Process Time =     1.794 = 00:00:01.794 = 302%
Global Time  =     0.593 = 00:00:00.593 = 100%



C:\!FreeArc\public\FastECC>timer ntt64g.exe n 20 32
MFA_NTT<2^20,32>: 159 ms = 202 MiB/s,  cpu 671 ms = 423%,  os 0 ms
Checksum mismatch: original 2679569933,  after NTT: 2655366899,  after NTT+iNTT 489034790

Kernel Time  =     0.015 = 00:00:00.015 =   3%
User Time    =     1.700 = 00:00:01.700 = 363%
Process Time =     1.716 = 00:00:01.716 = 366%
Global Time  =     0.468 = 00:00:00.468 = 100%



GF(0xFFF00001):

C:\!FreeArc\public\FastECC>ntt64g.exe o 16 8192
Rec_NTT<2^16,8192>: 3467 ms = 148 MiB/s,  cpu 3432 ms = 99%,  os 0 ms
Verified!  Original 4053815590,  after NTT: 3321235535
C:\!FreeArc\public\FastECC>ntt64g.exe n 16 8192
MFA_NTT<2^16,8192>: 3122 ms = 164 MiB/s,  cpu 3089 ms = 99%,  os 0 ms
Verified!  Original 4053815590,  after NTT: 3321235535
C:\!FreeArc\public\FastECC>ntt64g.exe b
Butterfly: 3955 ms = 2589 MiB/s,  cpu 3931 ms = 99%,  os 0 ms

C:\!FreeArc\public\FastECC>ntt32g.exe o 16 8192
Rec_NTT<2^16,8192>: 5349 ms = 96 MiB/s,  cpu 5320 ms = 99%,  os 0 ms
Verified!  Original 4053815590,  after NTT: 3321235535
C:\!FreeArc\public\FastECC>ntt32g.exe n 16 8192
MFA_NTT<2^16,8192>: 4623 ms = 111 MiB/s,  cpu 4571 ms = 99%,  os 0 ms
Verified!  Original 4053815590,  after NTT: 3321235535
C:\!FreeArc\public\FastECC>ntt32g.exe b
Butterfly: 6136 ms = 1669 MiB/s,  cpu 6053 ms = 99%,  os 0 ms


GF(0x10001), real results are 2x slower since we process 32-bit words with essentially 16-bit values:

C:\!FreeArc\public\FastECC>ntt64g.exe =o 16 8192
Rec_NTT<2^16,8192>: 1380 ms = 371 MiB/s,  cpu 1388 ms = 101%,  os 0 ms
Verified!  Original 1807322693,  after NTT: 4193637393
C:\!FreeArc\public\FastECC>ntt64g.exe =n 16 8192
MFA_NTT<2^16,8192>: 1312 ms = 390 MiB/s,  cpu 1295 ms = 99%,  os 0 ms
Verified!  Original 1807322693,  after NTT: 4193637393
C:\!FreeArc\public\FastECC>ntt64g.exe =b
Butterfly: 1493 ms = 6857 MiB/s,  cpu 1482 ms = 99%,  os 0 ms              = 2.65x faster

C:\!FreeArc\public\FastECC>ntt32g.exe =o 16 8192
Rec_NTT<2^16,8192>: 3709 ms = 138 MiB/s,  cpu 3619 ms = 98%,  os 0 ms
Verified!  Original 1807322693,  after NTT: 4193637393
C:\!FreeArc\public\FastECC>ntt32g.exe =n 16 8192
MFA_NTT<2^16,8192>: 3691 ms = 139 MiB/s,  cpu 3604 ms = 98%,  os 0 ms
Verified!  Original 1807322693,  after NTT: 4193637393
C:\!FreeArc\public\FastECC>ntt32g.exe =b
Butterfly: 3832 ms = 2673 MiB/s,  cpu 3806 ms = 99%,  os 0 ms
```
