
All tests are performed on i7-4770 with 2-channel DDR3-1600 memory, employing all CPU cores.
Speeds are measured in MiB/s (mebibytes/second), which is 5% higher than MB/s (megabytes/second).

Executables are compiled by (-DSIMD selects vectorizable code path):
- `*64g-avx2`: 64-bit GCC 6.3 with -DSIMD=AVX2 -mavx2
- `*64g-sse2`: 64-bit GCC 6.3 with -DSIMD=SSE2
- `*64g`: 64-bit GCC 6.3
- `*64m`: 64-bit MSVC 2017
- `*32g-avx2`: 32-bit GCC 6.3 with -DSIMD=AVX2 -mavx2
- `*32g-sse2`: 32-bit GCC 6.3 with -DSIMD=SSE2 -msse2
- `*32g`: 32-bit GCC 6.3
- `*32m`: 32-bit MSVC 2017 with -arch:SSE2

NOTE: currently, executables built by MSVC/ICL/Clang are slower than GCC-compiled ones.
In the final version, main loops will be implemented with SSE2/AVX2 intrinsincs and have the same performance, irrespective of compiler.
For NTT(2^20), we expect speed of 1 GB/s for SSE2 version, and 2 GB/s for AVX2 version.


### Reed-Solomon encoding

Reed-Solomon encoding (2^19 source blocks => 2^19 ECC blocks, 2052 bytes each) in GF(0xFFF00001):
```
rs64g-avx2:  2154 ms = 952 MiB/s,  cpu 15538 ms = 721%,  os 62 ms
rs64g-sse2:  2641 ms = 777 MiB/s,  cpu 18346 ms = 695%,  os 31 ms
rs64g:       3462 ms = 593 MiB/s,  cpu 24648 ms = 712%,  os 47 ms
rs64m:       3572 ms = 574 MiB/s,  cpu 27300 ms = 764%,  os 686 ms

rs32g-avx2:  2324 ms = 883 MiB/s,  cpu 17098 ms = 736%,  os 16 ms
rs32g-sse2:  2718 ms = 755 MiB/s,  cpu 19594 ms = 721%,  os 16 ms
rs32g:       6039 ms = 340 MiB/s,  cpu 44913 ms = 744%,  os 0 ms
rs32m:       5637 ms = 364 MiB/s,  cpu 43977 ms = 780%,  os 390 ms
```

### Number-theoretic transform

MFA NTT (2^19 blocks of 2052 bytes each) in GF(0xFFF00001):
```
ntt64g-avx2:  931 ms = 1102 MiB/s,  cpu 6536 ms = 702%,  os 31 ms
ntt64g-sse2:  1158 ms = 886 MiB/s,  cpu 8159 ms = 705%,  os 16 ms
ntt64g:       1608 ms = 638 MiB/s,  cpu 10811 ms = 673%,  os 0 ms
ntt64m:       1724 ms = 595 MiB/s,  cpu 12714 ms = 737%,  os 406 ms

ntt32g-avx2:  962 ms = 1067 MiB/s,  cpu 6661 ms = 693%,  os 0 ms
ntt32g-sse2:  1202 ms = 853 MiB/s,  cpu 8252 ms = 686%,  os 0 ms
ntt32g:       2473 ms = 415 MiB/s,  cpu 18486 ms = 747%,  os 0 ms
ntt32m:       3182 ms = 322 MiB/s,  cpu 24523 ms = 771%,  os 312 ms
```

### Butterfly

Raw Butterfly operation in GF(0xFFF00001) processing 20 GiB of data.
NTT(2^N) requires about N-1 Butterfly operations per element,
so ideal NTT(2^19) implementation should be about 18x slower than the raw Butterfly:
```
ntt64g-avx2 b:  783 ms = 24370 MiB/s,  cpu 5460 ms = 698%,  os 0 ms
ntt64g-sse2 b:  1096 ms = 17402 MiB/s,  cpu 8611 ms = 786%,  os 0 ms
ntt64g b:       1745 ms = 10932 MiB/s,  cpu 13775 ms = 790%,  os 0 ms
ntt64m b:       2080 ms = 9170 MiB/s,  cpu 16115 ms = 775%,  os 156 ms

ntt32g-avx2 b:  733 ms = 26035 MiB/s,  cpu 5725 ms = 781%,  os 0 ms
ntt32g-sse2 b:  1109 ms = 17198 MiB/s,  cpu 8408 ms = 758%,  os 0 ms
ntt32g b:       3054 ms = 6246 MiB/s,  cpu 23291 ms = 763%,  os 0 ms
ntt32m b:       3131 ms = 6091 MiB/s,  cpu 24477 ms = 782%,  os 94 ms
```

### Lucky number?

MFA NTT, recursive NTT and raw Butterfly operations in various GF(P) and rings:
```
ntt64g n 16 8192:  MFA_NTT<2^16,8192,P=0xFFF00001>: 737 ms = 694 MiB/s,  cpu 4820 ms = 654%,  os 16 ms
ntt64g o 16 8192:  Rec_NTT<2^16,8192,P=0xFFF00001>: 1585 ms = 323 MiB/s,  cpu 3214 ms = 203%,  os 0 ms
ntt64g b:  Butterfly: 921 ms = 11120 MiB/s,  cpu 6521 ms = 708%,  os 0 ms

ntt64m n 16 8192:  MFA_NTT<2^16,8192,P=0xFFF00001>: 820 ms = 624 MiB/s,  cpu 5491 ms = 669%,  os 406 ms
ntt64m o 16 8192:  Rec_NTT<2^16,8192,P=0xFFF00001>: 1596 ms = 321 MiB/s,  cpu 7051 ms = 442%,  os 546 ms
ntt64m b:  Butterfly: 1162 ms = 8809 MiB/s,  cpu 8268 ms = 711%,  os 250 ms

ntt32g n 16 8192:  MFA_NTT<2^16,8192,P=0xFFF00001>: 1126 ms = 455 MiB/s,  cpu 7628 ms = 678%,  os 31 ms
ntt32g o 16 8192:  Rec_NTT<2^16,8192,P=0xFFF00001>: 2135 ms = 240 MiB/s,  cpu 4352 ms = 204%,  os 0 ms
ntt32g b:  Butterfly: 1578 ms = 6488 MiB/s,  cpu 11435 ms = 725%,  os 0 ms

ntt32m n 16 8192:  MFA_NTT<2^16,8192,P=0xFFF00001>: 1281 ms = 400 MiB/s,  cpu 9329 ms = 728%,  os 328 ms
ntt32m o 16 8192:  Rec_NTT<2^16,8192,P=0xFFF00001>: 2186 ms = 234 MiB/s,  cpu 8798 ms = 402%,  os 250 ms
ntt32m b:  Butterfly: 1592 ms = 6432 MiB/s,  cpu 11700 ms = 735%,  os 78 ms


ntt64g =n 16 8192:  MFA_NTT<2^16,8192,P=65537>: 512 ms = 500 MiB/s,  cpu 3323 ms = 649%,  os 0 ms
ntt64g =o 16 8192:  Rec_NTT<2^16,8192,P=65537>: 902 ms = 284 MiB/s,  cpu 1919 ms = 213%,  os 0 ms
ntt64g =b:  Butterfly: 364 ms = 28124 MiB/s,  cpu 1919 ms = 527%,  os 0 ms

ntt64m =n 16 8192:  MFA_NTT<2^16,8192,P=65537>: 791 ms = 324 MiB/s,  cpu 5538 ms = 700%,  os 312 ms
ntt64m =o 16 8192:  Rec_NTT<2^16,8192,P=65537>: 1494 ms = 171 MiB/s,  cpu 6677 ms = 447%,  os 624 ms
ntt64m =b:  Butterfly: 1016 ms = 10081 MiB/s,  cpu 7301 ms = 719%,  os 203 ms

ntt32g =n 16 8192:  MFA_NTT<2^16,8192,P=65537>: 811 ms = 316 MiB/s,  cpu 5678 ms = 700%,  os 0 ms
ntt32g =o 16 8192:  Rec_NTT<2^16,8192,P=65537>: 1602 ms = 160 MiB/s,  cpu 3276 ms = 205%,  os 0 ms
ntt32g =b:  Butterfly: 895 ms = 11440 MiB/s,  cpu 5897 ms = 659%,  os 0 ms

ntt32m =n 16 8192:  MFA_NTT<2^16,8192,P=65537>: 770 ms = 332 MiB/s,  cpu 5476 ms = 711%,  os 250 ms
ntt32m =o 16 8192:  Rec_NTT<2^16,8192,P=65537>: 1518 ms = 169 MiB/s,  cpu 7332 ms = 483%,  os 250 ms
ntt32m =b:  Butterfly: 813 ms = 12599 MiB/s,  cpu 5834 ms = 718%,  os 140 ms


ntt64g -n 16 8192:  MFA_NTT<2^16,8192,P=2^32-1>: 529 ms = 967 MiB/s,  cpu 3104 ms = 586%,  os 0 ms
ntt64g -o 16 8192:  Rec_NTT<2^16,8192,P=2^32-1>: 908 ms = 564 MiB/s,  cpu 1888 ms = 208%,  os 0 ms
ntt64g -b:  Butterfly: 459 ms = 22308 MiB/s,  cpu 2449 ms = 534%,  os 0 ms

ntt64m -n 16 8192:  MFA_NTT<2^16,8192,P=2^32-1>: 725 ms = 706 MiB/s,  cpu 5023 ms = 693%,  os 359 ms
ntt64m -o 16 8192:  Rec_NTT<2^16,8192,P=2^32-1>: 1361 ms = 376 MiB/s,  cpu 6646 ms = 488%,  os 515 ms
ntt64m -b:  Butterfly: 797 ms = 12852 MiB/s,  cpu 5444 ms = 683%,  os 281 ms

ntt32g -n 16 8192:  MFA_NTT<2^16,8192,P=2^32-1>: 861 ms = 595 MiB/s,  cpu 5460 ms = 634%,  os 0 ms
ntt32g -o 16 8192:  Rec_NTT<2^16,8192,P=2^32-1>: 1612 ms = 318 MiB/s,  cpu 3370 ms = 209%,  os 0 ms
ntt32g -b:  Butterfly: 843 ms = 12148 MiB/s,  cpu 5585 ms = 663%,  os 0 ms

ntt32m -n 16 8192:  MFA_NTT<2^16,8192,P=2^32-1>: 788 ms = 650 MiB/s,  cpu 5678 ms = 721%,  os 94 ms
ntt32m -o 16 8192:  Rec_NTT<2^16,8192,P=2^32-1>: 1554 ms = 330 MiB/s,  cpu 7160 ms = 461%,  os 312 ms
ntt32m -b:  Butterfly: 706 ms = 14502 MiB/s,  cpu 4976 ms = 705%,  os 47 ms


ntt64g +n 16 8192:  MFA_NTT<2^16,8192,P=2^64-1>: 517 ms = 991 MiB/s,  cpu 3619 ms = 701%,  os 16 ms
ntt64g +o 16 8192:  Rec_NTT<2^16,8192,P=2^64-1>: 881 ms = 581 MiB/s,  cpu 1888 ms = 214%,  os 0 ms
ntt64g +b:  Butterfly: 306 ms = 33443 MiB/s,  cpu 1778 ms = 581%,  os 0 ms

ntt64m +n 16 8192:  MFA_NTT<2^16,8192,P=2^64-1>: 498 ms = 1029 MiB/s,  cpu 3182 ms = 640%,  os 250 ms
ntt64m +o 16 8192:  Rec_NTT<2^16,8192,P=2^64-1>: 1082 ms = 473 MiB/s,  cpu 6271 ms = 580%,  os 484 ms
ntt64m +b:  Butterfly: 341 ms = 30065 MiB/s,  cpu 1919 ms = 563%,  os 172 ms
```

### Slow NTT

Comparison of slow O(N^2) NTT with fast algorithms:

```
C:\>timer ntt64g.exe s 20 32
Slow_NTT<2^20,32,P=0xFFF00001>: 9123729 ms = 0 MiB/s,  cpu 60703889 ms = 665%,  os 51808 ms
Verified! Original 2679569933,  after NTT: 1187104119
Kernel Time  =   101.041 = 00:01:41.041 =   0%
User Time    = 121412.286 = 33:43:32.286 = 665%
Process Time = 121513.328 = 33:45:13.328 = 665%
Global Time  = 18257.312 = 05:04:17.312 = 100%

C:\>timer ntt64g.exe o 20 32
Rec_NTT<2^20,32,P=0xFFF00001>: 322 ms = 99 MiB/s,  cpu 1264 ms = 393%,  os 0 ms
Verified!  Original 2679569933,  after NTT: 1187104119
Kernel Time  =     0.015 = 00:00:00.015 =   2%
User Time    =     1.778 = 00:00:01.778 = 308%
Process Time =     1.794 = 00:00:01.794 = 310%
Global Time  =     0.577 = 00:00:00.577 = 100%

C:\>timer ntt64g.exe n 20 32
MFA_NTT<2^20,32,P=0xFFF00001>: 124 ms = 258 MiB/s,  cpu 515 ms = 415%,  os 0 ms
Verified!  Original 2679569933,  after NTT: 1187104119
Kernel Time  =     0.000 = 00:00:00.000 =   0%
User Time    =     1.591 = 00:00:01.591 = 377%
Process Time =     1.591 = 00:00:01.591 = 377%
Global Time  =     0.421 = 00:00:00.421 = 100%
```
