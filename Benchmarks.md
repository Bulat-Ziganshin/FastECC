
All tests are performed on i7-4770 with 2-channel DDR3-1600 memory, employing all CPU cores.
Speeds are measured in MiB/s (mebibytes/second), add 5% to convert into MB/s (megabytes/second).

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
rs64g-avx2:  1766 ms = 1162 MiB/s,  cpu 12932 ms = 732%,  os 31 ms
rs64g-sse2:  2354 ms = 872 MiB/s,  cpu 16677 ms = 708%,  os 62 ms
rs64g:       3165 ms = 648 MiB/s,  cpu 23728 ms = 750%,  os 31 ms
rs64m:       3289 ms = 624 MiB/s,  cpu 25381 ms = 772%,  os 390 ms

rs32g-avx2:  1975 ms = 1039 MiB/s,  cpu 14758 ms = 747%,  os 16 ms
rs32g-sse2:  2490 ms = 824 MiB/s,  cpu 18283 ms = 734%,  os 16 ms
rs32g:       4661 ms = 440 MiB/s,  cpu 35755 ms = 767%,  os 16 ms
rs32m:       5991 ms = 342 MiB/s,  cpu 47050 ms = 785%,  os 187 ms
```

### Number-theoretic transform

MFA NTT (2^19 blocks of 2052 bytes each) in GF(0xFFF00001):
```
ntt64g-avx2:  770 ms = 1333 MiB/s,  cpu 5460 ms = 710%,  os 31 ms
ntt64g-sse2:  1064 ms = 964 MiB/s,  cpu 7644 ms = 718%,  os 16 ms
ntt64g:       1520 ms = 675 MiB/s,  cpu 11014 ms = 724%,  os 0 ms
ntt64m:       1594 ms = 644 MiB/s,  cpu 12043 ms = 756%,  os 406 ms

ntt32g-avx2:  802 ms = 1279 MiB/s,  cpu 5725 ms = 714%,  os 0 ms
ntt32g-sse2:  1077 ms = 953 MiB/s,  cpu 7878 ms = 732%,  os 16 ms
ntt32g:       2463 ms = 417 MiB/s,  cpu 18112 ms = 735%,  os 0 ms
ntt32m:       2679 ms = 383 MiB/s,  cpu 20811 ms = 777%,  os 109 ms
```

### Butterfly

Raw Butterfly operation in GF(0xFFF00001) processing 20 GiB of data.
NTT(2^N) requires about N-1 Butterfly operations per element,
so ideal NTT(2^19) implementation should be about 18x slower than the raw Butterfly:
```
ntt64g-avx2 b:  778 ms = 24528 MiB/s,  cpu 5491 ms = 706%,  os 0 ms
ntt64g-sse2 b:  1095 ms = 17423 MiB/s,  cpu 8502 ms = 777%,  os 0 ms
ntt64g b:       1733 ms = 11006 MiB/s,  cpu 13666 ms = 789%,  os 0 ms
ntt64m b:       2058 ms = 9267 MiB/s,  cpu 15897 ms = 772%,  os 250 ms

ntt32g-avx2 b:  733 ms = 26036 MiB/s,  cpu 5756 ms = 786%,  os 0 ms
ntt32g-sse2 b:  1098 ms = 17368 MiB/s,  cpu 8455 ms = 770%,  os 0 ms
ntt32g b:       3043 ms = 6267 MiB/s,  cpu 23868 ms = 784%,  os 0 ms
ntt32m b:       2835 ms = 6728 MiB/s,  cpu 22433 ms = 791%,  os 78 ms
```

### Lucky number?

MFA NTT, recursive NTT and raw Butterfly operations in various GF(P) and rings:
```
ntt64g-avx2 n 16 8192:  MFA_NTT<2^16,8192,P=0xFFF00001>: 374 ms = 1368 MiB/s,  cpu 2168 ms = 579%,  os 0 ms
ntt64g-avx2 o 16 8192:  Rec_NTT<2^16,8192,P=0xFFF00001>: 848 ms = 604 MiB/s,  cpu 1716 ms = 202%,  os 0 ms
ntt64g-avx2 b:          Butterfly: 802 ms = 23794 MiB/s,  cpu 5476 ms = 683%,  os 0 ms

ntt64g-sse2 n 16 8192:  MFA_NTT<2^16,8192,P=0xFFF00001>: 462 ms = 1107 MiB/s,  cpu 3120 ms = 675%,  os 0 ms
ntt64g-sse2 o 16 8192:  Rec_NTT<2^16,8192,P=0xFFF00001>: 986 ms = 519 MiB/s,  cpu 1981 ms = 201%,  os 0 ms
ntt64g-sse2 b:          Butterfly: 1129 ms = 16897 MiB/s,  cpu 7940 ms = 703%,  os 0 ms

ntt64g n 16 8192:       MFA_NTT<2^16,8192,P=0xFFF00001>: 659 ms = 777 MiB/s,  cpu 4477 ms = 680%,  os 0 ms
ntt64g o 16 8192:       Rec_NTT<2^16,8192,P=0xFFF00001>: 1570 ms = 326 MiB/s,  cpu 3151 ms = 201%,  os 0 ms
ntt64g b:               Butterfly: 1769 ms = 10782 MiB/s,  cpu 13166 ms = 744%,  os 0 ms

ntt64m n 16 8192:       MFA_NTT<2^16,8192,P=0xFFF00001>: 716 ms = 715 MiB/s,  cpu 5210 ms = 728%,  os 187 ms
ntt64m o 16 8192:       Rec_NTT<2^16,8192,P=0xFFF00001>: 1545 ms = 331 MiB/s,  cpu 6770 ms = 438%,  os 437 ms
ntt64m b:               Butterfly: 2117 ms = 9008 MiB/s,  cpu 16021 ms = 757%,  os 296 ms



ntt32g-avx2 n 16 8192:  MFA_NTT<2^16,8192,P=0xFFF00001>: 373 ms = 1373 MiB/s,  cpu 2168 ms = 582%,  os 0 ms
ntt32g-avx2 o 16 8192:  Rec_NTT<2^16,8192,P=0xFFF00001>: 865 ms = 592 MiB/s,  cpu 1794 ms = 207%,  os 0 ms
ntt32g-avx2 b:          Butterfly: 777 ms = 24550 MiB/s,  cpu 5444 ms = 701%,  os 0 ms

ntt32g-sse2 n 16 8192:  MFA_NTT<2^16,8192,P=0xFFF00001>: 471 ms = 1086 MiB/s,  cpu 3011 ms = 639%,  os 0 ms
ntt32g-sse2 o 16 8192:  Rec_NTT<2^16,8192,P=0xFFF00001>: 975 ms = 525 MiB/s,  cpu 1981 ms = 203%,  os 0 ms
ntt32g-sse2 b:          Butterfly: 1108 ms = 17216 MiB/s,  cpu 8128 ms = 734%,  os 0 ms

ntt32g n 16 8192:       MFA_NTT<2^16,8192,P=0xFFF00001>: 998 ms = 513 MiB/s,  cpu 7316 ms = 733%,  os 0 ms
ntt32g o 16 8192:       Rec_NTT<2^16,8192,P=0xFFF00001>: 2099 ms = 244 MiB/s,  cpu 4274 ms = 204%,  os 0 ms
ntt32g b:               Butterfly: 3063 ms = 6226 MiB/s,  cpu 23837 ms = 778%,  os 0 ms

ntt32m n 16 8192:       MFA_NTT<2^16,8192,P=0xFFF00001>: 1165 ms = 440 MiB/s,  cpu 8611 ms = 739%,  os 125 ms
ntt32m o 16 8192:       Rec_NTT<2^16,8192,P=0xFFF00001>: 2271 ms = 225 MiB/s,  cpu 9407 ms = 414%,  os 265 ms
ntt32m b:               Butterfly: 2918 ms = 6536 MiB/s,  cpu 22745 ms = 779%,  os 47 ms

------------------------------------------------------------------------------------------------

ntt64g-avx2 =n 16 8192:  MFA_NTT<2^16,8192,P=65537>: 261 ms = 982 MiB/s,  cpu 1498 ms = 575%,  os 0 ms
ntt64g-avx2 =o 16 8192:  Rec_NTT<2^16,8192,P=65537>: 631 ms = 406 MiB/s,  cpu 1279 ms = 203%,  os 0 ms
ntt64g-avx2 =b:          Butterfly: 339 ms = 28118 MiB/s,  cpu 1778 ms = 524%,  os 0 ms

ntt64g-sse2 =n 16 8192:  MFA_NTT<2^16,8192,P=65537>: 326 ms = 786 MiB/s,  cpu 1825 ms = 560%,  os 16 ms
ntt64g-sse2 =o 16 8192:  Rec_NTT<2^16,8192,P=65537>: 757 ms = 338 MiB/s,  cpu 1544 ms = 204%,  os 0 ms
ntt64g-sse2 =b:          Butterfly: 660 ms = 14460 MiB/s,  cpu 4228 ms = 641%,  os 0 ms

ntt64g =n 16 8192:       MFA_NTT<2^16,8192,P=65537>: 337 ms = 760 MiB/s,  cpu 1950 ms = 579%,  os 0 ms
ntt64g =o 16 8192:       Rec_NTT<2^16,8192,P=65537>: 757 ms = 338 MiB/s,  cpu 1544 ms = 204%,  os 0 ms
ntt64g =b:               Butterfly: 635 ms = 15008 MiB/s,  cpu 4384 ms = 690%,  os 0 ms

ntt64m =n 16 8192:       MFA_NTT<2^16,8192,P=65537>: 667 ms = 384 MiB/s,  cpu 4649 ms = 697%,  os 312 ms
ntt64m =o 16 8192:       Rec_NTT<2^16,8192,P=65537>: 1469 ms = 174 MiB/s,  cpu 6802 ms = 463%,  os 484 ms
ntt64m =b:               Butterfly: 1885 ms = 5058 MiB/s,  cpu 14305 ms = 759%,  os 250 ms



ntt32g-avx2 =n 16 8192:  MFA_NTT<2^16,8192,P=65537>: 628 ms = 407 MiB/s,  cpu 4384 ms = 698%,  os 0 ms
ntt32g-avx2 =o 16 8192:  Rec_NTT<2^16,8192,P=65537>: 1484 ms = 172 MiB/s,  cpu 2917 ms = 197%,  os 0 ms
ntt32g-avx2 =b:          Butterfly: 1706 ms = 5591 MiB/s,  cpu 12761 ms = 748%,  os 0 ms

ntt32g-sse2 =n 16 8192:  MFA_NTT<2^16,8192,P=65537>: 637 ms = 402 MiB/s,  cpu 4337 ms = 681%,  os 0 ms
ntt32g-sse2 =o 16 8192:  Rec_NTT<2^16,8192,P=65537>: 1506 ms = 170 MiB/s,  cpu 2964 ms = 197%,  os 0 ms
ntt32g-sse2 =b:          Butterfly: 1691 ms = 5639 MiB/s,  cpu 12761 ms = 755%,  os 0 ms

ntt32g =n 16 8192:       MFA_NTT<2^16,8192,P=65537>: 718 ms = 356 MiB/s,  cpu 4992 ms = 695%,  os 0 ms
ntt32g =o 16 8192:       Rec_NTT<2^16,8192,P=65537>: 1530 ms = 167 MiB/s,  cpu 3104 ms = 203%,  os 0 ms
ntt32g =b:               Butterfly: 1681 ms = 5672 MiB/s,  cpu 12823 ms = 763%,  os 0 ms

ntt32m =n 16 8192:       MFA_NTT<2^16,8192,P=65537>: 663 ms = 386 MiB/s,  cpu 4696 ms = 708%,  os 94 ms
ntt32m =o 16 8192:       Rec_NTT<2^16,8192,P=65537>: 1476 ms = 173 MiB/s,  cpu 6989 ms = 473%,  os 187 ms
ntt32m =b:               Butterfly: 1505 ms = 6335 MiB/s,  cpu 11357 ms = 754%,  os 203 ms

------------------------------------------------------------------------------------------------

ntt64g-avx2 -n 16 8192:  MFA_NTT<2^16,8192,P=2^32-1>: 472 ms = 1084 MiB/s,  cpu 3073 ms = 651%,  os 0 ms
ntt64g-avx2 -o 16 8192:  Rec_NTT<2^16,8192,P=2^32-1>: 1071 ms = 478 MiB/s,  cpu 2153 ms = 201%,  os 0 ms
ntt64g-avx2 -b:          Butterfly: 1123 ms = 16988 MiB/s,  cpu 8050 ms = 717%,  os 0 ms

ntt64g-sse2 -n 16 8192:  MFA_NTT<2^16,8192,P=2^32-1>: 374 ms = 1369 MiB/s,  cpu 2465 ms = 659%,  os 16 ms
ntt64g-sse2 -o 16 8192:  Rec_NTT<2^16,8192,P=2^32-1>: 834 ms = 614 MiB/s,  cpu 1732 ms = 208%,  os 0 ms
ntt64g-sse2 -b:          Butterfly: 805 ms = 23696 MiB/s,  cpu 5741 ms = 713%,  os 0 ms

ntt64g -n 16 8192:       MFA_NTT<2^16,8192,P=2^32-1>: 379 ms = 1349 MiB/s,  cpu 2200 ms = 580%,  os 0 ms
ntt64g -o 16 8192:       Rec_NTT<2^16,8192,P=2^32-1>: 820 ms = 624 MiB/s,  cpu 1638 ms = 200%,  os 0 ms
ntt64g -b:               Butterfly: 829 ms = 23008 MiB/s,  cpu 5834 ms = 704%,  os 0 ms

ntt64m -n 16 8192:       MFA_NTT<2^16,8192,P=2^32-1>: 558 ms = 918 MiB/s,  cpu 3838 ms = 688%,  os 234 ms
ntt64m -o 16 8192:       Rec_NTT<2^16,8192,P=2^32-1>: 1352 ms = 379 MiB/s,  cpu 6505 ms = 481%,  os 671 ms
ntt64m -b:               Butterfly: 1344 ms = 14194 MiB/s,  cpu 10202 ms = 759%,  os 203 ms



ntt32g-avx2 -n 16 8192:  MFA_NTT<2^16,8192,P=2^32-1>: 491 ms = 1043 MiB/s,  cpu 2948 ms = 601%,  os 0 ms
ntt32g-avx2 -o 16 8192:  Rec_NTT<2^16,8192,P=2^32-1>: 1085 ms = 472 MiB/s,  cpu 2168 ms = 200%,  os 0 ms
ntt32g-avx2 -b:          Butterfly: 1131 ms = 16857 MiB/s,  cpu 8003 ms = 707%,  os 0 ms

ntt32g-sse2 -n 16 8192:  MFA_NTT<2^16,8192,P=2^32-1>: 369 ms = 1387 MiB/s,  cpu 2215 ms = 600%,  os 0 ms
ntt32g-sse2 -o 16 8192:  Rec_NTT<2^16,8192,P=2^32-1>: 822 ms = 623 MiB/s,  cpu 1654 ms = 201%,  os 0 ms
ntt32g-sse2 -b:          Butterfly: 843 ms = 22635 MiB/s,  cpu 5694 ms = 676%,  os 0 ms

ntt32g -n 16 8192:       MFA_NTT<2^16,8192,P=2^32-1>: 764 ms = 670 MiB/s,  cpu 5195 ms = 680%,  os 0 ms
ntt32g -o 16 8192:       Rec_NTT<2^16,8192,P=2^32-1>: 1524 ms = 336 MiB/s,  cpu 3167 ms = 208%,  os 0 ms
ntt32g -b:               Butterfly: 1617 ms = 11796 MiB/s,  cpu 12324 ms = 762%,  os 0 ms

ntt32m -n 16 8192:       MFA_NTT<2^16,8192,P=2^32-1>: 692 ms = 740 MiB/s,  cpu 4976 ms = 720%,  os 109 ms
ntt32m -o 16 8192:       Rec_NTT<2^16,8192,P=2^32-1>: 1530 ms = 335 MiB/s,  cpu 7098 ms = 464%,  os 250 ms
ntt32m -b:               Butterfly: 1250 ms = 15258 MiB/s,  cpu 9563 ms = 765%,  os 62 ms

------------------------------------------------------------------------------------------------

ntt64g-avx2 +n 16 8192:  MFA_NTT<2^16,8192,P=2^64-1>: 306 ms = 1674 MiB/s,  cpu 1700 ms = 556%,  os 16 ms
ntt64g-avx2 +o 16 8192:  Rec_NTT<2^16,8192,P=2^64-1>: 709 ms = 722 MiB/s,  cpu 1466 ms = 207%,  os 0 ms
ntt64g-avx2 +b:          Butterfly: 601 ms = 31746 MiB/s,  cpu 3510 ms = 584%,  os 0 ms

ntt64g-sse2 +n 16 8192:  MFA_NTT<2^16,8192,P=2^64-1>: 323 ms = 1587 MiB/s,  cpu 1654 ms = 513%,  os 0 ms
ntt64g-sse2 +o 16 8192:  Rec_NTT<2^16,8192,P=2^64-1>: 710 ms = 721 MiB/s,  cpu 1466 ms = 206%,  os 16 ms
ntt64g-sse2 +b:          Butterfly: 583 ms = 32737 MiB/s,  cpu 3838 ms = 659%,  os 0 ms

ntt64g +n 16 8192:       MFA_NTT<2^16,8192,P=2^64-1>: 308 ms = 1665 MiB/s,  cpu 1716 ms = 558%,  os 31 ms
ntt64g +o 16 8192:       Rec_NTT<2^16,8192,P=2^64-1>: 722 ms = 709 MiB/s,  cpu 1544 ms = 214%,  os 0 ms
ntt64g +b:               Butterfly: 593 ms = 32185 MiB/s,  cpu 3853 ms = 650%,  os 0 ms

ntt64m +n 16 8192:       MFA_NTT<2^16,8192,P=2^64-1>: 298 ms = 1719 MiB/s,  cpu 1825 ms = 613%,  os 265 ms
ntt64m +o 16 8192:       Rec_NTT<2^16,8192,P=2^64-1>: 1036 ms = 494 MiB/s,  cpu 5928 ms = 572%,  os 593 ms
ntt64m +b:               Butterfly: 529 ms = 36051 MiB/s,  cpu 3604 ms = 681%,  os 156 ms
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
