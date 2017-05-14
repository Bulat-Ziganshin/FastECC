
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

MFA NTT, recursive NTT and raw Butterfly operations in various GF(p) and rings:
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

### Small NTT

How NTT speed depends on the order and Ring base. These tests are single-threaded (`set OMP_NUM_THREADS=1`)

```
C:\>for %s in (20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt64g-avx2 s %s
MFA_NTT<2^20,2052,P=0xFFF00001>*1: 7325 ms = 280 MiB/s,  cpu 7301 ms = 100%,  os 16 ms
MFA_NTT<2^19,2052,P=0xFFF00001>*2: 7042 ms = 291 MiB/s,  cpu 7020 ms = 100%,  os 31 ms
MFA_NTT<2^18,2052,P=0xFFF00001>*4: 6717 ms = 306 MiB/s,  cpu 6724 ms = 100%,  os 0 ms
MFA_NTT<2^17,2052,P=0xFFF00001>*8: 6144 ms = 334 MiB/s,  cpu 6131 ms = 100%,  os 0 ms
MFA_NTT<2^16,2052,P=0xFFF00001>*16: 5753 ms = 357 MiB/s,  cpu 5741 ms = 100%,  os 16 ms
MFA_NTT<2^15,2052,P=0xFFF00001>*32: 5383 ms = 381 MiB/s,  cpu 5398 ms = 100%,  os 0 ms
MFA_NTT<2^14,2052,P=0xFFF00001>*64: 5033 ms = 408 MiB/s,  cpu 5008 ms = 99%,  os 16 ms
MFA_NTT<2^13,2052,P=0xFFF00001>*128: 4577 ms = 448 MiB/s,  cpu 4540 ms = 99%,  os 31 ms
MFA_NTT<2^12,2052,P=0xFFF00001>*256: 4057 ms = 506 MiB/s,  cpu 4056 ms = 100%,  os 0 ms
MFA_NTT<2^11,2052,P=0xFFF00001>*512: 3592 ms = 571 MiB/s,  cpu 3572 ms = 99%,  os 16 ms
MFA_NTT<2^10,2052,P=0xFFF00001>*1024: 3207 ms = 640 MiB/s,  cpu 3198 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=0xFFF00001>*2048: 2807 ms = 731 MiB/s,  cpu 2792 ms = 99%,  os 16 ms
MFA_NTT<2^8,2052,P=0xFFF00001>*4096: 2384 ms = 861 MiB/s,  cpu 2371 ms = 99%,  os 16 ms
MFA_NTT<2^7,2052,P=0xFFF00001>*8192: 2017 ms = 1017 MiB/s,  cpu 2012 ms = 100%,  os 0 ms
MFA_NTT<2^6,2052,P=0xFFF00001>*16384: 1734 ms = 1183 MiB/s,  cpu 1700 ms = 98%,  os 31 ms
MFA_NTT<2^5,2052,P=0xFFF00001>*32768: 1354 ms = 1515 MiB/s,  cpu 1342 ms = 99%,  os 0 ms
MFA_NTT<2^4,2052,P=0xFFF00001>*65536: 991 ms = 2070 MiB/s,  cpu 983 ms = 99%,  os 0 ms
MFA_NTT<2^3,2052,P=0xFFF00001>*131072: 667 ms = 3077 MiB/s,  cpu 655 ms = 98%,  os 0 ms
MFA_NTT<2^2,2052,P=0xFFF00001>*262144: 377 ms = 5447 MiB/s,  cpu 374 ms = 99%,  os 0 ms
MFA_NTT<2^1,2052,P=0xFFF00001>*524288: 168 ms = 12187 MiB/s,  cpu 172 ms = 102%,  os 0 ms

C:\>for %s in (20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt64g-sse2 s %s
MFA_NTT<2^20,2052,P=0xFFF00001>*1: 9470 ms = 217 MiB/s,  cpu 9469 ms = 100%,  os 0 ms
MFA_NTT<2^19,2052,P=0xFFF00001>*2: 9074 ms = 226 MiB/s,  cpu 9079 ms = 100%,  os 0 ms
MFA_NTT<2^18,2052,P=0xFFF00001>*4: 8565 ms = 240 MiB/s,  cpu 8533 ms = 100%,  os 47 ms
MFA_NTT<2^17,2052,P=0xFFF00001>*8: 8052 ms = 255 MiB/s,  cpu 8050 ms = 100%,  os 16 ms
MFA_NTT<2^16,2052,P=0xFFF00001>*16: 7453 ms = 275 MiB/s,  cpu 7457 ms = 100%,  os 0 ms
MFA_NTT<2^15,2052,P=0xFFF00001>*32: 6899 ms = 297 MiB/s,  cpu 6895 ms = 100%,  os 0 ms
MFA_NTT<2^14,2052,P=0xFFF00001>*64: 6446 ms = 318 MiB/s,  cpu 6443 ms = 100%,  os 16 ms
MFA_NTT<2^13,2052,P=0xFFF00001>*128: 5910 ms = 347 MiB/s,  cpu 5897 ms = 100%,  os 16 ms
MFA_NTT<2^12,2052,P=0xFFF00001>*256: 5335 ms = 385 MiB/s,  cpu 5273 ms = 99%,  os 47 ms
MFA_NTT<2^11,2052,P=0xFFF00001>*512: 4807 ms = 427 MiB/s,  cpu 4805 ms = 100%,  os 0 ms
MFA_NTT<2^10,2052,P=0xFFF00001>*1024: 4233 ms = 485 MiB/s,  cpu 4228 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=0xFFF00001>*2048: 3719 ms = 552 MiB/s,  cpu 3728 ms = 100%,  os 0 ms
MFA_NTT<2^8,2052,P=0xFFF00001>*4096: 3237 ms = 634 MiB/s,  cpu 3229 ms = 100%,  os 0 ms
MFA_NTT<2^7,2052,P=0xFFF00001>*8192: 2769 ms = 741 MiB/s,  cpu 2761 ms = 100%,  os 0 ms
MFA_NTT<2^6,2052,P=0xFFF00001>*16384: 2308 ms = 889 MiB/s,  cpu 2262 ms = 98%,  os 47 ms
MFA_NTT<2^5,2052,P=0xFFF00001>*32768: 1785 ms = 1150 MiB/s,  cpu 1778 ms = 100%,  os 0 ms
MFA_NTT<2^4,2052,P=0xFFF00001>*65536: 1316 ms = 1559 MiB/s,  cpu 1326 ms = 101%,  os 0 ms
MFA_NTT<2^3,2052,P=0xFFF00001>*131072: 881 ms = 2329 MiB/s,  cpu 889 ms = 101%,  os 0 ms
MFA_NTT<2^2,2052,P=0xFFF00001>*262144: 514 ms = 3992 MiB/s,  cpu 515 ms = 100%,  os 0 ms
MFA_NTT<2^1,2052,P=0xFFF00001>*524288: 249 ms = 8235 MiB/s,  cpu 265 ms = 106%,  os 0 ms



C:\>for %s in (20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt32g-avx2 s %s
MFA_NTT<2^20,2032,P=0xFFF00001>*1: 7750 ms = 262 MiB/s,  cpu 7722 ms = 100%,  os 16 ms
MFA_NTT<2^19,2052,P=0xFFF00001>*2: 7398 ms = 277 MiB/s,  cpu 7394 ms = 100%,  os 0 ms
MFA_NTT<2^18,2052,P=0xFFF00001>*4: 7002 ms = 293 MiB/s,  cpu 7004 ms = 100%,  os 16 ms
MFA_NTT<2^17,2052,P=0xFFF00001>*8: 6508 ms = 315 MiB/s,  cpu 6505 ms = 100%,  os 0 ms
MFA_NTT<2^16,2052,P=0xFFF00001>*16: 6054 ms = 339 MiB/s,  cpu 6053 ms = 100%,  os 0 ms
MFA_NTT<2^15,2052,P=0xFFF00001>*32: 5646 ms = 363 MiB/s,  cpu 5647 ms = 100%,  os 0 ms
MFA_NTT<2^14,2052,P=0xFFF00001>*64: 5193 ms = 395 MiB/s,  cpu 5179 ms = 100%,  os 0 ms
MFA_NTT<2^13,2052,P=0xFFF00001>*128: 4801 ms = 427 MiB/s,  cpu 4789 ms = 100%,  os 0 ms
MFA_NTT<2^12,2052,P=0xFFF00001>*256: 4280 ms = 479 MiB/s,  cpu 4274 ms = 100%,  os 0 ms
MFA_NTT<2^11,2052,P=0xFFF00001>*512: 3743 ms = 548 MiB/s,  cpu 3728 ms = 100%,  os 16 ms
MFA_NTT<2^10,2052,P=0xFFF00001>*1024: 3400 ms = 603 MiB/s,  cpu 3401 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=0xFFF00001>*2048: 2998 ms = 684 MiB/s,  cpu 2980 ms = 99%,  os 16 ms
MFA_NTT<2^8,2052,P=0xFFF00001>*4096: 2603 ms = 788 MiB/s,  cpu 2590 ms = 99%,  os 0 ms
MFA_NTT<2^7,2052,P=0xFFF00001>*8192: 2192 ms = 936 MiB/s,  cpu 2200 ms = 100%,  os 0 ms
MFA_NTT<2^6,2052,P=0xFFF00001>*16384: 1827 ms = 1123 MiB/s,  cpu 1794 ms = 98%,  os 31 ms
MFA_NTT<2^5,2052,P=0xFFF00001>*32768: 1471 ms = 1395 MiB/s,  cpu 1466 ms = 100%,  os 0 ms
MFA_NTT<2^4,2052,P=0xFFF00001>*65536: 1063 ms = 1931 MiB/s,  cpu 1061 ms = 100%,  os 0 ms
MFA_NTT<2^3,2052,P=0xFFF00001>*131072: 703 ms = 2920 MiB/s,  cpu 702 ms = 100%,  os 0 ms
MFA_NTT<2^2,2052,P=0xFFF00001>*262144: 409 ms = 5013 MiB/s,  cpu 406 ms = 99%,  os 0 ms
MFA_NTT<2^1,2052,P=0xFFF00001>*524288: 203 ms = 10094 MiB/s,  cpu 203 ms = 100%,  os 0 ms

C:\>for %s in (20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt32g-sse2 s %s
MFA_NTT<2^20,2032,P=0xFFF00001>*1: 9425 ms = 216 MiB/s,  cpu 9422 ms = 100%,  os 0 ms
MFA_NTT<2^19,2052,P=0xFFF00001>*2: 9120 ms = 225 MiB/s,  cpu 9095 ms = 100%,  os 16 ms
MFA_NTT<2^18,2052,P=0xFFF00001>*4: 8677 ms = 236 MiB/s,  cpu 8658 ms = 100%,  os 0 ms
MFA_NTT<2^17,2052,P=0xFFF00001>*8: 8109 ms = 253 MiB/s,  cpu 8081 ms = 100%,  os 31 ms
MFA_NTT<2^16,2052,P=0xFFF00001>*16: 7474 ms = 275 MiB/s,  cpu 7472 ms = 100%,  os 0 ms
MFA_NTT<2^15,2052,P=0xFFF00001>*32: 6969 ms = 294 MiB/s,  cpu 6973 ms = 100%,  os 0 ms
MFA_NTT<2^14,2052,P=0xFFF00001>*64: 6446 ms = 318 MiB/s,  cpu 6427 ms = 100%,  os 31 ms
MFA_NTT<2^13,2052,P=0xFFF00001>*128: 5944 ms = 345 MiB/s,  cpu 5928 ms = 100%,  os 0 ms
MFA_NTT<2^12,2052,P=0xFFF00001>*256: 5369 ms = 382 MiB/s,  cpu 5366 ms = 100%,  os 16 ms
MFA_NTT<2^11,2052,P=0xFFF00001>*512: 4849 ms = 423 MiB/s,  cpu 4820 ms = 99%,  os 0 ms
MFA_NTT<2^10,2052,P=0xFFF00001>*1024: 4262 ms = 481 MiB/s,  cpu 4259 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=0xFFF00001>*2048: 3774 ms = 544 MiB/s,  cpu 3791 ms = 100%,  os 0 ms
MFA_NTT<2^8,2052,P=0xFFF00001>*4096: 3313 ms = 619 MiB/s,  cpu 3292 ms = 99%,  os 16 ms
MFA_NTT<2^7,2052,P=0xFFF00001>*8192: 2870 ms = 715 MiB/s,  cpu 2855 ms = 99%,  os 16 ms
MFA_NTT<2^6,2052,P=0xFFF00001>*16384: 2404 ms = 853 MiB/s,  cpu 2356 ms = 98%,  os 47 ms
MFA_NTT<2^5,2052,P=0xFFF00001>*32768: 1802 ms = 1139 MiB/s,  cpu 1810 ms = 100%,  os 0 ms
MFA_NTT<2^4,2052,P=0xFFF00001>*65536: 1347 ms = 1524 MiB/s,  cpu 1342 ms = 100%,  os 0 ms
MFA_NTT<2^3,2052,P=0xFFF00001>*131072: 923 ms = 2224 MiB/s,  cpu 905 ms = 98%,  os 0 ms
MFA_NTT<2^2,2052,P=0xFFF00001>*262144: 553 ms = 3712 MiB/s,  cpu 546 ms = 99%,  os 0 ms
MFA_NTT<2^1,2052,P=0xFFF00001>*524288: 282 ms = 7289 MiB/s,  cpu 281 ms = 100%,  os 0 ms

------------------------------------------------------------------------------------------------

C:\>for %s in (16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt64g-avx2 =s %s
MFA_NTT<2^16,2052,P=65537>*16: 2712 ms = 378 MiB/s,  cpu 2714 ms = 100%,  os 0 ms
MFA_NTT<2^15,2052,P=65537>*32: 2545 ms = 403 MiB/s,  cpu 2543 ms = 100%,  os 16 ms
MFA_NTT<2^14,2052,P=65537>*64: 2380 ms = 431 MiB/s,  cpu 2387 ms = 100%,  os 0 ms
MFA_NTT<2^13,2052,P=65537>*128: 2193 ms = 468 MiB/s,  cpu 2184 ms = 100%,  os 0 ms
MFA_NTT<2^12,2052,P=65537>*256: 1889 ms = 543 MiB/s,  cpu 1872 ms = 99%,  os 16 ms
MFA_NTT<2^11,2052,P=65537>*512: 1671 ms = 614 MiB/s,  cpu 1669 ms = 100%,  os 0 ms
MFA_NTT<2^10,2052,P=65537>*1024: 1503 ms = 682 MiB/s,  cpu 1498 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=65537>*2048: 1316 ms = 780 MiB/s,  cpu 1310 ms = 100%,  os 0 ms
MFA_NTT<2^8,2052,P=65537>*4096: 1174 ms = 874 MiB/s,  cpu 1170 ms = 100%,  os 0 ms
MFA_NTT<2^7,2052,P=65537>*8192: 994 ms = 1032 MiB/s,  cpu 998 ms = 100%,  os 0 ms
MFA_NTT<2^6,2052,P=65537>*16384: 850 ms = 1208 MiB/s,  cpu 827 ms = 97%,  os 31 ms
MFA_NTT<2^5,2052,P=65537>*32768: 625 ms = 1641 MiB/s,  cpu 624 ms = 100%,  os 0 ms
MFA_NTT<2^4,2052,P=65537>*65536: 498 ms = 2062 MiB/s,  cpu 499 ms = 100%,  os 0 ms
MFA_NTT<2^3,2052,P=65537>*131072: 367 ms = 2798 MiB/s,  cpu 359 ms = 98%,  os 0 ms
MFA_NTT<2^2,2052,P=65537>*262144: 222 ms = 4611 MiB/s,  cpu 234 ms = 105%,  os 0 ms
MFA_NTT<2^1,2052,P=65537>*524288: 73 ms = 14027 MiB/s,  cpu 78 ms = 107%,  os 0 ms

C:\>for %s in (16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt64g-sse2 =s %s
MFA_NTT<2^16,2052,P=65537>*16: 4939 ms = 208 MiB/s,  cpu 4930 ms = 100%,  os 0 ms
MFA_NTT<2^15,2052,P=65537>*32: 4605 ms = 223 MiB/s,  cpu 4618 ms = 100%,  os 0 ms
MFA_NTT<2^14,2052,P=65537>*64: 4230 ms = 243 MiB/s,  cpu 4228 ms = 100%,  os 0 ms
MFA_NTT<2^13,2052,P=65537>*128: 3874 ms = 265 MiB/s,  cpu 3884 ms = 100%,  os 0 ms
MFA_NTT<2^12,2052,P=65537>*256: 3367 ms = 305 MiB/s,  cpu 3323 ms = 99%,  os 31 ms
MFA_NTT<2^11,2052,P=65537>*512: 3005 ms = 341 MiB/s,  cpu 3011 ms = 100%,  os 0 ms
MFA_NTT<2^10,2052,P=65537>*1024: 2690 ms = 381 MiB/s,  cpu 2699 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=65537>*2048: 2417 ms = 425 MiB/s,  cpu 2418 ms = 100%,  os 0 ms
MFA_NTT<2^8,2052,P=65537>*4096: 2137 ms = 480 MiB/s,  cpu 2122 ms = 99%,  os 0 ms
MFA_NTT<2^7,2052,P=65537>*8192: 1784 ms = 575 MiB/s,  cpu 1763 ms = 99%,  os 16 ms
MFA_NTT<2^6,2052,P=65537>*16384: 1507 ms = 681 MiB/s,  cpu 1513 ms = 100%,  os 0 ms
MFA_NTT<2^5,2052,P=65537>*32768: 1178 ms = 871 MiB/s,  cpu 1186 ms = 101%,  os 0 ms
MFA_NTT<2^4,2052,P=65537>*65536: 889 ms = 1154 MiB/s,  cpu 889 ms = 100%,  os 0 ms
MFA_NTT<2^3,2052,P=65537>*131072: 606 ms = 1694 MiB/s,  cpu 593 ms = 98%,  os 0 ms
MFA_NTT<2^2,2052,P=65537>*262144: 351 ms = 2925 MiB/s,  cpu 343 ms = 98%,  os 0 ms
MFA_NTT<2^1,2052,P=65537>*524288: 126 ms = 8150 MiB/s,  cpu 125 ms = 99%,  os 0 ms



C:\>for %s in (16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt32g-avx2 =s %s
MFA_NTT<2^16,2052,P=65537>*16: 10044 ms = 102 MiB/s,  cpu 10046 ms = 100%,  os 0 ms
MFA_NTT<2^15,2052,P=65537>*32: 9347 ms = 110 MiB/s,  cpu 9344 ms = 100%,  os 0 ms
MFA_NTT<2^14,2052,P=65537>*64: 8529 ms = 120 MiB/s,  cpu 8502 ms = 100%,  os 0 ms
MFA_NTT<2^13,2052,P=65537>*128: 7715 ms = 133 MiB/s,  cpu 7706 ms = 100%,  os 0 ms
MFA_NTT<2^12,2052,P=65537>*256: 6914 ms = 148 MiB/s,  cpu 6911 ms = 100%,  os 0 ms
MFA_NTT<2^11,2052,P=65537>*512: 6109 ms = 168 MiB/s,  cpu 6100 ms = 100%,  os 16 ms
MFA_NTT<2^10,2052,P=65537>*1024: 5847 ms = 175 MiB/s,  cpu 5819 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=65537>*2048: 5018 ms = 204 MiB/s,  cpu 5008 ms = 100%,  os 0 ms
MFA_NTT<2^8,2052,P=65537>*4096: 4286 ms = 239 MiB/s,  cpu 4290 ms = 100%,  os 0 ms
MFA_NTT<2^7,2052,P=65537>*8192: 3516 ms = 292 MiB/s,  cpu 3510 ms = 100%,  os 16 ms
MFA_NTT<2^6,2052,P=65537>*16384: 2835 ms = 362 MiB/s,  cpu 2839 ms = 100%,  os 0 ms
MFA_NTT<2^5,2052,P=65537>*32768: 2455 ms = 418 MiB/s,  cpu 2449 ms = 100%,  os 0 ms
MFA_NTT<2^4,2052,P=65537>*65536: 1765 ms = 581 MiB/s,  cpu 1763 ms = 100%,  os 0 ms
MFA_NTT<2^3,2052,P=65537>*131072: 1079 ms = 951 MiB/s,  cpu 1076 ms = 100%,  os 0 ms
MFA_NTT<2^2,2052,P=65537>*262144: 501 ms = 2046 MiB/s,  cpu 484 ms = 96%,  os 0 ms
MFA_NTT<2^1,2052,P=65537>*524288: 87 ms = 11755 MiB/s,  cpu 78 ms = 89%,  os 0 ms

C:\>for %s in (16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt32g-sse2 =s %s
MFA_NTT<2^16,2052,P=65537>*16: 10073 ms = 102 MiB/s,  cpu 10046 ms = 100%,  os 16 ms
MFA_NTT<2^15,2052,P=65537>*32: 9287 ms = 110 MiB/s,  cpu 9282 ms = 100%,  os 0 ms
MFA_NTT<2^14,2052,P=65537>*64: 8489 ms = 121 MiB/s,  cpu 8471 ms = 100%,  os 16 ms
MFA_NTT<2^13,2052,P=65537>*128: 7754 ms = 132 MiB/s,  cpu 7753 ms = 100%,  os 0 ms
MFA_NTT<2^12,2052,P=65537>*256: 6971 ms = 147 MiB/s,  cpu 6942 ms = 100%,  os 0 ms
MFA_NTT<2^11,2052,P=65537>*512: 6134 ms = 167 MiB/s,  cpu 6100 ms = 99%,  os 31 ms
MFA_NTT<2^10,2052,P=65537>*1024: 5671 ms = 181 MiB/s,  cpu 5678 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=65537>*2048: 4992 ms = 206 MiB/s,  cpu 4961 ms = 99%,  os 0 ms
MFA_NTT<2^8,2052,P=65537>*4096: 4293 ms = 239 MiB/s,  cpu 4274 ms = 100%,  os 0 ms
MFA_NTT<2^7,2052,P=65537>*8192: 3618 ms = 284 MiB/s,  cpu 3588 ms = 99%,  os 31 ms
MFA_NTT<2^6,2052,P=65537>*16384: 2948 ms = 348 MiB/s,  cpu 2902 ms = 98%,  os 47 ms
MFA_NTT<2^5,2052,P=65537>*32768: 2495 ms = 411 MiB/s,  cpu 2480 ms = 99%,  os 0 ms
MFA_NTT<2^4,2052,P=65537>*65536: 1789 ms = 573 MiB/s,  cpu 1794 ms = 100%,  os 0 ms
MFA_NTT<2^3,2052,P=65537>*131072: 1139 ms = 901 MiB/s,  cpu 1123 ms = 99%,  os 0 ms
MFA_NTT<2^2,2052,P=65537>*262144: 560 ms = 1831 MiB/s,  cpu 577 ms = 103%,  os 0 ms
MFA_NTT<2^1,2052,P=65537>*524288: 123 ms = 8317 MiB/s,  cpu 125 ms = 101%,  os 0 ms

------------------------------------------------------------------------------------------------

C:\>for %s in (16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt64g-avx2 -s %s
MFA_NTT<2^16,2052,P=2^32-1>*16: 7815 ms = 263 MiB/s,  cpu 7784 ms = 100%,  os 0 ms
MFA_NTT<2^15,2052,P=2^32-1>*32: 7264 ms = 283 MiB/s,  cpu 7254 ms = 100%,  os 0 ms
MFA_NTT<2^14,2052,P=2^32-1>*64: 6724 ms = 305 MiB/s,  cpu 6708 ms = 100%,  os 0 ms
MFA_NTT<2^13,2052,P=2^32-1>*128: 6201 ms = 331 MiB/s,  cpu 6209 ms = 100%,  os 0 ms
MFA_NTT<2^12,2052,P=2^32-1>*256: 5623 ms = 365 MiB/s,  cpu 5616 ms = 100%,  os 0 ms
MFA_NTT<2^11,2052,P=2^32-1>*512: 5082 ms = 404 MiB/s,  cpu 5039 ms = 99%,  os 31 ms
MFA_NTT<2^10,2052,P=2^32-1>*1024: 4629 ms = 443 MiB/s,  cpu 4618 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=2^32-1>*2048: 4109 ms = 499 MiB/s,  cpu 4103 ms = 100%,  os 0 ms
MFA_NTT<2^8,2052,P=2^32-1>*4096: 3600 ms = 570 MiB/s,  cpu 3572 ms = 99%,  os 0 ms
MFA_NTT<2^7,2052,P=2^32-1>*8192: 3109 ms = 660 MiB/s,  cpu 3042 ms = 98%,  os 62 ms
MFA_NTT<2^6,2052,P=2^32-1>*16384: 2618 ms = 784 MiB/s,  cpu 2605 ms = 100%,  os 0 ms
MFA_NTT<2^5,2052,P=2^32-1>*32768: 2176 ms = 943 MiB/s,  cpu 2168 ms = 100%,  os 0 ms
MFA_NTT<2^4,2052,P=2^32-1>*65536: 1644 ms = 1248 MiB/s,  cpu 1654 ms = 101%,  os 0 ms
MFA_NTT<2^3,2052,P=2^32-1>*131072: 1156 ms = 1775 MiB/s,  cpu 1154 ms = 100%,  os 0 ms
MFA_NTT<2^2,2052,P=2^32-1>*262144: 693 ms = 2960 MiB/s,  cpu 686 ms = 99%,  os 0 ms
MFA_NTT<2^1,2052,P=2^32-1>*524288: 303 ms = 6778 MiB/s,  cpu 296 ms = 98%,  os 0 ms

C:\>for %s in (16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt64g-sse2 -s %s
MFA_NTT<2^16,2052,P=2^32-1>*16: 5627 ms = 365 MiB/s,  cpu 5600 ms = 100%,  os 0 ms
MFA_NTT<2^15,2052,P=2^32-1>*32: 5171 ms = 397 MiB/s,  cpu 5148 ms = 100%,  os 0 ms
MFA_NTT<2^14,2052,P=2^32-1>*64: 4845 ms = 424 MiB/s,  cpu 4820 ms = 99%,  os 16 ms
MFA_NTT<2^13,2052,P=2^32-1>*128: 4449 ms = 461 MiB/s,  cpu 4399 ms = 99%,  os 47 ms
MFA_NTT<2^12,2052,P=2^32-1>*256: 4021 ms = 510 MiB/s,  cpu 4025 ms = 100%,  os 0 ms
MFA_NTT<2^11,2052,P=2^32-1>*512: 3612 ms = 568 MiB/s,  cpu 3557 ms = 98%,  os 47 ms
MFA_NTT<2^10,2052,P=2^32-1>*1024: 3241 ms = 633 MiB/s,  cpu 3245 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=2^32-1>*2048: 2883 ms = 712 MiB/s,  cpu 2870 ms = 100%,  os 0 ms
MFA_NTT<2^8,2052,P=2^32-1>*4096: 2525 ms = 813 MiB/s,  cpu 2527 ms = 100%,  os 0 ms
MFA_NTT<2^7,2052,P=2^32-1>*8192: 2201 ms = 932 MiB/s,  cpu 2184 ms = 99%,  os 16 ms
MFA_NTT<2^6,2052,P=2^32-1>*16384: 1889 ms = 1086 MiB/s,  cpu 1888 ms = 100%,  os 0 ms
MFA_NTT<2^5,2052,P=2^32-1>*32768: 1494 ms = 1373 MiB/s,  cpu 1498 ms = 100%,  os 0 ms
MFA_NTT<2^4,2052,P=2^32-1>*65536: 1127 ms = 1820 MiB/s,  cpu 1123 ms = 100%,  os 0 ms
MFA_NTT<2^3,2052,P=2^32-1>*131072: 801 ms = 2561 MiB/s,  cpu 796 ms = 99%,  os 0 ms
MFA_NTT<2^2,2052,P=2^32-1>*262144: 489 ms = 4197 MiB/s,  cpu 484 ms = 99%,  os 0 ms
MFA_NTT<2^1,2052,P=2^32-1>*524288: 217 ms = 9445 MiB/s,  cpu 218 ms = 101%,  os 0 ms



C:\>for %s in (16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt32g-avx2 -s %s
MFA_NTT<2^16,2052,P=2^32-1>*16: 8344 ms = 246 MiB/s,  cpu 8065 ms = 97%,  os 0 ms
MFA_NTT<2^15,2052,P=2^32-1>*32: 7532 ms = 272 MiB/s,  cpu 7519 ms = 100%,  os 0 ms
MFA_NTT<2^14,2052,P=2^32-1>*64: 6928 ms = 296 MiB/s,  cpu 6911 ms = 100%,  os 0 ms
MFA_NTT<2^13,2052,P=2^32-1>*128: 6393 ms = 321 MiB/s,  cpu 6380 ms = 100%,  os 16 ms
MFA_NTT<2^12,2052,P=2^32-1>*256: 5809 ms = 353 MiB/s,  cpu 5803 ms = 100%,  os 16 ms
MFA_NTT<2^11,2052,P=2^32-1>*512: 5226 ms = 393 MiB/s,  cpu 5210 ms = 100%,  os 16 ms
MFA_NTT<2^10,2052,P=2^32-1>*1024: 4748 ms = 432 MiB/s,  cpu 4742 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=2^32-1>*2048: 4253 ms = 482 MiB/s,  cpu 4243 ms = 100%,  os 0 ms
MFA_NTT<2^8,2052,P=2^32-1>*4096: 3718 ms = 552 MiB/s,  cpu 3713 ms = 100%,  os 0 ms
MFA_NTT<2^7,2052,P=2^32-1>*8192: 3255 ms = 630 MiB/s,  cpu 3245 ms = 100%,  os 0 ms
MFA_NTT<2^6,2052,P=2^32-1>*16384: 2733 ms = 751 MiB/s,  cpu 2699 ms = 99%,  os 31 ms
MFA_NTT<2^5,2052,P=2^32-1>*32768: 2234 ms = 919 MiB/s,  cpu 2231 ms = 100%,  os 0 ms
MFA_NTT<2^4,2052,P=2^32-1>*65536: 1704 ms = 1205 MiB/s,  cpu 1700 ms = 100%,  os 0 ms
MFA_NTT<2^3,2052,P=2^32-1>*131072: 1190 ms = 1724 MiB/s,  cpu 1186 ms = 100%,  os 0 ms
MFA_NTT<2^2,2052,P=2^32-1>*262144: 722 ms = 2841 MiB/s,  cpu 718 ms = 99%,  os 0 ms
MFA_NTT<2^1,2052,P=2^32-1>*524288: 321 ms = 6399 MiB/s,  cpu 312 ms = 97%,  os 0 ms

C:\>for %s in (16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt32g-sse2 -s %s
MFA_NTT<2^16,2052,P=2^32-1>*16: 5620 ms = 365 MiB/s,  cpu 5616 ms = 100%,  os 0 ms
MFA_NTT<2^15,2052,P=2^32-1>*32: 5242 ms = 391 MiB/s,  cpu 5257 ms = 100%,  os 0 ms
MFA_NTT<2^14,2052,P=2^32-1>*64: 4856 ms = 423 MiB/s,  cpu 4836 ms = 100%,  os 16 ms
MFA_NTT<2^13,2052,P=2^32-1>*128: 4516 ms = 454 MiB/s,  cpu 4493 ms = 99%,  os 16 ms
MFA_NTT<2^12,2052,P=2^32-1>*256: 4102 ms = 500 MiB/s,  cpu 4087 ms = 100%,  os 16 ms
MFA_NTT<2^11,2052,P=2^32-1>*512: 3694 ms = 556 MiB/s,  cpu 3682 ms = 100%,  os 16 ms
MFA_NTT<2^10,2052,P=2^32-1>*1024: 3271 ms = 627 MiB/s,  cpu 3276 ms = 100%,  os 0 ms
MFA_NTT<2^9,2052,P=2^32-1>*2048: 2933 ms = 700 MiB/s,  cpu 2917 ms = 99%,  os 16 ms
MFA_NTT<2^8,2052,P=2^32-1>*4096: 2585 ms = 794 MiB/s,  cpu 2590 ms = 100%,  os 0 ms
MFA_NTT<2^7,2052,P=2^32-1>*8192: 2255 ms = 910 MiB/s,  cpu 2246 ms = 100%,  os 16 ms
MFA_NTT<2^6,2052,P=2^32-1>*16384: 1937 ms = 1059 MiB/s,  cpu 1934 ms = 100%,  os 0 ms
MFA_NTT<2^5,2052,P=2^32-1>*32768: 1500 ms = 1368 MiB/s,  cpu 1498 ms = 100%,  os 0 ms
MFA_NTT<2^4,2052,P=2^32-1>*65536: 1153 ms = 1780 MiB/s,  cpu 1154 ms = 100%,  os 0 ms
MFA_NTT<2^3,2052,P=2^32-1>*131072: 815 ms = 2518 MiB/s,  cpu 827 ms = 101%,  os 0 ms
MFA_NTT<2^2,2052,P=2^32-1>*262144: 505 ms = 4060 MiB/s,  cpu 515 ms = 102%,  os 0 ms
MFA_NTT<2^1,2052,P=2^32-1>*524288: 232 ms = 8849 MiB/s,  cpu 234 ms = 101%,  os 0 ms

------------------------------------------------------------------------------------------------

C:\>for %s in (16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt64g-avx2 +s %s
MFA_NTT<2^16,2048,P=2^64-1>*16: 3675 ms = 557 MiB/s,  cpu 3666 ms = 100%,  os 0 ms
MFA_NTT<2^15,2048,P=2^64-1>*32: 3431 ms = 597 MiB/s,  cpu 3432 ms = 100%,  os 0 ms
MFA_NTT<2^14,2048,P=2^64-1>*64: 3093 ms = 662 MiB/s,  cpu 3089 ms = 100%,  os 0 ms
MFA_NTT<2^13,2048,P=2^64-1>*128: 2796 ms = 733 MiB/s,  cpu 2808 ms = 100%,  os 0 ms
MFA_NTT<2^12,2048,P=2^64-1>*256: 2417 ms = 847 MiB/s,  cpu 2402 ms = 99%,  os 16 ms
MFA_NTT<2^11,2048,P=2^64-1>*512: 2073 ms = 988 MiB/s,  cpu 2075 ms = 100%,  os 0 ms
MFA_NTT<2^10,2048,P=2^64-1>*1024: 1946 ms = 1052 MiB/s,  cpu 1934 ms = 99%,  os 16 ms
MFA_NTT<2^9,2048,P=2^64-1>*2048: 1691 ms = 1211 MiB/s,  cpu 1669 ms = 99%,  os 16 ms
MFA_NTT<2^8,2048,P=2^64-1>*4096: 1469 ms = 1394 MiB/s,  cpu 1451 ms = 99%,  os 0 ms
MFA_NTT<2^7,2048,P=2^64-1>*8192: 1182 ms = 1733 MiB/s,  cpu 1186 ms = 100%,  os 0 ms
MFA_NTT<2^6,2048,P=2^64-1>*16384: 954 ms = 2146 MiB/s,  cpu 952 ms = 100%,  os 0 ms
MFA_NTT<2^5,2048,P=2^64-1>*32768: 796 ms = 2572 MiB/s,  cpu 811 ms = 102%,  os 0 ms
MFA_NTT<2^4,2048,P=2^64-1>*65536: 553 ms = 3704 MiB/s,  cpu 546 ms = 99%,  os 0 ms
MFA_NTT<2^3,2048,P=2^64-1>*131072: 361 ms = 5666 MiB/s,  cpu 359 ms = 99%,  os 0 ms
MFA_NTT<2^2,2048,P=2^64-1>*262144: 199 ms = 10308 MiB/s,  cpu 203 ms = 102%,  os 0 ms
MFA_NTT<2^1,2048,P=2^64-1>*524288: 93 ms = 22026 MiB/s,  cpu 94 ms = 101%,  os 0 ms

C:\>for %s in (16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1) do @ntt64g-sse2 +s %s
MFA_NTT<2^16,2048,P=2^64-1>*16: 4123 ms = 497 MiB/s,  cpu 4118 ms = 100%,  os 0 ms
MFA_NTT<2^15,2048,P=2^64-1>*32: 3785 ms = 541 MiB/s,  cpu 3791 ms = 100%,  os 0 ms
MFA_NTT<2^14,2048,P=2^64-1>*64: 3506 ms = 584 MiB/s,  cpu 3510 ms = 100%,  os 0 ms
MFA_NTT<2^13,2048,P=2^64-1>*128: 3226 ms = 635 MiB/s,  cpu 3198 ms = 99%,  os 0 ms
MFA_NTT<2^12,2048,P=2^64-1>*256: 2837 ms = 722 MiB/s,  cpu 2824 ms = 100%,  os 16 ms
MFA_NTT<2^11,2048,P=2^64-1>*512: 2495 ms = 821 MiB/s,  cpu 2480 ms = 99%,  os 16 ms
MFA_NTT<2^10,2048,P=2^64-1>*1024: 2252 ms = 909 MiB/s,  cpu 2246 ms = 100%,  os 0 ms
MFA_NTT<2^9,2048,P=2^64-1>*2048: 1980 ms = 1035 MiB/s,  cpu 1950 ms = 99%,  os 16 ms
MFA_NTT<2^8,2048,P=2^64-1>*4096: 1781 ms = 1150 MiB/s,  cpu 1778 ms = 100%,  os 0 ms
MFA_NTT<2^7,2048,P=2^64-1>*8192: 1497 ms = 1368 MiB/s,  cpu 1482 ms = 99%,  os 0 ms
MFA_NTT<2^6,2048,P=2^64-1>*16384: 1262 ms = 1622 MiB/s,  cpu 1217 ms = 96%,  os 47 ms
MFA_NTT<2^5,2048,P=2^64-1>*32768: 951 ms = 2153 MiB/s,  cpu 952 ms = 100%,  os 0 ms
MFA_NTT<2^4,2048,P=2^64-1>*65536: 727 ms = 2818 MiB/s,  cpu 718 ms = 99%,  os 0 ms
MFA_NTT<2^3,2048,P=2^64-1>*131072: 524 ms = 3907 MiB/s,  cpu 515 ms = 98%,  os 0 ms
MFA_NTT<2^2,2048,P=2^64-1>*262144: 346 ms = 5919 MiB/s,  cpu 343 ms = 99%,  os 0 ms
MFA_NTT<2^1,2048,P=2^64-1>*524288: 194 ms = 10582 MiB/s,  cpu 187 ms = 97%,  os 0 ms

```

### Slow NTT

Comparison of slow O(N^2) NTT with fast algorithms:

```
C:\>timer ntt64g.exe q 20 32
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
