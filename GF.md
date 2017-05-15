
### Program usage

`prime N` - check whether N is prime, print divisors of N, and search, starting at N+1, for prime numbers as well as numbers with only small divisors.

This program is useful for researching ring properties, in particular maximal order. Some other ring/field properties are checked by the [NTT program](NTT.md).


### Lucky number: choosing the best base for computations

Since GF(2^n) doesn't have much roots of unity, efficient NTT-based Reed-Solomon codes implementation can't perform computations in this field.
Instead, we need to use other Galois Field, or even Ring modulo some number. GF(p^n) has a maximal order of p^n-1.
For rings, the maximal order is defined by complex formula that you can find in chapter `39.7 Composite modulus: the ring Z=mZ` of [FxtBook](http://www.jjj.de/fxt/fxtbook.pdf).

Good candidates for the base have the form p=2^a+b where b is a small positive or negative number.
Such bases allow efficient data storage and fast radix conversion from/to 2^a base (see [below](#data-packing)),
while computations become more efficient when b=1 and especially when b=-1.

Good candidates for the base are:
- GF(0xFFF00001) - my current favourite. `0xFFF00000 = 2^20*3*3*5*7*13` has 504 divisors overall, and for random 32-bit N we can find a divisor that is only a few percents larger.
This means that when we need to process N source blocks, we can perform NTT using only a few percents more memory than the source data occupy. Source blocks up to 4 KB
can be converted into 1024 numbers in the 0..0xFFF00000 range plus a single bit.
- GF(0x10001) - computations may be a bit faster, but NTT order may be only 2,4..65536.
Compact memory storage require to recode data into base-0x10000 plus one overflow bit per 32K values, that may slowdown the NTT operation.
- GF(0x10001^2) - may be slightly faster than GF(0xFFF00001), maximal NTT order is `0x10001^2-1 = 2^17*3*3*11*331`, the same storage problems.
- Mod(2^32-1) - 2x faster, but NTT order may be only 2,4..65536. May be used as fast algorithm for block counts equal to 2^N or slightly lower, for N<=16.
- GF(2^31-1) - also 2x faster, max. order is large, but its divisors `p-1 = 2*3*3*7*11*31*151*331` doesn't look fascinating.
- GF(p^2) for p=2^31-1 - again 2x faster, max order `p^2-1 = 2^32*3*3*7*11*31*151*331` so the divisors are almost as dense as for GF(0xFFF00001).
It may be the best base, but its efficient implementation will require extra work.
- GF(2^61-1) - fastest for pure (non-SIMD) x64 code, but `p-1 = 2*3*3*5*5*7*11*13*31*41*61*151*331*1321` has not too much divisors
- GF(p^2) for p=2^61-1 may be also interesting since it's almost as fast as GF(2^61-1) and `p^2-1 = 2^62*3*3*5*5*7*11*13*31*41*61*151*331*1321`,
providing ideal coverage of integer space by divisors. I think that it may be 3-4x faster than GF(0xFFF00001).
It may be the best base for x64, but its efficient implementation will require extra work.
- Mod(2^64-1) - among fastest variants for x64, but NTT order should be a divisor of `2^16*3*5*17449`, so it doesn't provide too much choice.

Intermediate data can be stored unnormalized, i.e. as arbitrary 32/64-bit value.
Normalization required only when operation result may overflow its register size, and it can be partial - only packing the result back to the register size.
Full normalization is required only on the final data.
For example, 32-bit code for operations Mod(2^32-1) is as simple as:

```
EAX += EBX
    add eax, ebx
    adc eax, 0

EAX -= EBX
    not ebx
    add eax, ebx
    adc eax, 0

EAX *= EBX
    mul ebx
    add eax,edx
    adc eax, 0
```

64-bit code for operations Mod(2^64-1) is essentially the same.
The entire Butterfly computation consists of these three operations, so (on Skylake) it can be performed in 3 CPU cycles (limited by ADC throughput).
One Butterfly operation processes 8 bytes for Mod(2^32-1) or 16 bytes for Mod(2^64-1), so it will process 10/20 GB/s per core.
NTT(2^N) require N passes over data, so its speed will be 10/N or 20/N GB/s per core.
F.e. NTT(2^20) using Mod(2^64-1) operations will run at 1 GB/s per core, 4 GB/s overall!!!

Unfortunately, it seems that while we can perform NTT/iNTT and thus RS encoding in arbitrary rings at O(N*log(N)) speed,
RS decoding require full division support and therefore can be implemented only in Galois Fields.

Great overview of GF(p) Butterfly optimizations was provided by David Harvey in the talk
[Faster arithmetic for number-theoretic transforms](http://web.maths.unsw.edu.au/~davidharvey/talks/fastntt-2-talk.pdf).


<a name="data-packing"/>

### Efficient data packing

I also found a way to recode 4 KBytes to the 0xFFF00001 base using just 1 extra bit.
Idea is the following: we have 1024 32-bit values and we don't use the highest value 0xFFF00000 in our encoding,
so we just need to recode 1024 12-bit values from base 4096 to base 4095.

Decoding algorithm (translating data from the base 4095 to the base 4096): if the extra bit is 0, then keep input data intact.
If this bit is 1, it means that there is at least one 0xFFF value in output data.
In that case, first entry of the input data holds index of the first 0xFFF in output data (10 bits), plus the flag (1 - there are more 0xFFF in output data).
If the flag is 1, then the next input item, again, contains index of the next 0xFFF in output data plus continuation flag.
After the flag 0, remaining input items contains values of remaining output elements.

---

Once input (source) data are recoded in this way, we need to store the extra bit in the way which ensure that the bit can be restored
in any situation when the data block can be restored. The best way to ensure this, that I found, is to save the extra bit as one more (1025'th)
source word. So, all operations are performed on 4100-byte blocks, and parity sectors stored are 4100-byte long. Sad, but I don't see better choice.
Remaining bits of the extra word can be used to store block checksum, although i don't see much gain in that.

Of course, when 64-bit base and/or GF(p^2) field are used, extra data will be increased to 8-16 bytes.
