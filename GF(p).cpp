/// Implementation of the GF(P) Galois field
#if defined(_M_X64) || defined(_M_AMD64) || defined(__x86_64__)
#define MY_CPU_AMD64
#endif

#if defined(MY_CPU_AMD64) || defined(_M_IA64)
#define MY_CPU_64BIT
#endif

#if _MSC_VER
#define MSVC_ONLY
#endif

#if __GNUC__ && defined(MY_CPU_64BIT)
#include <inttypes.h>
typedef unsigned __int128 uint128_t;
#elif _MSC_VER && defined(MY_CPU_64BIT)
#include <intrin.h>
#define constexpr  /* incompatible with __umulh :( */
#endif

#define unlikely /* to do... */


/***********************************************************************************************************************
*** Main operations in GF(P) *******************************************************************************************
************************************************************************************************************************/

// Twice wider type to hold intermediate MUL results
template<typename Type> struct Double           {};
template<>              struct Double<uint32_t> {typedef uint64_t T;};


template <typename T, T P>
constexpr T GF_Sub (T X, T Y)
{
    T res = X - Y;
    return res + (res>X)*P;   // res>X? res+P : res
}

template <typename T, T P>
constexpr T GF_Add (T X, T Y)
{
    return GF_Sub<T,P> (X, P-Y);
}


#if __GNUC__ && defined(MY_CPU_64BIT)
// Alternative GF_Mul64 implementation for 64-bit GCC
template<>              struct Double<uint64_t>    {typedef uint128_t T;};

// 4x wider type to hold intermediate MUL results
template<typename Type> struct Quadruple           {};
template<>              struct Quadruple<uint32_t> {typedef uint128_t T;};

// The binary logarithm, rounded down
template <typename T>
static constexpr int trunc_log2 (T x)
{return x<=1? 0 : 1+trunc_log2(x/2);}

template <typename T, T P>
constexpr T GF_Mul64 (T X, T Y)
{
    using DoubleT = typename Double<T>::T;
    using QuadT   = typename Quadruple<T>::T;

    // See chapter "16.9 Division" in the http://www.agner.org/optimize/optimizing_assembly.pdf
    constexpr int     BITS  = trunc_log2(P) + 8*sizeof(DoubleT);
    constexpr QuadT   invP2 = (QuadT(2) << BITS) / P;  // double the invP value
    constexpr DoubleT invP  = (invP2+1) / 2;           // rounded invP value
    constexpr DoubleT extra = 1 - (invP2 & 1);         // 1 if invP was rounded down, 0 otherwise

    DoubleT res = DoubleT(X)*Y;
    res -= DoubleT(((res+extra)*QuadT(invP)) >> BITS) * P;
    return T(res);
}

#elif MSVC_ONLY && defined(MY_CPU_64BIT)
// Alternative GF_Mul64 implementation made with MSVC intrinsics

template <typename T, T P>
constexpr T GF_Mul64 (T X, T Y)
{
    using DoubleT = typename Double<T>::T;
    constexpr DoubleT estInvP = ((DoubleT(1)<<63) / P) << 1;
    constexpr DoubleT invP    = (estInvP*P > (estInvP+1)*P? estInvP : estInvP+1);

    DoubleT res = DoubleT(X)*Y;
    res -= __umulh(res,invP) * P;
    return T(res>=P? res-P : res);
}

#else

// GF_Mul64 is optimized for 64-bit CPUs
template <typename T, T P>
constexpr T GF_Mul64 (T X, T Y)
{
    using DoubleT = typename Double<T>::T;
    return T((DoubleT(X)*Y) % P);
}

#endif

// GF_Mul32 is optimized for 32-bit CPUs, SIMD and GPUs
template <typename T, T P>
constexpr T GF_Mul32 (T X, T Y)
{
    using DoubleT = typename Double<T>::T;
    // invP32 := (2**64)/P - 2**32  :  if 2**31<P<2**32, then 2**32 < (2**64)/P < 2**33, and invP32 is a 32-bit value
    const DoubleT estInvP = ((DoubleT(1)<<63) / P) << 1;                        // == invP & (~1)
    const T       invP32  = T(estInvP*P > (estInvP+1)*P? estInvP : estInvP+1);  // we can't use 1<<64 for exact invP computation so we add the posible 1 in other way

    DoubleT res = DoubleT(X)*Y;
    res  -=  ((res + (res>>32)*invP32) >> 32) * P;    // The same as res -= ((res*invP) >> 64) * P, where invP = (2**64)/P, but optimized for 32-bit computations
    return T(res>=P? res-P : res);
}

#if defined(MY_CPU_64BIT) && !defined(USE_GF_Mul32)
#define GF_Mul GF_Mul64
#else
#define GF_Mul GF_Mul32
#endif


/***********************************************************************************************************************
*** Optimized operations for various P *********************************************************************************
************************************************************************************************************************/

#if 0 && MSVC_ONLY && defined(MY_CPU_64BIT)
// Optimized operations for P=0xFFF00001 (or any other 32-bit P)
template <> constexpr uint32_t GF_Add<uint32_t,0xFFF00001> (uint32_t X, uint32_t Y)
{
    uint32_t res, temp;
    auto carry = _addcarry_u32(0,X,Y,&temp);
    return temp - (carry? 0xFFF00001 : 0);
}

template <> constexpr uint32_t GF_Sub<uint32_t,0xFFF00001> (uint32_t X, uint32_t Y)
{
    uint32_t res, temp;
    auto carry = _subborrow_u32(0,X,Y,&temp);
    return temp + (carry? 0xFFF00001 : 0);
}
#endif


// Optimized operations for P=0x10001
template <> constexpr uint32_t GF_Add<uint32_t,0x10001> (uint32_t X, uint32_t Y)
{
    uint32_t res = X+Y;
    return res - (res>=0x10001)*0x10001;   // res>P? res-P : res
}
template <> constexpr uint32_t GF_Sub<uint32_t,0x10001> (uint32_t X, uint32_t Y)
{
    uint32_t res = X-Y;
    return res + (int32_t(res)<0)*0x10001;   // res<0? res+P : res
}
template <> constexpr uint32_t GF_Mul<uint32_t,0x10001> (uint32_t X, uint32_t Y)
{
    uint32_t res = X*Y;
    if (unlikely(res==0) && X && Y)
        return 1;  // 65536*65536
    res = (res&0xFFFF) - (res>>16);
#ifdef MY_CPU_64BIT
    return res + (int32_t(res)<0)*0x10001;        // faster for gcc64
#else
    return int32_t(res)<0? res+0x10001 : res;     // faster for gcc32
#endif
}


// Optimized operations for P=0xFFFFFFFF
// Note: finally data should be normalized, i.e. 0xFFFFFFFF mapped to 0
template <> constexpr uint32_t GF_Add<uint32_t,0xFFFFFFFF> (uint32_t X, uint32_t Y)
{
    uint64_t res = uint64_t(X)+Y;
    return uint32_t(res) + uint32_t(res>>32);
}
template <> constexpr uint32_t GF_Sub<uint32_t,0xFFFFFFFF> (uint32_t X, uint32_t Y)
{
    return GF_Add<uint32_t,0xFFFFFFFF> (X,~Y);
}
template <> constexpr uint32_t GF_Mul<uint32_t,0xFFFFFFFF> (uint32_t X, uint32_t Y)
{
    uint64_t res = uint64_t(X) * Y;
    res = (res&0xFFFFFFFF) + (res>>32);
    return uint32_t(res) + uint32_t(res>>32);
}


// Optimized operations for P=0xFFFFFFFFFFFFFFFF
// Note: finally data should be normalized, i.e. 0xFFFFFFFFFFFFFFFF mapped to 0
#if __GNUC__ && defined(MY_CPU_64BIT)

template <> constexpr uint64_t GF_Add<uint64_t,0xFFFFFFFFFFFFFFFF> (uint64_t X, uint64_t Y)
{
    uint64_t res = X+Y;
    return res + (res<X);
}
template <> constexpr uint64_t GF_Sub<uint64_t,0xFFFFFFFFFFFFFFFF> (uint64_t X, uint64_t Y)
{
    uint64_t res = X-Y;
    return res - (res>X);
}
template <> constexpr uint64_t GF_Mul<uint64_t,0xFFFFFFFFFFFFFFFF> (uint64_t X, uint64_t Y)
{
    uint128_t res = uint128_t(X) * Y;
    uint64_t a = res,  b = res>>64;
    uint64_t c = a+b;
    return c + (c<a);
}

#elif MSVC_ONLY && defined(MY_CPU_64BIT)

template <> constexpr uint64_t GF_Add<uint64_t,0xFFFFFFFFFFFFFFFF> (uint64_t X, uint64_t Y)
{
    uint64_t res, temp;
    auto carry = _addcarry_u64(0,X,Y,&temp);
    _addcarry_u64(carry,temp,0,&res);
    return res;
}
template <> constexpr uint64_t GF_Sub<uint64_t,0xFFFFFFFFFFFFFFFF> (uint64_t X, uint64_t Y)
{
    uint64_t res, temp;
    auto carry = _subborrow_u64(0,X,Y,&temp);
    _subborrow_u64(carry,temp,0,&res);
    return res;
}
template <> constexpr uint64_t GF_Mul<uint64_t,0xFFFFFFFFFFFFFFFF> (uint64_t X, uint64_t Y)
{
    uint64_t a, b = _umul128(X,Y,&a), res, temp;
    auto carry = _addcarry_u64(0,a,b,&temp);
    _addcarry_u64(carry,temp,0,&res);
    return res;
}
#endif


/***********************************************************************************************************************
*** Extra operations in GF(P) ******************************************************************************************
************************************************************************************************************************/

template <typename T, T P>
constexpr T GF_Pow (T X, T N)
{
    T res = 1;
    for ( ; N; N/=2)
    {
        if (N&1)  res = GF_Mul<T,P> (res,X);
        X = GF_Mul<T,P> (X,X);
    }
    return res;
}


// Some primary root of 1 of power N
template <typename T, T P>
constexpr T GF_Root (T N)
{
    static_assert (P==0xFFF00001 || P==0x10001, "Only GF(0xFFF00001) and GF(0x10001) are supported by generic GF_Root");
    T main_root  =  (P==0x10001? 3 : 19);  // root of power P-1 in the GF(P)

    //assert ((P-1) % N  ==  0);
    return GF_Pow<T,P> (main_root, (P-1) / N);
}
template <> constexpr uint32_t GF_Root<uint32_t,0xFFFFFFFF> (uint32_t N)
{
    uint32_t main_root = 7;  // root of power 65536 in Z/mZ(0xFFFFFFFF)
    //assert (65536 % N  ==  0);
    return GF_Pow<uint32_t,0xFFFFFFFF> (main_root, 65536 / N);
}
#ifdef MY_CPU_64BIT
template <> constexpr uint64_t GF_Root<uint64_t,0xFFFFFFFFFFFFFFFF> (uint64_t N)
{
    uint64_t main_root = 7;  // root of power 2^16*3*5*17449 in Z/mZ(0xFFFFFFFFFFFFFFFF)
    //assert (uint64_t(65536)*3*5*17449) % N  ==  0);
    return GF_Pow<uint64_t,0xFFFFFFFFFFFFFFFF> (main_root, (uint64_t(65536)*3*5*17449) / N);
}
#endif        


template <typename T, T P>
constexpr T GF_Inv (T X)
{
    return GF_Pow<T,P> (X, P==0xFFFFFFFFFFFFFFFF? 0xFFFFFFFF : P==0xFFFFFFFF? 0xFFFF : P-2);
}


template <typename T, T P>
constexpr T GF_Div (T X, T Y)
{
    T InvY = GF_Inv<T,P> (Y);
    return GF_Mul<T,P> (X, InvY);
}


// Normalize value, i.e. return X%P
// Required after optimized operations with P=2^32-1 or 2^64-1
template <typename T, T P>
constexpr T GF_Normalize (T X)
{
    return X % P;
}
