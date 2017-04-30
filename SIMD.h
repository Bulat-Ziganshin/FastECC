/// Constants for specification of enabled SIMD extension in the compiler cmdline
/// Example: gcc -mavx2 -DSIMD=AVX2  / cl -arch:AVX2 -DSIMD=AVX2
///          #if SIMD >= AVX2

#define SSE     1
#define SSE2    2
#define SSE3    3
#define SSSE3   4
#define SSE41   5
#define SSE42   6
#define AVX     7
#define AVX2    8
