@echo off
set name=ntt
set main=main.cpp
if "%1"=="rs" set name=rs
if "%1"=="rs" set main=rs.cpp
if "%1"=="rs" shift
set options_ms=-MP -Gy -GR- -nologo -Fo%TEMP%\ -Fp%TEMP%\ %main% user32.lib shell32.lib ole32.lib advapi32.lib %1 %2 %3 %4 %5 %6 %7 %8 %9 -link -LARGEADDRESSAWARE
set options_ms=-GL %options_ms%
set options_ms_cl=-O2 -EHsc %options_ms%
set options_ms_x86=-MACHINE:x86 -SUBSYSTEM:CONSOLE,5.01
set options_ms_x64=-MACHINE:x64 -SUBSYSTEM:CONSOLE,5.02

::call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86_amd64
::call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\Tools\VsDevCmd.bat" -arch=x64 -host_arch=x64 -no_logo
cl -Fe%name%64m.exe -Fa%name%64.asm %options_ms_cl% %options_ms_x64%

::call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86
::call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars32.bat"
call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\Tools\VsDevCmd.bat" -arch=x86 -host_arch=x64 -no_logo
cl -Fe%name%32m.exe -Fa%name%32.asm -arch:SSE2 %options_ms_cl% %options_ms_x86%

::g++ -std=c++14 -m64 -O3 %main% -static -fopenmp -o%name%64g-avx2 -mavx2 -DSIMD=AVX2
::g++ -std=c++14 -m64 -O3 %main% -static -fopenmp -o%name%64g-sse2        -DSIMD=SSE2
::g++ -std=c++14 -m64 -O3 %main% -static -fopenmp -o%name%64g
::g++ -std=c++14 -m32 -O3 %main% -static -fopenmp -o%name%32g-avx2 -mavx2 -DSIMD=AVX2 -Xlinker --large-address-aware
::g++ -std=c++14 -m32 -O3 %main% -static -fopenmp -o%name%32g-sse2 -msse2 -DSIMD=SSE2 -Xlinker --large-address-aware
::g++ -std=c++14 -m32 -O3 %main% -static -fopenmp -o%name%32g      -mmmx              -Xlinker --large-address-aware

::cl -Feprime.exe -O2 -EHsc prime.cpp -link %options_ms_x86%

del *.exe.bak *.obj *.res >nul 2>nul
