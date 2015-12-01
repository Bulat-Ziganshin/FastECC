@echo off
set name=rs
set options=%name%.cpp
set options_ms=-MP -Gy -GR- -nologo -Fo%TEMP%\ -Fp%TEMP%\ %options% user32.lib shell32.lib ole32.lib advapi32.lib %* -link -LARGEADDRESSAWARE
set options_ms=-GL %options_ms%
set options_ms_cl=-O2 -EHsc %options_ms%
set options_ms_x86=-MACHINE:x86 -SUBSYSTEM:CONSOLE,5.01
set options_ms_x64=-MACHINE:x64 -SUBSYSTEM:CONSOLE,5.02

call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86_amd64
cl -Fe%name%64m.exe -Fa %options_ms_cl% %options_ms_x64%

call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86
cl -Fe%name%32m.exe %options_ms_cl% %options_ms_x86%

del *.exe.bak *.obj *.res >nul 2>nul
