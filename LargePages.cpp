
#include <stdio.h>
#include <windows.h>

typedef unsigned int   uint;
typedef unsigned long long qword;

const uint g_PageFlag0 = 0x1000;
const uint g_PageMask0 = 0x1000-1;
uint g_PageFlag = g_PageFlag0;
uint g_PageMask = g_PageMask0;


#pragma comment(lib, "advapi32.lib")

void InitLargePages( void ) 
{
    HANDLE hToken = 0;
    LUID luid;
    TOKEN_PRIVILEGES tp;
    OpenProcessToken( GetCurrentProcess(), TOKEN_ADJUST_PRIVILEGES, &hToken );
    LookupPrivilegeValue( NULL, "SeLockMemoryPrivilege", &luid );
    tp.PrivilegeCount = 1;
    tp.Privileges[0].Luid = luid;
    tp.Privileges[0].Attributes = SE_PRIVILEGE_ENABLED;
    AdjustTokenPrivileges( hToken, FALSE, &tp, sizeof(tp), 0, 0 ); 
    CloseHandle(hToken);

    uint size = 0x1000;
    typedef uint (__stdcall *GetLPMP)( void );
    GetLPMP LPM = (GetLPMP)GetProcAddress(GetModuleHandle(TEXT("kernel32.dll")),"GetLargePageMinimum");
    if( LPM ) size = LPM();
    if( (size==0) || ((size&(size-1))!=0) ) size=0x1000;
    g_PageMask = size-1;
    g_PageFlag = g_PageFlag0;
    if( size>0x1000 ) g_PageFlag |= 0x20000000;
}


template< class T >
T* VAlloc( qword size ) 
{
    void* r;
    size *= sizeof(T);
    qword s = (size+g_PageMask) & (~qword(g_PageMask));
    Retry:
#ifndef MY_CPU_64BIT
    if( s>=(qword(1)<<32) ) r=0; else
#endif
    r = VirtualAlloc(0, s, g_PageFlag, PAGE_READWRITE );
    if( (r==0) && (g_PageMask!=g_PageMask0) ) {
        g_PageFlag = g_PageFlag0;
        g_PageMask = g_PageMask0; 
        s = size;
        goto Retry;
    }
    return (T*)r;
}


template< class T >
void VFree( T* p ) 
{
    VirtualFree( p, 0, MEM_RELEASE );
}



int main( void ) {

  InitLargePages();

  uint* p = VAlloc<uint>( 1<<20 );

  printf( "p=%08X Flags=%08X\n", p, g_PageFlag );

}
