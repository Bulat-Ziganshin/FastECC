#include <windows.h>

const DWORD    g_PageFlag0 = 0x1000;
const uint64_t g_PageMask0 = 0x1000-1;
DWORD          g_PageFlag  = g_PageFlag0;
uint64_t       g_PageMask  = g_PageMask0;


#pragma comment(lib, "advapi32.lib")

void InitLargePages()
{
    HANDLE hToken = 0;
    LUID luid;
    TOKEN_PRIVILEGES tp;
    OpenProcessToken (GetCurrentProcess(), TOKEN_ADJUST_PRIVILEGES, &hToken);
    LookupPrivilegeValue (NULL, "SeLockMemoryPrivilege", &luid);
    tp.PrivilegeCount = 1;
    tp.Privileges[0].Luid = luid;
    tp.Privileges[0].Attributes = SE_PRIVILEGE_ENABLED;
    AdjustTokenPrivileges (hToken, FALSE, &tp, sizeof(tp), 0, 0);
    CloseHandle(hToken);

    uint64_t size = 0x1000;
    typedef SIZE_T (__stdcall *GetLPMP)();
    GetLPMP LPM = (GetLPMP)GetProcAddress(GetModuleHandle(TEXT("kernel32.dll")),"GetLargePageMinimum");
    if (LPM)  size = LPM();
    if (size==0 || (size&(size-1))!=0)  size=0x1000;
    g_PageMask = size-1;
    g_PageFlag = g_PageFlag0;
    if (size>0x1000)  g_PageFlag |= 0x20000000;
}


template< class T >
T* VAlloc (uint64_t size)
{
    size *= sizeof(T);
    uint64_t PageMask = g_PageMask;
    uint64_t s = (size+PageMask) & (~PageMask);
    if (s > size_t(-1)) {
        return 0;
    }

    void* r = 0;
    if (size > g_PageMask) {
        r = VirtualAlloc(0, s, g_PageFlag, PAGE_READWRITE);       // alloc using 2 MB pages only if the size is no less than one page
    }
    if (r==0 && PageMask!=g_PageMask0) {
        PageMask = g_PageMask0;
        s = (size+PageMask) & (~PageMask);
        r = VirtualAlloc(0, size, g_PageFlag0, PAGE_READWRITE);   // alloc using 4 KB pages, if preceding attempt failed
    }
    printf("Allocated %.0lf MiB with %s\n", s/1048576.0, PageMask==0x1fffff? "2MiB pages": PageMask==0xfff? "4KiB pages": "unknown pagesize");
    return (T*)r;
}


template< class T >
void VFree( T* p )
{
    VirtualFree(p, 0, MEM_RELEASE);
}
