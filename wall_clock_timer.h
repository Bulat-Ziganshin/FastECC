/**
 * Copyright 1993-2013 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#ifndef TIMER_H
#define TIMER_H

#include <stdlib.h>

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <sys/time.h>
#endif

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
double PCFreq = 0.0;
__int64 timerStart = 0;
#else
struct timeval timerStart;
#endif

void StartTimer()
{
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    LARGE_INTEGER li;

    if (!QueryPerformanceFrequency(&li))
    {
        printf("QueryPerformanceFrequency failed!\n");
    }

    PCFreq = (double)li.QuadPart/1000.0;
    QueryPerformanceCounter(&li);
    timerStart = li.QuadPart;
#else
    gettimeofday(&timerStart, NULL);
#endif
}

// Elapsed time, milliseconds
double GetTimer()
{
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    return (double)(li.QuadPart-timerStart)/PCFreq;
#else
    struct timeval timerStop, timerElapsed;
    gettimeofday(&timerStop, NULL);
    timersub(&timerStop, &timerStart, &timerElapsed);
    return timerElapsed.tv_sec*1000.0+timerElapsed.tv_usec/1000.0;
#endif
}

// Returns number of seconds spent in the current process
void GetProcessKernelUserTimes (double *KernelTime, double *UserTime)
{
    FILETIME kt, ut, x, y;
    int ok = GetProcessTimes(GetCurrentProcess(),&x,&y,&kt,&ut);
    *KernelTime = !ok? -1 : ((((unsigned long long)(kt.dwHighDateTime) << 32) + kt.dwLowDateTime)) / 1e7;
    *UserTime   = !ok? -1 : ((((unsigned long long)(ut.dwHighDateTime) << 32) + ut.dwLowDateTime)) / 1e7;
}
#endif // TIMER_H

void time_it (double size, const char* name, std::function<void()> Code)
{
    static int _ = (StartTimer(),0);
    double start = GetTimer(), KernelTime[2], UserTime[2];
    GetProcessKernelUserTimes (KernelTime, UserTime);
    Code();
    GetProcessKernelUserTimes (KernelTime+1, UserTime+1);
    double wall_time = GetTimer()-start;
    double cpu_time  = (UserTime[1] - UserTime[0]) * 1000;
    printf("%s: %.0lf ms = %.0lf MiB/s,  cpu %.0lf ms = %.0lf%%,  os %.0lf ms\n",
        name, wall_time, (size / wall_time)*1000 / (1<<20),
        cpu_time, cpu_time/wall_time*100,
        (KernelTime[1]-KernelTime[0])*1000);
}
