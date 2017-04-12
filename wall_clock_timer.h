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

// time elapsed in ms
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
