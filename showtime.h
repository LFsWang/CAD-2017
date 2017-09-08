#ifndef SHOW_TIME_H
#define SHOW_TIME_H
#include<ctime>
#include<omp.h>
#include<cstdio>

inline void _showclock(const char *str=nullptr)
{
#ifdef _WIN32
    long long CL_PER_SEC = CLOCKS_PER_SEC;//1000;
#else
    long long CL_PER_SEC = CLOCKS_PER_SEC;//1000000;//test on centos
#endif
    static long long last = 0;
    auto show_time = [&](long long time){
        long long ms = time%CL_PER_SEC; time/=CL_PER_SEC;
        long long sec = time%60;  time/=60;
        long long min = time%60;  time/=60;
        printf("%2llu:%02llu:%02llu %06llu",time,min,sec,ms);
    };
    
    long long now = std::clock();
    if(str)printf("%s ,",str);
    printf("Time:");show_time(now);printf("\t(");
    long long diff = now - last;
    show_time(diff);printf(")\n");
    last = now;
}

inline void showclock(const char *str=nullptr)
{
	static double begin = omp_get_wtime();
    static double last = begin;
    auto show_time = [&](double time){
        printf(" %fs ",time);
    };
    
    double now = omp_get_wtime();//std::clock();
    if(str)printf("%s ,",str);
    printf("Time:");show_time(now-begin);printf("\t(");
    double diff = now - last;
    show_time(diff);printf(")\n");
    last = now;
}
#endif