#pragma once
#include <time.h>

double calculate_elapsed_time(struct timeval start, struct timeval stop)
{
    return (stop.tv_sec + stop.tv_usec * 1e-6) - (start.tv_sec + start.tv_usec * 1e-6);
}
