#ifndef _RAND_UTIL_H_
#define _RAND_UTIL_H_

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

uint64_t rand_uniform(uint64_t max);
double rand_normal(double mean, double sd);
double rand_zipfian(double s, double max, uint64_t source);

#endif
