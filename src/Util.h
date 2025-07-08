#ifndef UTIL_H
#define UTIL_H

#include <stdint.h>
typedef struct Color Color;

uint64_t time_diff(double (*func)(Color *, Color *), Color col1, Color col2,
		   const uint64_t RUNS);

#endif
