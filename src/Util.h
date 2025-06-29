#ifndef UTIL_H
#define UTIL_H

typedef struct Color Color;

unsigned long long time_diff(double (*func)(Color, Color), Color col1, Color col2,
		const long int RUNS);

#endif
