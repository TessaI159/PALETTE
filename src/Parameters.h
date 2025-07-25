#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdbool.h>

enum Space { OK, CIE, SRGB };
enum DiffFunc { EUCLIDEAN, D94, D2000 };

struct Parameters {
	bool	      cw;
	enum Space    space;
	enum DiffFunc diff_fn;
};

extern struct Parameters g_params;

#endif
