#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdbool.h>

enum Space { OK, CIE, SRGB };
enum DiffFunc { EUCLIDEAN, D94, D2000 };

struct Parameters {
  bool cw; /*<< A global parameter indicating whether or not to use chromatic
              weighting */
  enum Space space; /*<< A global parameter indicating what space colors are
                      stored in */
  enum DiffFunc
      diff_fn; /*<< A global parameter indicating which difference
                 function is to be utilized. Only specific up to
                 Euclidean, D94, D2000. Do not use to indicate AVX/SSEetc.*/
};

/* g being global */
extern struct Parameters g_params;

#endif
