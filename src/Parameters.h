#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdbool.h>

enum WorkingSpace {
  WORKING_OKLAB,
  WORKING_CIELAB
};

struct Parameters {
  bool chroma_weight;
  enum WorkingSpace working_colorspace;
};

extern struct Parameters global_parameters;

#endif
