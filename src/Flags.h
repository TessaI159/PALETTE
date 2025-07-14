#ifndef FLAGS_H
#define FLAGS_H
#include <stdbool.h>

struct Flags {
  bool chroma_average;
};

extern struct Flags global_flags;

#endif
