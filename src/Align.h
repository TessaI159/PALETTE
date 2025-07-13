#ifndef ALIGN_H
#define ALIGN_H

#if defined(_MSC_VER)
#define ALIGN16 __declspec(align(16))
#else
#define ALIGN16 __attribute__ ((aligned(16)))
#endif


#endif
