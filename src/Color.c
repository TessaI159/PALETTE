#include "Color.h"
#include "Align.h"
#include <string.h>

bool Color_alloc(struct Color *c, uint32_t num)
{
    if (!c)
        return false;
    c->num = 0;
    c->alpha = NULL;
    c->beta = NULL;
    c->gamma = NULL;
    if (num == 0)
        return true;

    size_t size = (size_t)num * sizeof(float);
    c->alpha = (float *)aligned_malloc(32, size);
    c->beta = (float *)aligned_malloc(32, size);
    c->gamma = (float *)aligned_malloc(32, size);

    if (!c->alpha || !c->beta || !c->gamma) {
        Color_free(c);
        return false;
    }

    c->num = num;
    return true;
}

void Color_free(struct Color *c)
{
    if (!c)
        return;
    if (c->alpha) {
        aligned_free(c->alpha);
        c->alpha = NULL;
    }
    if (c->beta) {
        aligned_free(c->beta);
        c->beta = NULL;
    }
    if (c->gamma) {
        aligned_free(c->gamma);
        c->gamma = NULL;
    }
    c->num = 0;
}

bool Color_clone(struct Color *dest, const struct Color *src)
{
    if (!dest || !src)
        return false;
    if (!Color_alloc(dest, src->num))
        return false;

    if (src->num > 0) {
        size_t size = (size_t)src->num * sizeof(float);
        memcpy(dest->alpha, src->alpha, size);
        memcpy(dest->beta, src->beta, size);
        memcpy(dest->gamma, src->gamma, size);
    }

    dest->num = src->num;
    return true;
}
