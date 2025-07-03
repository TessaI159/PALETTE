#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "Queue.h"

static inline size_t nearest_pow_2(int n) {
	return 1ULL << (int)ceil(log2(n));
}

bool Queue_create(struct Queue *q, size_t cap) {
	memset(q, 0, sizeof(struct Queue));
	cap	  = nearest_pow_2(cap);
	q->buffer = malloc(sizeof(void *) * cap);
	if (!q->buffer) {
		fprintf(stderr, "Could not malloc queue.\n");
		return false;
	}
	q->cap	= cap;
	q->head = 0;
	q->tail = 0;
	q->size = 0;
	q->mask = cap - 1;
	return true;
}

void Queue_destroy(struct Queue *q) {
	if (q->buffer) {
		free(q->buffer);
	}
	memset(q, 0, sizeof(struct Queue));
}

bool Enqueue(struct Queue *q, void *obj) {
	if (q->size == q->cap)
		return false;
	q->buffer[q->tail] = obj;
	q->tail		   = (q->tail + 1) & (q->mask);
	q->size++;
	return true;
}

void *Dequeue(struct Queue *q) {
	if (q->size == 0)
		return NULL;

	void *rt = q->buffer[q->head];
	q->head	 = (q->head + 1) & (q->mask);
	q->size--;
	return rt;
}

bool Queue_full(struct Queue *q) {
	return ((q->tail + 1) & (q->mask)) == q->head;
}
