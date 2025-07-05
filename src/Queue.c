#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "Queue.h"

/* TODO (Tess): Account for overflow*/
static inline size_t nearest_pow_2(int n) {
	if (n <= 2)
		return 2;
	return 1ULL << (int)ceil(log2(n));
}

bool Queue_create(struct Queue *q, size_t cap) {
	if (cap <= 0) {
		cap = 2;
	}

	cap	  = nearest_pow_2(cap);
	q->buffer = malloc(sizeof(void *) * cap);
	if (!q->buffer) {
		fprintf(stderr, "Could not malloc queue.\n");
		return false;
	}
	memset(q->buffer, 0, sizeof(void *) * cap);
	q->cap	= cap;
	q->head = 0;
	q->tail = 0;
	q->mask = cap - 1;
	return true;
}

void Queue_destroy(struct Queue *q) {
	while (!Queue_is_empty(q)) {
		void *elem = Queue_pop(q);
		if (elem != NULL)
			free(elem);
	}
	if (q->buffer) {
		free(q->buffer);
	}
	memset(q, 0, sizeof(struct Queue));
}

void Queue_destroy_with(struct Queue *q, void (*destroy)(void *)) {
	while (!Queue_is_empty(q)) {
		void *elem = Queue_pop(q);
		if (elem != NULL)
			destroy(elem);
	}
	if (q->buffer) {
		free(q->buffer);
	}
	memset(q, 0, sizeof(struct Queue));
}

bool Queue_push(struct Queue *q, void *obj) {
	if (obj == NULL)
		return false;
	if (((q->tail + 1) & (q->mask)) == q->head)
		return false;
	q->buffer[q->tail] = obj;
	q->tail		   = (q->tail + 1) & (q->mask);
	return true;
}

void *Queue_pop(struct Queue *q) {
	if (q->tail == q->head)
		return NULL;
	void *rt = q->buffer[q->head];
	q->head	 = (q->head + 1) & (q->mask);
	return rt;
}

bool Queue_is_full(const struct Queue *q) {
	return ((q->tail + 1) & (q->mask)) == q->head;
}

bool Queue_is_empty(const struct Queue *q) {
	return (q->tail == q->head);
}

size_t Queue_empty_slots(const struct Queue *q) {
	return ((q->head - q->tail - 1) & q->mask);
}

size_t Queue_size(const struct Queue *q) {
	return ((q->tail - q->head) & q->mask);
}
