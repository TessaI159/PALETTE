#ifndef QUEUE_H
#define QUEUE_H
#include <stdbool.h>
#include <stddef.h>

#include <stdint.h>

struct Queue {
	size_t cap;
	void **buffer;
	size_t head;
	size_t tail;
	size_t size;
	size_t mask;
};

/* Creates a queue with capacity cap. Cap is rounded to the nearest power of 2
   >= cap. */
bool Queue_create(struct Queue *q, size_t cap);
/* Destroys the queue */
void Queue_destroy(struct Queue *q);
/* Adds an object to the queue if there is space. Returns true on success,
   false on failure */
bool Enqueue(struct Queue *q, void *obj);
/* Removes an item from the queue, if there is an item to return. Returns NULL
   if the queue is empty */
void *Dequeue(struct Queue *q);
bool  Queue_full(struct Queue *q);
bool  Queue_empty(struct Queue *q);
/* Returns the number of elements currently in the queue. */
size_t Queue_size(struct Queue *q);

#endif
