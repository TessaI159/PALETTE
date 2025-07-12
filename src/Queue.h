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
	size_t mask;
};

/**
 *@brief Initializes a circular queue with the given capacity
 *
 * This function sets up the internal structure of the queue and allocates
 * memory for its backing buffer. If the provided capacity is zero or negative a
 * default minimum capacity of 2 is used. Internally, the capacity is rounded up
 * to the nearest power of two - 1.
 *
 *@param q A pointer to a pre-allocated Queue struct to initializes
 *@param The desired capacity of the queue
 *
 *@return true if the queue was successfully initialized, false if not
 *
 *@note After a successful call:
 *       - 'q->cap' will be the least power of 2 >= the requested capacity
 *       - 'q->size' will be 0
 *       - The queue will be empty and ready for use with Queue_push
 *       - The memory allocated for q->buffer must be freed via Queue_destroy
 *
 *@note The contents of 'q' before the call are overwritten
 *
 *@note The true capacity of the queue is cap - 1 due to its circular nature
 *
 *@warning The behavior is undefined if 'q == NULL'
 **/
bool Queue_create(struct Queue *q, size_t cap);

/**
 *@brief Destroys a circular queue
 *
 * This function destroys the queue, as well as any items that may be left in
 * the queue when called. If the queue elements require a custom destructor, use
 * Queue_destroy_with instead
 *
 *@param q A pointer to the Queue struct to be destroyed
 *
 *@return Nothing
 **/
void Queue_destroy(struct Queue *q);

/**
 *@brief Destroys a circular queue with a customer destructor function
 *
 * This function destroys the queue, as well as any items that may be left in
 * the queue when called. Uses destroy to release any queue elements.
 *
 *@param q A pointer to the Queue struct to be destroyed
 *@param (*destroy) A custom destructor function
 *
 *@return Nothing
 *
 *@note calling Queue_destroy_with(q, free) is equivalent to calling
 *Queue_destroy(q);
 **/
void Queue_destroy_with(struct Queue *q, void (*destroy)(void *));

/**
 *@brief Push an item onto the queue
 *
 * This function pushes a single item onto the queue, assuming there is space
 *
 *@param q A pointer to the queue to which you want to push
 *@param obj A pointer to the object which you want to push to the queue
 *
 *@return true on success, false on failure
 *
 *@note obj will be cast to a void pointer
 **/
bool Queue_push(struct Queue *q, void *obj);

/**
 *@brief Pop an item from the queue
 *
 * This function pops a single item from the queue, assuming the queue is
 * non-empty
 *
 *@param q A pointer to the queue from which you want to pop
 *
 *@return A void pointer to the popped object on success, or NULL on failure
 *
 *@note Returns a void pointer, so the caller is responsible for casting to the
 * correct type
 **/
void *Queue_pop(struct Queue *q);

/**
*@brief Check if the queue is full
*
*@param q A pointer to the queue which you want to check
*
*@return true if the queue is full, false elsewise
*
**/
bool Queue_is_full(const struct Queue *q);

/**
*@brief Check if the queue is empty
*
*@param q A pointer to the queue which you want to check
*
*@return true if the queue is empty, false elsewise
*
**/
bool  Queue_is_empty(const struct Queue *q);

/**
*@brief Check the size of the queue
*
*@param q A pointer to the queue which you want to check
*
*@return A size_t which correlates to the number of elements on the queue
*
**/
size_t Queue_size(const struct Queue *q);

/**
*@brief Check the number of empty slots on the queue
*
*@param q A pointer to the queue which you want to check
*
*@return A size_t which correlates to the number of empty slots on the queue
*
**/
size_t Queue_empty_slots(const struct Queue *q);

#endif
