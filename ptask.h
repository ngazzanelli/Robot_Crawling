#ifndef	PTHASK_H
#define	PTHASK_H

#include <time.h>
#include <pthread.h>

//	Time management functions
extern void time_copy (struct timespec *td, struct timespec ts);
extern void time_add_us (struct timespec *t, int ms);
extern int time_cmp (struct timespec t1, struct timespec t2);
int time_diff_nsec (struct timespec t1, struct timespec t2);

//	Task management functions
extern void pt_ptask_init(int scheduler);
extern int pt_task_create(void* (*task) (void*), int i, int period, int drel, int prio);
extern int pt_get_index(void* arg);
extern int pt_get_period(int i);
extern int pt_get_deadline(int i);
extern void pt_set_period(int i, int per);
extern void pt_set_deadline(int i, int drel);
extern int pt_get_dmiss(int i);
extern void pt_set_activation(int i);
extern int pt_deadline_miss(int i);
extern void pt_wait_for_period (int i);
extern void pt_wait_for_end(int i);
#endif