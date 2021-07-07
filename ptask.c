#include "ptask.h"
#include <sched.h>
#include <time.h>
#include <pthread.h>

#define	NT	50

struct task_par {
	int		arg;
	long	wcet;
	int		period;
	int 	deadline;
	int 	priority;
	int		dmiss;
	struct timespec	at;
	struct timespec dl;		
};

struct task_par		tp[NT];
pthread_t			tid[NT];
int					policy = SCHED_FIFO;
//----------------------------
//	Time management functions
//----------------------------

void time_copy(struct timespec *td, struct timespec ts)
{
	td->tv_sec = ts.tv_sec;
	td->tv_nsec = ts.tv_nsec;
}

void time_add_ms (struct timespec *t, int ms)
{
	t->tv_sec += ms/1000;
	t->tv_nsec += (ms%1000)*1000000;
	if (t->tv_nsec > (1000000000)) {
		t->tv_nsec -= 1000000000;
		t->tv_sec += 1;
	}
}

int time_cmp (struct timespec t1, struct timespec t2)
{
	if (t1.tv_sec > t2.tv_sec) return 1;
	if (t1.tv_sec < t2.tv_sec) return -1;
	if (t1.tv_nsec > t2.tv_nsec) return 1;
	if (t1.tv_nsec > t2.tv_nsec) return -1;
	return 0;
}

//----------------------------
//	Task management functions
//----------------------------

void pt_ptask_init (int scheduler) 
{
	policy = scheduler;
}

int pt_task_create(void* (*task) (void*), int i, int period, int drel, int prio)
{
pthread_attr_t		myatt;
struct sched_param	mypar;
int					tret;
	
	if (i>=NT) return -1;
	
	tp[i].arg = i;
	tp[i].period = period;
	tp[i].deadline = drel;
	tp[i].priority = prio;
	tp[i].dmiss = 0;
	
	pthread_attr_init(&myatt);
	pthread_attr_setinheritsched(&myatt, PTHREAD_EXPLICIT_SCHED);
	pthread_attr_setschedpolicy(&myatt, SCHED_FIFO);
	mypar.sched_priority = tp[i].priority;
	pthread_attr_setschedparam(&myatt, &mypar);
	tret = pthread_create(&tid[i], &myatt, task, (void *)(&tp[i]));
	return tret;
	
}

int pt_get_index(void* arg)
{
struct task_par	*tpar;

	tpar = (struct task_par *)arg;
	return tpar->arg;
}

int pt_get_period(int i)
{
	return tp[i].period;
}

int pt_get_dmiss(int i)
{
	return tp[i].dmiss;
}

void pt_set_activation(int i)
{
struct timespec	t;
	
	clock_gettime(CLOCK_MONOTONIC, &t);
	time_copy(&(tp[i].at), t);
	time_copy(&(tp[i].dl), t);
	time_add_ms(&(tp[i].at), tp[i].period);
	time_add_ms(&(tp[i].dl), tp[i].deadline);	
}

int pt_deadline_miss(int i)
{
struct timespec	now;

	clock_gettime(CLOCK_MONOTONIC, &now);
	if (time_cmp(now, tp[i].dl) > 0) {
		tp[i].dmiss++;
		return 1;
	}
	return 0;
}

void pt_wait_for_period (int i)
{
	clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &(tp[i].at), NULL);
	time_add_ms(&(tp[i].at), tp[i].period);
	time_add_ms(&(tp[i].dl), tp[i].period);
}

void pt_wait_for_end (int i)
{
	pthread_join(tid[i], NULL);
}