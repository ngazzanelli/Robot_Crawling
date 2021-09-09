#include "ptask.h"
#include <sched.h>
#include <time.h>
#include <pthread.h>
#include <math.h>

#define	NT	50			//Maximum Number of Tasks

struct task_par {
	int		arg;
	long	wcet;
	int		period;		// [micro sec]
	int 	deadline;	// [micro sec]
	int 	priority;
	int		dmiss;
	struct timespec	at;
	struct timespec dl;		
};

// Mutex 
static pthread_mutex_t mux_period[NT] = PTHREAD_MUTEX_INITIALIZER;

static struct task_par		tp[NT];
pthread_t					tid[NT];
int							policy = SCHED_FIFO;


//----------------------------
//	Time management functions
//----------------------------
void time_copy(struct timespec *td, struct timespec ts){

	td->tv_sec = ts.tv_sec;
	td->tv_nsec = ts.tv_nsec;
}

void time_add_us (struct timespec *t, int us){

	t->tv_sec += us/1000000;
	t->tv_nsec += (us%1000000)*1000;

	if (t->tv_nsec > (1000000000)) {
		t->tv_nsec -= 1000000000;
		t->tv_sec += 1;
	}
}

int time_cmp (struct timespec t1, struct timespec t2){

	if (t1.tv_sec > t2.tv_sec) return 1;
	if (t1.tv_sec < t2.tv_sec) return -1;
	if (t1.tv_nsec > t2.tv_nsec) return 1;
	if (t1.tv_nsec > t2.tv_nsec) return -1;
	return 0;
}

int time_diff_nsec (struct timespec t1, struct timespec t2){

	return (t1.tv_sec-t2.tv_sec)*pow(10, 9)+(t1.tv_nsec-t2.tv_nsec);
}


//----------------------------
//	Task management functions
//----------------------------
void pt_ptask_init (int scheduler){
	policy = scheduler;
}

int pt_task_create(void* (*task) (void*), int i, int period, int drel, int prio){
	pthread_attr_t		myatt;
	struct sched_param	mypar;
	int					tret;
	
	if (i>=NT) return -1;
	
	tp[i].arg = i;
	pthread_mutex_lock(&mux_period[i]);
	tp[i].period = period;
	tp[i].deadline = drel;
	pthread_mutex_unlock(&mux_period[i]);
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

int pt_get_index(void* arg){
	struct task_par	*tpar;

	tpar = (struct task_par *)arg;
	return tpar->arg;
}

int pt_get_period(int i){
	int per;

	pthread_mutex_lock(&mux_period[i]);
	per = tp[i].period;
	pthread_mutex_unlock(&mux_period[i]);
	return per;
}

void pt_set_period(int i, int per){
	pthread_mutex_lock(&mux_period[i]);
	tp[i].period = per;
	pthread_mutex_unlock(&mux_period[i]);
}

void pt_set_deadline(int i, int drel){	
	pthread_mutex_lock(&mux_period[i]);
	tp[i].deadline = drel;
	pthread_mutex_unlock(&mux_period[i]);
}

int pt_get_deadline(int i){
	int drel;

	pthread_mutex_lock(&mux_period[i]);
	drel = tp[i].deadline;
	pthread_mutex_unlock(&mux_period[i]);
	return drel;
}

int pt_get_dmiss(int i){
	return tp[i].dmiss;
}

void pt_set_activation(int i){
	struct timespec	t;
	int per, drel;

	clock_gettime(CLOCK_MONOTONIC, &t);
	time_copy(&(tp[i].at), t);
	time_copy(&(tp[i].dl), t);
	pthread_mutex_lock(&mux_period[i]);
	per = tp[i].period;
	drel = tp[i].deadline;
	pthread_mutex_unlock(&mux_period[i]);
	time_add_us(&(tp[i].at), per);
	time_add_us(&(tp[i].dl), drel);	
}

int pt_deadline_miss(int i){
	struct timespec	now;

	clock_gettime(CLOCK_MONOTONIC, &now);
	if (time_cmp(now, tp[i].dl) > 0) {
		tp[i].dmiss++;
		return 1;
	}
	return 0;
}

void pt_wait_for_period (int i){
	int per;

	clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &(tp[i].at), NULL);
	pthread_mutex_lock(&mux_period[i]);
	per = tp[i].period;
	pthread_mutex_unlock(&mux_period[i]);
	time_add_us(&(tp[i].at), per);
	time_add_us(&(tp[i].dl), per);
}

void pt_wait_for_end (int i){
	pthread_join(tid[i], NULL);
}