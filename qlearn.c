
#include "qlearn.h"
#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>
#include <math.h>

// Global constants
#define MAXSTA	100		// max # of states
#define MAXACT	100		// max # of actions
#define ALPHA0	1.0		// default learning rate
#define EPSINI	0.9 	// initial exploration factor
#define EPSFIN	0.1 	// final exploration factor
#define GAMMA0	0.9 	// default discount factor
#define DECAY0	0.95	// default epsilon decay rate

// QL matrices
//static int T[MAXSTA][MAXACT];	// transition matrix	INUTILE
//static int R[MAXSTA][MAXACT];	// reward matrix		INUTILE
static float Q[MAXSTA][MAXACT];	// Q matrix
static float Q_copy[MAXSTA][MAXACT];	// Q matrix for graphic
										// (si usa per non accedere
										// a troppe sezioni critiche 
										// ogni volta che il crawler
										// sta imparando)
static pthread_mutex_t mux_Qmatrix = PTHREAD_MUTEX_INITIALIZER;

// Global variables
static int nsta;		// actual # of states
static int nact;		// actual number of actions
//static int goalstate;	// store the goal state			INUTILE
static float alpha;		// learning rate
static float gam;		// discount factor
static float epsilon; 	// actual exploration probability
static float decay;		// decay rate for epsilon
static float eps_norm;	// normalized exploration probability
static float eps_ini;	// initial exploration probability
static float eps_fin;	// final exploration probability


// Auxiliary functions
float frand(float xmin, float xmax)
{
float range;
	range = (xmax - xmin);
	return (xmin + range*(float)rand()/RAND_MAX);
}

void ql_copy_Q(){
	int s,a;
	pthread_mutex_lock(&mux_Qmatrix);
	for (s=0; s<nsta; s++)
		for(a=0; a<nact; a++)
			Q_copy[s][a] = Q[s][a];
	pthread_mutex_unlock(&mux_Qmatrix);
}

void ql_get_Q(float *dest){//da verificare che in dest ci sia abbastanza spazio
	int s,a;
	pthread_mutex_lock(&mux_Qmatrix);
	for (s=0; s<nsta; s++)
		for(a=0; a<nact; a++)
			dest[s*4+a] = Q_copy[s][a];
	pthread_mutex_unlock(&mux_Qmatrix);
}


// Initialization
void ql_init(int ns, int na)
{
int s, a;
	
	nsta = ns;
	nact = na;

	if(nsta > MAXSTA){
		printf("NUmber of states too big\n");
		exit(1);
	}
	if(nact > MAXACT){
		printf("Number of actions too big\n");
		exit(1);
	}

	alpha = ALPHA0;
	gam = GAMMA0;
	eps_norm = 1.0;
	eps_ini = EPSINI;
	eps_fin = EPSFIN;
	decay = DECAY0;

	for (s=0; s<nsta; s++)
		for(a=0; a<nact; a++)
			Q[s][a] = 0;
}

// Set parameters functions
void ql_set_learning_rate(float a)
{
	alpha = a;
}
void ql_set_discount_factor(float g)
{
	gam = g;
}
void ql_set_expl_range(float e_ini, float e_fin)
{
	eps_ini = e_ini;
	eps_fin = e_fin;
}
void ql_set_epsini(float e_ini)
{
	eps_ini = e_ini;
}
void ql_set_epsfin(float e_fin)
{
	eps_fin = e_fin; 
}
void ql_set_expl_decay(float d)
{
	decay = d;
}

// Get parameters functions
float ql_get_learning_rate()
{
	return alpha;
}
float ql_get_discount_factor()
{
	return gam;
}
float ql_get_epsini()
{
	return eps_ini; 
}
float ql_get_epsfin()
{
	return eps_fin; 
}
float ql_get_expl_decay()
{
	return decay;
}

// Reduce exploration
void ql_reduce_exploration()
{
	eps_norm = decay*eps_norm;
	epsilon = eps_fin + eps_norm*(eps_ini - eps_fin);
	//printf("Ho ridotto epsilon = %f\n", epsilon);

}

// Maximum Q value in a given state s
float ql_maxQ(int s)
{
int a;
float m;

	m = Q[s][0];	// initialized with Q value for action 0
	for (a=1; a<nact; a++)
		if (Q[s][a] > m)
			m = Q[s][a];
	return m;
}

// Best action in a given state s
float ql_best_action(int s)
{
int a, ba;
float m;

	ba = 0;			// initialized best action with action 0
	m = Q[s][0];	// initialized with Q value for action 0

	for (a=1; a<nact; a++)
		if (Q[s][a] > m){
			m = Q[s][a];
			ba = a;
		}
	return ba;
}

// Epsilon-greedy policy in a given state s
int ql_egreedy_policy (int s)
{
int ra, ba;
float x;

	ba = ql_best_action(s);
	ra = rand()%nact;
	x = frand(0, 1);
	if (x < epsilon) return ra;
	else return ba;
}

// Update Q value 
float ql_updateQ(int s, int a, int r, int snew)
{
float Qtarget, TDerr;

	//printf("r = %d, s = %d, a = %d, snew = %d, gamma = %f\n", r, s, a, snew, gam);
	Qtarget = r + gam*ql_maxQ(snew);
	//printf("Qtarget = %f\n", Qtarget);
	TDerr = Qtarget - Q[s][a];
	//printf("TDerr = %f\n", TDerr);
	Q[s][a] = Q[s][a] + alpha*TDerr;
	//printf("Q[%d][%d] = %f\n", s, a, Q[s][a]);
	return fabs(TDerr);
}

// Show the Q matrix
void ql_print_Qmatrix()
{
int i, j;
	printf("Q matrix:\n");
	for(i=0; i<nsta; i++){
		for(j=0; j<nact; j++)
			printf("%f ", Q[i][j]);
		printf("\n");
	}
}