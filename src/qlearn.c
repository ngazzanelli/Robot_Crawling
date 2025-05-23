#include "qlearn.h"
#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>
#include <math.h>

// Global constants
#define MAXSTA	100		// max # of states
#define MAXACT	100		// max # of actions
#define ALPHA0	0.1		// default learning rate
#define EPSINI	0.9 	// initial exploration factor
#define EPSFIN	0.01 	// final exploration factor
#define GAMMA0	0.9 	// default discount factor
#define DECAY0	0.95	// default epsilon decay rate


//Mutex
static pthread_mutex_t mux_Qmatrix = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_eps = PTHREAD_MUTEX_INITIALIZER;

// Global variables
static int nsta;						// actual # of states
static int nact;						// actual number of actions
static float alpha;						// learning rate
static float gam;						// discount factor
static float epsilon; 					// actual exploration probability
static float decay;						// decay rate for epsilon
static float eps_norm;					// normalized exploration probability
static float eps_ini;					// initial exploration probability
static float eps_fin;					// final exploration probability
static float Q[MAXSTA][MAXACT];			// Q matrix
static float Q_copy[MAXSTA][MAXACT];	// Q matrix for graphic task, used for 
										// avoiding useless critical sections 


//-----------------------------------------------------
// The following Function initializes all the variables
// needed to implemetn QLearning Algorithm
//-----------------------------------------------------
void ql_init(int ns, int na){
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
	epsilon = EPSINI;

	for (s = 0; s < nsta; s++)
		for(a = 0; a < nact; a++)
			Q[s][a] = 0;

	// IDEAL CYCLE TEST 
	 /*Q[9][0] = 10;
	 Q[16][0] = 10;
	 Q[23][2] = 10;
	 Q[24][2] = 10;
	 Q[25][1] = 10;
	 Q[18][1] = 10;
	 Q[11][3] = 10;
	 Q[10][3] = 10;*/
}


//-----------------------------------------------------
// The following Functions manage Q matrix and its copy
//-----------------------------------------------------
void ql_copy_Q(){
	int s,a;
	pthread_mutex_lock(&mux_Qmatrix);
	for (s = 0; s < nsta; s++)
		for(a = 0; a < nact; a++)
			Q_copy[s][a] = Q[s][a];
	pthread_mutex_unlock(&mux_Qmatrix);
}

void ql_get_Q(float *dest){
	int s, a;
	pthread_mutex_lock(&mux_Qmatrix);
	for (s = 0; s < nsta; s++)
		for(a = 0; a < nact; a++)
			dest[s*4 + a] = Q_copy[s][a];
	pthread_mutex_unlock(&mux_Qmatrix);
}


//-----------------------------------------------------
// The following Functions set the parameters 
// used in the QLearning Algorithm
//-----------------------------------------------------
void ql_set_learning_rate(float a){
	alpha = a;
}

void ql_set_discount_factor(float g){
	gam = g;
}

void ql_set_expl_range(float e_ini, float e_fin){
	eps_ini = e_ini;
	eps_fin = e_fin;
}

void ql_set_epsini(float e_ini){
	eps_ini = e_ini;
}

void ql_set_epsfin(float e_fin){
	eps_fin = e_fin; 
}

void ql_set_epsilon(float e){
	epsilon = e;
}

void ql_set_expl_decay(float d){
	decay = d;
}


//-----------------------------------------------------
// The following Functions get the parameters 
// used in the QLearning Algorithm
//-----------------------------------------------------
float ql_get_learning_rate(){
	return alpha;
}

float ql_get_discount_factor(){
	return gam;
}

float ql_get_epsini(){
	return eps_ini; 
}

float ql_get_epsfin(){
	return eps_fin; 
}

float ql_get_epsilon(){
	float ret;
	pthread_mutex_lock(&mux_eps);
	ret = epsilon;
	pthread_mutex_unlock(&mux_eps);
	return ret;
}

float ql_get_expl_decay(){
	return decay;
}


//-----------------------------------------------------
// The following Function reduces exploration factor 
// epsilon, depending on decay, inital and final
// epsilon values
//-----------------------------------------------------
void ql_reduce_exploration(){

	eps_norm = decay*eps_norm;
	pthread_mutex_lock(&mux_eps);
	epsilon = eps_fin + eps_norm*(eps_ini - eps_fin);
	pthread_mutex_unlock(&mux_eps);
	//printf("Ho ridotto epsilon = %f\n", epsilon);
	ql_get_epsilon();

}


//-----------------------------------------------------
// The following Function computes the Maximum Q
// value in a given state s
//-----------------------------------------------------
float ql_maxQ(int s){
int a;
float m;

	m = Q[s][0];	// initialized with Q value for action 0
	for (a = 1; a < nact; a++)
		if (Q[s][a] > m)
			m = Q[s][a];
	return m;
}


//-----------------------------------------------------
// The following Function computes the best action
// in a given state s
//-----------------------------------------------------
float ql_best_action(int s){
	int a, ba;
	float m;

	ba = 0;			// initialized best action with action 0
	m = Q[s][0];	// initialized with Q value for action 0

	for (a = 1; a < nact; a++)
		if (Q[s][a] > m){
			m = Q[s][a];
			ba = a;
		}
	return ba;
}


//-----------------------------------------------------
// The following Function computes Epsilon-greedy
// policy in a given state s
//-----------------------------------------------------
int ql_egreedy_policy (int s){
int ra, ba;
float x, eps;

	ba = ql_best_action(s);
	ra = rand()%nact;
	//printf("L'azione casuale è %d\n", ra);
	//printf("Epsilon = %f\n", epsilon);
	x = frand(0, 1);
	pthread_mutex_lock(&mux_eps);
	eps = epsilon;
	pthread_mutex_unlock(&mux_eps);
	if (x < eps) return ra;
	else return ba;
}


//-----------------------------------------------------
// The following Function updates Q value 
//-----------------------------------------------------
float ql_updateQ(int s, int a, int r, int snew){
	float Qtarget, TDerr;

	Qtarget = r + gam*ql_maxQ(snew);
	TDerr = Qtarget - Q[s][a];
	Q[s][a] = Q[s][a] + alpha*TDerr;
	return fabs(TDerr);
}


//-----------------------------------------------------
// The following Function shows the Q matrix
//-----------------------------------------------------
void ql_print_Qmatrix()
{
	int i, j;
	printf("Q matrix:\n");
	for(i = 0; i < nsta; i++){
		for(j = 0; j < nact; j++)
			printf("%f ", Q[i][j]);
		printf("\n");
	}
}


//-----------------------------------------------------
// The following Function saves Q-Matrix in the file
// named "filename"; if this file doesn't exist, it's
// created.
//-----------------------------------------------------
void ql_Q_to_file(char* filename)
{
	int i, j;
	FILE* fptr;

	fptr = fopen(filename, "w");

	if(fptr == NULL){
		printf("QLEARN: Error in opening file\n");
		return;
	}

	for(i = 0; i < nsta; i++)
		for(j = 0; j < nact; j++)
			fprintf(fptr, "%d %f\n", i*nact + j, Q[i][j]);

	fclose(fptr);

	return;

}


//-----------------------------------------------------
// The following Function loads Q-Matrix from the file
// named "filename"; if this doesn't exist, the function
// prints an Error Message.
//-----------------------------------------------------
int ql_Q_from_file(char* filename)
{
	int i, j, index;
	FILE* fptr;
		

	fptr = fopen(filename, "r");

	if(fptr == NULL){
		printf("QLEARN: Error, file: %s doesn't exist!\n", filename);
		return 0;
	}

	for(i = 0; i < nsta; i++)
		for(j = 0; j < nact; j++)
			fscanf(fptr, "%d %f\n", &index, &Q[i][j]);
			  
	fclose(fptr);

	return 1;
}


//-----------------------------------------------------
// THe following Function computes a random number 
//-----------------------------------------------------
float frand(float xmin, float xmax){
	float range;
	range = (xmax - xmin);
	return (xmin + range*(float)rand()/RAND_MAX);
}
