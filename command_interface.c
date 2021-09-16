#include <allegro.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include "ptask.h"
#include "qlearn.h"
#include "matrices.h"

// Filename for saving/loading Q-Matrix
#define Q_MATRIX_FILE	"./q_matrix.txt"

// System State Constants
#define RESET   0
#define PLAY    1
#define PAUSE   2
#define STOP    3

// Q-learning Parameters Change Constants
#define NPARAM  5                  // Total Number of Possible Learning Parameters
#define STEP    0.01               // Increase/Decrease Step of Learning Parameters
#define ALPHA		0
#define GAMMA		1
#define	DECAY		2
#define EPS_MAX	3
#define EPS_MIN	4

// Functions from other modules
extern void init_state();
extern int next_desired_state(int a);
extern void set_dyn_dt(float dt);
extern float get_dyn_dt();


// Mutexes
static pthread_mutex_t mux_sys_state			= PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_parameter_values     = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_parameter_selected	= PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_pause_graphic		= PTHREAD_MUTEX_INITIALIZER;


// Static Variables
static int		sys_state;
static int		pause_graphic;
static int		parameter_selected;		// Selected Parameter
static float	values[5];				// Values of QLearning Parameters
//  value[0] -> alpha
//  value[1] -> gamma
//  value[2] -> decay
//  value[3] -> initial epsilon
//  value[4] -> final epsilon


//-----------------------------------------------------
// The following Functions manage the Qlearning 
// Parameters Selection from user
//-----------------------------------------------------
void inc_parameter_selected()
{
  pthread_mutex_lock(&mux_parameter_selected);
  parameter_selected = (parameter_selected + 1) % NPARAM;
  pthread_mutex_unlock(&mux_parameter_selected);
}

void dec_parameter_selected()
{
  pthread_mutex_lock(&mux_parameter_selected);
  parameter_selected = (parameter_selected + 4) % NPARAM;
  pthread_mutex_unlock(&mux_parameter_selected);
}

int get_parameter_selected()
{
  int ret;
  pthread_mutex_lock(&mux_parameter_selected);
  ret = parameter_selected;
  pthread_mutex_unlock(&mux_parameter_selected);
  return ret;
}


//-----------------------------------------------------
// The following Functions manage the change of 
// Qlearning Parameters from user
//-----------------------------------------------------
void inc_parameter_value(int p)
{
	pthread_mutex_lock(&mux_parameter_values);
	//printf("GRAPHIC: il parametro %d valeva %f ", p, values[p]);
	values[p] += STEP;
	if(values[p] > 1) 
		values[p] = 1;
	if(p == EPS_MIN){	//epsilon min
		if(values[p] > values[EPS_MAX])
			values[p] = values[EPS_MAX];
	}
	//printf("e ora vale %f\n", values[p]);
	pthread_mutex_unlock(&mux_parameter_values);
}

void dec_parameter_value(int p)
{
	pthread_mutex_lock(&mux_parameter_values);
	//printf("il parametro %d valeva %f ", p, values[p]);
	values[p] -= STEP;
	if(values[p] < 0) 
		values[p] = 0;
	if(p == EPS_MAX){//epsilon max
		if(values[p] < values[EPS_MIN])
			values[p] = values[EPS_MIN];
	}
	//printf("e ora vale %f\n", values[p]);
	pthread_mutex_unlock(&mux_parameter_values);
}

void get_parameter_values(float *buff)
{
	pthread_mutex_lock(&mux_parameter_values);
	vector_copy(values, buff, 5);
	pthread_mutex_unlock(&mux_parameter_values);
}

void init_parameter_values()
{
	values[ALPHA] = ql_get_learning_rate();
	values[GAMMA] = ql_get_discount_factor();
	values[DECAY] = ql_get_expl_decay();
	values[EPS_MAX] = ql_get_epsini();
	values[EPS_MIN] = ql_get_epsfin();
}

void set_qlearning_values()
{
	pthread_mutex_lock(&mux_parameter_values);
	ql_set_learning_rate(values[ALPHA]);
	ql_set_discount_factor(values[GAMMA]);
	ql_set_expl_decay(values[DECAY]);
	ql_set_epsini(values[EPS_MAX]);
	ql_set_epsfin(values[EPS_MIN]);
	ql_set_epsilon(values[EPS_MAX]);
	pthread_mutex_unlock(&mux_parameter_values);
}


//-----------------------------------------------------
// The following Functions manage the
// application acceleration
//-----------------------------------------------------
void change_pause_graphic()
{
	pthread_mutex_lock(&mux_pause_graphic);
	pause_graphic = (pause_graphic)?0:1;
	pthread_mutex_unlock(&mux_pause_graphic);
}

int get_pause_graphic()
{
	int temp;

	pthread_mutex_lock(&mux_pause_graphic);
	temp = pause_graphic;
	pthread_mutex_unlock(&mux_pause_graphic);
	return temp;
}

//-----------------------------------------------------
// The following Functions manage and track 
// the System State
//-----------------------------------------------------
void set_sys_state(int i)
{
	if(i >= 0 && i < 4){
		pthread_mutex_lock(&mux_sys_state);
		sys_state = i;
		pthread_mutex_unlock(&mux_sys_state);
	}
}
int get_sys_state(int *s)
{
	pthread_mutex_lock(&mux_sys_state);
	*s = sys_state;
	pthread_mutex_unlock(&mux_sys_state);
	return *s;
}


//------------------------------------------------------
//  The following function waits for a key pressed and 
//  extracts the corresponding ascii code and scan code
//------------------------------------------------------
char get_scancode()
{
	if(keypressed())
		return readkey()>>8;
	else    
		return 0;
}


//-------------------------------------------------------------
// The following Function manages the interaction with the user
// and the state of the system. 
//------------------------------------------------------------- 
void key_manager(int exec)
{
	int p;          // Tracks Qlearning parameters change
	int pg;         // Tracks pause_graphics
	char cm;

	cm = get_scancode();

	switch(cm){

		/* Writing/Reading Q-Matrix to/from file */
		case KEY_F:
		printf("INTERPRETER: hai premuto il tasto F\n");
			if(exec == PAUSE){
				ql_Q_to_file(Q_MATRIX_FILE);
				printf("INTERPRETER: Salvata matrice Q su file\n");
			}
			break;

		case KEY_L:
		printf("INTERPRETER: hai premuto il tasto L\n");
			if(exec == RESET){
				if(ql_Q_from_file(Q_MATRIX_FILE) == 1)
					printf("INTERPRETER: Letta matrice Q da file\n");
			}
			break;
		/*---------------------------------------*/ 	

		case KEY_R:
			printf("INTERPRETER: hai premuto il tasto R\n");
			set_sys_state(RESET);
			init_state();
			break;

		case KEY_S:
			printf("INTERPRETER: hai premuto il tasto S\n");
			if(exec == RESET)
				set_sys_state(PLAY);
			break;

		case KEY_P:
			printf("INTERPRETER: hai premuto il tasto P\n");
			if(exec == PLAY)
				set_sys_state(PAUSE);
			else if(exec == PAUSE)
				set_sys_state(PLAY);
			break;

		case KEY_E:
			printf("INTERPRETER: hai premuto il tasto E\n");
			set_sys_state(STOP);
			break;

		case KEY_UP:
			if(sys_state == RESET){
				printf("INTERPRETER: hai premuto il tasto UP\n");
				dec_parameter_selected();
			}
			break;

		case KEY_DOWN:
			if(sys_state == RESET){
				printf("INTERPRETER: hai premuto il tasto DOWN\n");
				inc_parameter_selected();
			}
			break;

		case KEY_RIGHT:
			if(sys_state == RESET){
				printf("INTERPRETER: hai premuto il tasto RIGHT\n");
				p = get_parameter_selected();
				inc_parameter_value(p);
			}
			break;

		case KEY_LEFT:
			if(sys_state == RESET){
				printf("INTERPRETER: hai premuto il tasto LEFT\n");
				p = get_parameter_selected();
				dec_parameter_value(p);
			}
			break;

		case KEY_B:
			printf("INTERPRETER: hai premuto il tasto G\n");
			change_pause_graphic();
			pg = get_pause_graphic(); 

			if(pg == 0){    			// accelerator disabled
				set_dyn_dt(0.001);  // setting dynamic integration dt
				// setting qlearn task parameters
				pt_set_period(3, 100000);  
				pt_set_deadline(3, 100000);
				// setting model task parameters
				pt_set_period(4, 1000); 
				pt_set_deadline(4, 1000);

			} else{          			// accelerator enabled
				set_dyn_dt(0.001);  // setting dynamic integration dt
				// setting qlearn task parameters
				pt_set_period(3, 60000);
				pt_set_deadline(3, 6000);
				// setting model task parameters
				pt_set_period(4, 600); 
				pt_set_deadline(4, 600); 
				  
			}
		default: break;
	}
}


//---------------------------------------
// User Interface Task
// --------------------------------------
void* interpreter(void * arg)
{
	printf("INTERPRETER: task started\n");
	int i,  exec;

	i = pt_get_index(arg);
	pt_set_activation(i);

	set_sys_state(RESET);
	init_parameter_values();

	while(get_sys_state(&exec) != STOP){
		key_manager(exec);
		pt_deadline_miss(i);
		pt_wait_for_period(i);
	}

	printf("INTERPRETER: task finished\n");
	return NULL;
}


//---------------------------------------
// Manual Mode Functions
// --------------------------------------
#define TH1UP   0       // Action move up link 1
#define TH1DW   1       // Action move down link 1
#define TH2UP   2       // Action move up link 2
#define TH2DW   3       // Action move down link 2

void key_manager_manual(int exec)
{        
	
	char cm;
	cm = get_scancode();

	switch(cm){

		case KEY_R:
			printf("INTERPRETER: hai premuto il tasto R\n");
			set_sys_state(RESET);
			break;

		case KEY_S:
			printf("INTERPRETER: hai premuto il tasto S\n");

			if(exec == RESET){
				init_state();
				set_sys_state(1);
			}
			break;

		case KEY_P:
			printf("INTERPRETER: hai premuto il tasto P\n");

			if(exec == PAUSE)
				set_sys_state(PLAY);
			else if(exec == PLAY)
				set_sys_state(PAUSE);

			break;

		case KEY_E:
			printf("INTERPRETER: hai premuto il tasto E\n");
			set_sys_state(STOP);  
			break;

		case KEY_UP:
			printf("INTERPRETER: hai premuto il tasto UP\n");
			next_desired_state(TH1UP);
			break;

		case KEY_DOWN:
			printf("INTERPRETER: hai premuto il tasto DOWN\n");
			next_desired_state(TH1DW);
			break;

		case KEY_RIGHT:
			printf("INTERPRETER: hai premuto il tasto RIGHT\n");
			next_desired_state(TH2UP);
			break;

		case KEY_LEFT:
			printf("INTERPRETER: hai premuto il tasto LEFT\n");
			next_desired_state(TH2DW);
			break;

		default: break;
	}
}


void* manual_interpreter(void* arg)
{
	int i,  exec;

	printf("INTERPRETER: task started\n");
	i = pt_get_index(arg);
	pt_set_activation(i);

	set_sys_state(RESET);
	init_parameter_values();

	while(get_sys_state(&exec) != STOP){
		key_manager_manual(exec);
		pt_deadline_miss(i);
		pt_wait_for_period(i);

	}

	printf("INTERPRETER: task finished\n");
	return NULL;
}
