#include <allegro.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include "ptask.h"
#include "qlearn.h"
#include "matrices.h"


// System State Constants
#define RESET   0
#define PLAY    1
#define PAUSE   2
#define STOP    3

// Q-learning Parameters Change Constants
#define NPARAM  5                  // Total Number of Possible Learning Parameters
#define STEP    0.1                // Increase/Decrease Step of Learning Parameters


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
	if(values[p] < 1)
		values[p] += STEP;
	pthread_mutex_unlock(&mux_parameter_values);
}

void dec_parameter_value(int p)
{
	pthread_mutex_lock(&mux_parameter_values);
	if(values[p]>0.01)
	   values[p] -= STEP;
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
	values[0] = ql_get_learning_rate();
	values[1] = ql_get_discount_factor();
	values[2] = ql_get_expl_decay();
	values[3] = ql_get_epsini();
	values[4] = ql_get_epsfin();
}

void set_qlearning_values()
{
	pthread_mutex_lock(&mux_parameter_values);
	ql_set_learning_rate(values[0]);
	ql_set_discount_factor(values[1]);
	ql_set_expl_decay(values[2]);
	ql_set_epsini(values[3]);
	ql_set_epsfin(values[4]);
	ql_set_epsilon(values[3]);
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
	if(i>=0 && i<4){
		pthread_mutex_lock(&mux_sys_state);
		sys_state=i;
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
			if(pg == 0){    // acceleratore non attivo
				set_dyn_dt(0.001);   // settiamo il passo di integrazione della dinamica
				//scale_fact = 100 è il fattore di scala per il periodo del qlearning rispetto al periodo del model 
				pt_set_period(4, 1000); // settiamo il periodo del model
				pt_set_period(3, 100000);  // settiamo il periodo del qlearning
				pt_set_deadline(4, 1000); // settiamo la deadline relativa del model
				pt_set_deadline(3, 100000);  // settiamo la deadline relativa del qlearning
			}else{          // acceleratore attivo
				set_dyn_dt(0.001);  // settiamo il passo di integrazione della dinamica
				//scale_fact = 10 è il fattore di scala per il periodo del qlearning rispetto al periodo del model 
				pt_set_period(4, 700); // settiamo il periodo del model
				pt_set_period(3, 70000);  // settiamo il periodo del qlearning
				pt_set_deadline(4, 700); // settiamo la deadline relativa del model
				pt_set_deadline(3, 7000);  // settiamo la deadline relativa del qlearning
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
