#include <pthread.h>
#include "qlearn.h"
#include "ptask.h"

#define NSTATES     36
#define NACTIONS    4

// Global definitions
#define TH1UP   0       // action move up link 1
#define TH1DW   1       // action move down link 1
#define TH2UP   2       // action move up link 2
#define TH2DW   3       // action move down link 2
#define TH1MAX  1.55    // max theta1 angle [rad]
#define TH1MIN  -0.55   // min theta1 angle [rad]
#define TH2MAX  1.05    // max theta2 angle [rad]
#define TH2MIN  -1.05   // min theta2 angle [rad]
#define DTH1    0.35    // theta1 quantization step [rad]
#define DTH2    0.35    // theta2 quantization step [rad]

// Global variables (inutili se tanto ho le costanti, no?)
static float t1min, t2min;
static float t1max, t2max;
static float dt1, dt2;
static int n2;          // # of theta2 quantiations 

void init_global_variables()
{
    t1min = TH1MIN;
    t1max = TH1MAX;
    t2min = TH2MIN;
    t2max = TH2MAX;
    dt1 = DTH1;
    dt2 = DTH2;
    n2 = (t2max - t2min)/dt2 + 1;
}
//-------------------------------------------

// Struttura da condividere con model.c
typedef struct {
    float q1d;
    float q2d;
    int flag;
} target;
static target qd;

static pthread_mutex_t mux = PTHREAD_MUTEX_INITIALIZER;

void set_desired_joint();
void get_desired_joint(target* t){
  pthread_mutex_lock(&mux);
  t->q1d = qd.q1d;
  t->q2d = qd.q2d;
  t->flag = qd.flag;
  qd.flag = 0;
  pthread_mutex_unlock(&mux);
}
//-----------------------------------

// Funzione definita in model.c
extern float get_energy();
//-----------------------------------

// Struttura per lo stato del robot
typedef struct {
    float q1;
    float q2;
    float q3;
    float q4;
    float q5;
    float q6;
    float energy;
} state;

// Angles to state (non serve la matrice T di transizione dello stato)
int angles2state(float t1, float t2)
{
int i, j;

    i = (t1 - t1min)/dt1;
    j = (t2 - t2min)/dt2;
    //printf("i=%d, j=%d, n2=%d, ret = %d\n", i, j, n2, i*n2+j);
    return i*n2 + j;
}

extern int get_stop();
extern get_state(state* s);

// Learning loop 
void* qlearning(void* arg)
{
int i; 
int s, a, r, snew;
long step = 0;
float newerr, err = 0;
state robot;

    i = pt_get_index(arg);
    pt_set_activation(i);

    get_state(&robot);
    s = angles2state(robot.q4, robot.q5);
    //printf("Inizio il ciclo while dallo stato s = %d\n", s);
    while (!get_stop()){

        step++;
        a = ql_egreedy_policy(s);
        //printf("Ottenuta l'azione\n");
        snew = next_state(a);
        //printf("Ottenuto il nuovo stato\n");
        
        pt_deadline_miss(i);
        pt_wait_for_period(i);
        
        r = get_reward(s, snew);
        //printf("Ottenuto il reward\n");
        newerr = ql_updateQ(s, a, r, snew);
        //printf("Aggioranta matrice Q\n");
        err +=  (newerr - err)/step;
        s = snew;
        if (step%1000 == 0){
            ql_reduce_exploration();
            update_info(step, rob.space);
        }
    }
    return NULL;
}


int main()
{
int i;

    srand(time(NULL));
    init_global_variables();
    init_state();
    ql_init(NSTATES, NACTIONS);  //36 states, 4 actions
    task_create(dynamics, 1, PER, DL, PRI);
    task_create(graphics, 2, PER, DL, PRI);
    task_create(interpreter, 3, PER, DL, PRI);
    task_create(qlearning, 4, PER, DL, PRI);
    for(i=1; i<=4; i++){
		  wait_for_task_end(i);
		  //printf("fine ciclo %d\n", i);
    }
}