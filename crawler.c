#include <stdlib.h>
#include <stdio.h>
#include "qlearn.h"
#include "ptask.h"

#define NSTATES     49
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
#define PRI     10      // tasks priority 
#define DL      20      // tasks deadline
#define PER     20      // tasks period
#define RHIT    -10     // reward for hitting limit angles
#define RSCALE  100     // reward scale factor 

// Global variables (inutili se tanto ho le costanti, no?)
static float t1min, t2min;
static float t1max, t2max;
static float dt1, dt2;
static int n2;          // # of theta2 quantiations 

static pthread_mutex_t mux = PTHREAD_MUTEX_INITIALIZER;

// Struttura per lo stato reale del robot
typedef struct {
    float q1;
    float q2;
    float q3;
    float q4;
    float q5;
    float q6;
    float energy;   //energia spesa utile per il reward
    float dth3;     //variazione dell'angolo della ruota per il reward  
} state;

// Struttura per lo stato desiderato del robot
typedef struct {
    float t1d;
    float t2d;
    int flag;
} target;
static target qd;


extern int get_stop();
extern void get_state(state* s);


void init_global_variables(){
    t1min = TH1MIN;
    t1max = TH1MAX;
    t2min = TH2MIN;
    t2max = TH2MAX;
    dt1 = DTH1;
    dt2 = DTH2;
    n2 = (t2max - t2min)/dt2 + 1;
}

void set_desired_joint(int s){
    int i, j;    
    
    i = s/n2;
    j = s%n2;
    pthread_mutex_lock(&mux);
    qd.t1d = t1min + i*dt1;
    qd.t2d = t2min + j*dt2;
    qd.flag = 1;
    pthread_mutex_unlock(&mux);

}

void get_desired_joint(target* t){
  pthread_mutex_lock(&mux);
  t->t1d = qd.t1d;
  t->t2d = qd.t2d;
  t->flag = qd.flag;
  qd.flag = 0;
  pthread_mutex_unlock(&mux);
}


// Angles to state (non serve la matrice T di transizione dello stato)
int angles2state(float t1, float t2){

    int i, j;

    i = (t1 - t1min)/dt1;
    j = (t2 - t2min)/dt2;
    //printf("i=%d, j=%d, n2=%d, ret = %d\n", i, j, n2, i*n2+j);
    return i*n2 + j;
}

int get_reward(int s, int snew, state robot){

    int r;

    r = (int)(robot.dth3 * RSCALE - robot.energy);
    if (snew == s) 
        r += RHIT;        // hit the limit angle

    return r;
}

int next_desired_state(int a){
    switch(a){
        case TH1UP: 
            qd.t1d += DTH1;
            break;
        case TH1DW:
            qd.t1d -= DTH1;
            break;
        case TH2UP:
            qd.t2d += DTH2;
            break;
        case TH2DW: 
            qd.t2d -= DTH2;
            break;
        default: 
            break;
    }
    if (qd.t1d > TH1MAX) qd.t1d = TH1MAX;
    if (qd.t1d < TH1MIN) qd.t1d = TH1MIN;
    if (qd.t2d > TH2MAX) qd.t2d = TH2MAX;
    if (qd.t2d < TH2MIN) qd.t2d = TH2MIN;

    return angles2state(qd.t1d, qd.t2d);
}


// Learning loop 
void* qlearning(void* arg){
    printf("qlearning task started\n");
    int i;      // thread index
    int s, a, r, snew;
    long step = 0;
    float newerr, err = 0;
    state robot;

    i = pt_get_index(arg);
    pt_set_activation(i);


    while (1 /*!get_stop()*/){
        step++;

        get_state(&robot);
        s = angles2state(robot.q4, robot.q5);
        a = ql_egreedy_policy(s);
        //printf("Ottenuta l'azione\n");
        snew = next_desired_state(a);
        //printf("Ottenuto il nuovo stato\n");

        pt_deadline_miss(i);
        pt_wait_for_period(i);
        
        r = get_reward(s, snew, robot);
        //printf("Ottenuto il reward\n");
        newerr = ql_updateQ(s, a, r, snew);
        //printf("Aggioranta matrice Q\n");
        //err +=  (newerr - err)/step;
        if (step % 1000 == 0)
            ql_reduce_exploration();
    }
    printf("qlearning task finshed\n");
    return NULL;
}

//extern void* interpreter(void* arg);
extern void* dynamics(void*arg);
//extern void* graphics(void* arg);

int main(){
       
    int i;

    srand(time(NULL));
    init_global_variables();
    //init_state();
    ql_init(NSTATES, NACTIONS);  //49 states, 4 actions

    pt_task_create(dynamics, 1, PER, DL, PRI);
    //pt_task_create(graphics, 2, PER, DL, PRI);
    //pt_task_create(interpreter, 3, PER, DL, PRI);
    pt_task_create(qlearning, 4, PER, DL, PRI);

    for(i = 1; i <= 4; i++){
		  pt_wait_for_end(i);
		  printf("fine ciclo %d\n", i);
    }
}