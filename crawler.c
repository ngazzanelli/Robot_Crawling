#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "qlearn.h"
#include "ptask.h"


// Qlearning Constants
#define NSTATES     49
#define NACTIONS    4

// System State Constants
#define RESET   0
#define PLAY    1
#define PAUSE   2
#define STOP    3

// Task identifier Constants
#define INTERPRETER   1
#define GRAPHIC     2
#define CRAWLER     3
#define MODEL       4

// Global definitions
#define TH1UP   0       // Action move up link 1
#define TH1DW   1       // Action move down link 1
#define TH2UP   2       // Action move up link 2
#define TH2DW   3       // Action move down link 2
#define TH1MAX  1.40    // max theta1 angle [rad]
#define TH1MIN  -0.7    // min theta1 angle [rad]
#define TH2MAX  1.05    // max theta2 angle [rad]
#define TH2MIN  -1.05   // min theta2 angle [rad]
#define DTH1    0.35    // theta1 quantization step [rad]
#define DTH2    0.35    // theta2 quantization step [rad]
#define PRI     10      // Tasks priority 
#define DL      20      // Tasks deadline [ms]
#define PER     20      // Tasks period [ms]
#define PER_D   1       // Dynamic task peirod [ms] 
#define RHIT    -10     // Reward for hitting limit angles
#define RSCALE  10      // Reward scale factor 


// Actual Robot State  (SI PUÒ INCLUDERE ANCHE MATRICES.H)
typedef struct {
    float q1;
    float q2;
    float q3;
    float q4;
    float q5;
    float q6;
    float energy;   
    float dt3;      
} state;

// Desired Robot State
typedef struct {
    float t1d;
    float t2d;
    int flag;
} target;


// Reward Management Struct
typedef struct {
    int state;
    int reward;
    int flag;
} rs_for_plot;


// Functions from other modules
extern int get_sys_state(int *s);
extern void set_qlearning_values();
extern void init_parameter_values();
extern void init_state();
extern void get_state(state* s);

// Application Tasks: main will activate them
extern void* dynamics(void*arg);
extern void* interpreter(void* arg);
extern void* manual_interpreter(void* arg);
extern void* update_graphic(void* arg);
extern void* update_graphic_DC(void* arg);


//Mutexes
static pthread_mutex_t mux_desired_joint = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_reward = PTHREAD_MUTEX_INITIALIZER;


// Static Variables
static float t1min, t2min;
static float t1max, t2max;
static float dt1, dt2;
static int n2;              // # of theta2 quantiations 
static target qd;
static rs_for_plot rw;


//-----------------------------------------------------
// The following Function initializes global Variables
//-----------------------------------------------------
void init_global_variables(){
    t1min = TH1MIN;
    t1max = TH1MAX;
    t2min = TH2MIN;
    t2max = TH2MAX;
    dt1 = DTH1;
    dt2 = DTH2;
    n2 = (t2max - t2min)/dt2 + 1;
}


//-----------------------------------------------------
// The following Functions manage desired joint angles
// decided from Qlearning algorithm
//-----------------------------------------------------
void reset_desired_joint()
{
    pthread_mutex_lock(&mux_desired_joint);
    qd.t1d = 0;
    qd.t2d = 0;
    qd.flag = 0;
    pthread_mutex_unlock(&mux_desired_joint);
}

void get_desired_joint(target* t){
  pthread_mutex_lock(&mux_desired_joint);
  t->t1d = qd.t1d;
  t->t2d = qd.t2d;
  t->flag = qd.flag;
  qd.flag = 0;
  pthread_mutex_unlock(&mux_desired_joint);
}


//-----------------------------------------------------
//?????????????????????????????????????????????????????
//-----------------------------------------------------
void get_rs_for_plot(rs_for_plot* t){
  pthread_mutex_lock(&mux_reward);
  t->state = rw.state;
  t->reward = rw.reward;  
  t->flag = rw.flag;
  rw.flag = 0;
  pthread_mutex_unlock(&mux_reward);
}

void set_rs_for_plot(int r,int s){
  pthread_mutex_lock(&mux_reward);
  rw.state = s;
  rw.reward = r;  
  rw.flag = 1;
  pthread_mutex_unlock(&mux_reward);
}


//-----------------------------------------------------
// The following Function converts joint angles to 
// a state so that Transition Matrix is not needed
//-----------------------------------------------------
int angles2state(float t1, float t2){

    int i, j;
    //printf("le variabilidi giunto valgono %f e %f\n", t1, t2);
    i = round((t1 - t1min)/dt1);
    j = round((t2 - t2min)/dt2);
    //printf("le variabilidi giunto sono approssimate a %f e %f\n", i*dt1+TH1MIN, j*dt2+TH2MIN);
    //printf("i=%d, j=%d, n2=%d, ret = %d\n", i, j, n2, i*n2+j);
    return i*n2 + j;
}


//-----------------------------------------------------
// The following Function computes reward for the
// actual state
//-----------------------------------------------------
int get_reward(int s, int snew, state robot){

    int r = 0;

    if(robot.dt3 > 0)
        r = round(13 * robot.dt3);
    else if(robot.dt3 < 0)
        r = round(20 * robot.dt3);
    
    else{
        if(robot.energy > 0) //old_dt3
            r = round(5 * robot.energy);
    }

    r -= 5;

    if (snew == s) 
        r += RHIT;        // hit the limit angle

    return r;
}


//-----------------------------------------------------
// The following Function computes the next state given
// the action a
//-----------------------------------------------------
int next_desired_state(int a){
    pthread_mutex_lock(&mux_desired_joint);
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
    qd.flag = 1;
    pthread_mutex_unlock(&mux_desired_joint);
    return angles2state(qd.t1d, qd.t2d);
}


//-----------------------------------------------------
// The following Function implements the Learning Loop
//----------------------------------------------------- 
void* qlearning(void* arg){

    printf("QLEARN: task started\n");
    int i;                              // thread index
    int s, snew, a, exec, r = 0, old_exec = 0;
    long step = 0;
    //float newerr = 0;
    state robot;

    i = pt_get_index(arg);
    pt_set_activation(i);
    ql_init(NSTATES, NACTIONS);
    init_global_variables();
    init_parameter_values();
    //printf("QLEARN: inizio ciclo while\n");
    float old_q1 = 0;
    float old_dt3 = 0;
    
    while (get_sys_state(&exec) != STOP){

        if(exec == PLAY){

            if(old_exec == RESET){
                set_qlearning_values();
                old_exec = PLAY;
            }

            //float period = pt_get_period(3);
            //printf("QLEARN: il mio periodo vale %f\n", period);
            //printf("QLEARN: sono dentro all'if\n");
            step++;

            get_state(&robot);

            old_dt3 = robot.dt3;
            robot.energy = old_dt3;

            robot.dt3 = robot.q1 - old_q1;
            old_q1 = robot.q1;

            //printf("Ottenuto stato attuale variabili di giunto del robot\n");
            s = angles2state(robot.q4, robot.q5);
            set_rs_for_plot(r,s);
            //printf("Quantizzato lo stato: s = %d\n", s);
            a = ql_egreedy_policy(s);
            //printf("Ottenuta l'azione\n");
            snew = next_desired_state(a);   //Questa funzione aggiorna 
                                            //anche le variabili di giunto desiderate "qd"
            //printf("QLEARN: Ottenuto il nuovo stato\n");
        }
            
        pt_deadline_miss(i);
        pt_wait_for_period(i);

        if(exec == PLAY){  
            r = get_reward(s, snew, robot);
            //printf("Ottenuto il reward r = %d\n", r);
            /*newerr =*/ ql_updateQ(s, a, r, snew);
            ql_copy_Q();
            //printf("Aggioranta matrice Q\n");
            //err +=  (newerr - err)/step;
            //if (step % 100 == 0)
                //ql_print_Qmatrix();
            if (step % 100 == 0)
                ql_reduce_exploration();
              
        }

        if(exec == RESET){
            old_exec = exec;
            ql_init(NSTATES, NACTIONS);
            ql_copy_Q();
            reset_desired_joint();
        }
        
    }
    printf("QLEARN: task finshed\n");
    return NULL;
}


int main(){
       
    int i;      //ris;
    int mode;	// Manual or Qlearning mode decided from user


    printf("MAIN: scegliere la modalità: 0 -> controllo manuale, 1 -> qlearning\n");
    scanf("%d", &mode);

    init_state();
    init_global_variables();

    
    printf("MAIN: creo il task di gestione della grafica\n");
    pt_task_create( update_graphic, GRAPHIC, PER*1000, DL*1000, PRI);
    //printf("con il risultato %d\n",ris);
    
    if(mode == 1){
        printf("MAIN: creo il task di interpretazione dei comandi\n");
        pt_task_create( interpreter, INTERPRETER, PER*1000, DL*1000, PRI);
        //printf("con il risultato %d\n",ris);
        printf("MAIN: creo il task di qlearning\n");
        pt_task_create( qlearning, CRAWLER, 100*PER_D*1000, 100*PER_D*1000, PRI); //occhio a quanto valgono T e DT in model.c
        //printf("con il risultato %d\n",ris);
    }
    else{
        printf("MAIN: creo il task di interpretazione dei comandi\n");
        pt_task_create( manual_interpreter, INTERPRETER, PER*1000, DL*1000, PRI);
        //printf("con il risultato %d\n",ris);
    }
    printf("MAIN: creo il task per la risoluzione della dinamica\n");
    pt_task_create( dynamics, MODEL, PER_D*1000, PER_D*1000, PRI);
    //printf("con il risultato %d\n",ris);
    ql_Q_from_file("./prova.txt");    

    for(i = 1; i <= 4; i++){
          pt_wait_for_end(i);
          printf("MAIN: fine ciclo %d\n", i);
    }
    
    return 0;
}