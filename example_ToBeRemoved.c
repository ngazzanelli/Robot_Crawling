#include <stdlib.h> 
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "qlearn.h"

// Global definitions
#define L1      0.5     // link1 length [m]
#define L2      0.5     // link2 length [m]
#define RH      0.4     // robot height [m]
#define RL      0.8     // robot length [m]
#define TH1UP   0       // action move up link 1
#define TH1DW   1       // action move down link 1
#define TH2UP   2       // action move up link 2
#define TH2DW   3       // action move down link 2
#define TH1MAX  1.55    // max theta1 angle [rad]
#define TH1MIN  -0.20   // min theta1 angle [rad]
#define TH2MAX  -0.7    // max theta2 angle [rad]
#define TH2MIN  -2.45   // min theta2 angle [rad]
#define DTH1    0.35    // theta1 quantization step [rad]
#define DTH2    0.35    // theta2 quantization step [rad]

#define MAXSTEP 10000000L    // max number of steps

//Global variables modified by graphics.c
int stop;
int view;

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

// State variables
typedef struct{
    float x, y;     // current end effector position [m, m]
    float ox, oy;   // old end effector position [m, m]
    float z;        // robot position [m]
    float dz;       // position increment [m]
    float th1;      // angle of link 1 [rad]
    float th2;      // angle of link 2 [rad]
    float space;    // space covered [m]
} state;
static state rob;

void init_state()
{
float t1, t12;

    rob.th1 = TH1MIN;
    rob.th2 = TH2MIN;
    t1 = rob.th1;
    t12 = rob.th1  + rob.th2;
    rob.x = rob.ox = rob.z + L1*cos(t1) + L2*cos(t12);
    rob.y = rob.oy = RH + L1*sin(t1) + L2*sin(t12);
    rob.z = rob.dz = rob.space = 0.0;
}

void compute_end_point()
{
float t1, t12;

    t1 = rob.th1;
    t12 = rob.th1  + rob.th2;
    rob.ox = rob.x;
    rob.oy = rob.y;
    rob.x = rob.z + L1*cos(t1) + L2*cos(t12);
    rob.y = RH + L1*sin(t1) + L2*sin(t12);
}

// Angles to state (non serve la matrice T di transizione dello stato)
int angles2state(float t1, float t2)
{
int i, j;

    i = (t1 - t1min)/dt1;
    j = (t2 - t2min)/dt2;
    //printf("i=%d, j=%d, n2=%d, ret = %d\n", i, j, n2, i*n2+j);
    return i*n2 + j;
}

// Action to angles: aggiorna gli angoli del robot in base all'azione
void action2angles(int a)
{
    switch(a){
        case TH1UP: rob.th1 += DTH1; break;
        case TH1DW: rob.th1 -= DTH1; break;
        case TH2UP: rob.th2 += DTH2; break;
        case TH2DW: rob.th2 -= DTH2; break;
        default: break;
    }
    if (rob.th1 > TH1MAX) rob.th1 = TH1MAX;
    if (rob.th1 < TH1MIN) rob.th1 = TH1MIN;
    if (rob.th2 > TH2MAX) rob.th2 = TH2MAX;
    if (rob.th2 < TH2MIN) rob.th2 = TH2MIN;
}

// Next state: aggiorna lo stato in base all'azione e ritorna il nuovo stato
int next_state(int a)
{
int s;

    action2angles(a);
    compute_end_point();
    rob.dz = 0;
    if ((rob.y <= 0) || (rob.oy <= 0)){
        rob.dz = rob.ox - rob.x;
        rob.space += rob.dz;
        //printf("ho fatto un passo %f\n", rob.dz);
        //printf("ho fatto uno spazio %f\n", rob.space);
    }
    s = angles2state(rob.th1, rob.th2);
    return s;
}

// Reward (non serve la matrice R dei reward )
#define HARDL   -0.02   // hard level coordinate [m]
#define RHIT    -4      // when hitting angle limit
#define RLOW    -20     // when going too low level
#define RMOVE   -1      // for each move

int get_reward(int s, int snew)
{
int r;  
    if (snew == s) r = RHIT;        // hit the limit angle
    else{
        //if(rob.dz != 0) printf("Sono avanzato di %f\n", rob.dz);
        r = rob.dz*100;             // proportional to progress
    }
    if (rob.y < HARDL){ 
        //printf("Sono andato troppo basso\n");
        r += RLOW;                  // going too low
    }
    r += RMOVE;                     // for each move
    //if(r > 0) printf("reward=%d\n", r); //per debug vediamo quando dà reward positivi
    return r;           
}

// Extern graphic functions
extern void init_graphics(state s);
extern void display_links(state s);
extern void terminate_graphics();
extern void read_key();
extern void update_info(long step, float space);

// Learning loop 
float qlearn()
{
int s, a, r, snew;
long step = 0;
float newerr, err = 0;
    s = angles2state(rob.th1, rob.th2);
    //printf("Inizio il ciclo while dallo stato s = %d\n", s);
    while (!stop){
        step++;
        a = ql_egreedy_policy(s);
        //printf("Ottenuta l'azione\n");
        snew = next_state(a);
        //printf("Ottenuto il nuovo stato\n");
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
        // Monitoriamo la presenza del ciclo sugli stati
        //if (step > MAXSTEP*0.99)
        //    printf("Stato s = %d\n", s);
        //printf("Finito ciclo %ld\n", step);
        if(view){
            display_links(rob);
            //usleep(1000);
        }
        read_key();
    }
    return err;
}

int main()
{
    srand(time(NULL));
    init_global_variables();
    init_state();
    init_graphics(rob);
    ql_init(36, 4);  //36 states, 4 actions
    ql_set_learning_rate(0.5);
    ql_set_discount_factor(0.9);
    ql_set_expl_range(1.0, 0.1);
    ql_set_expl_decay(0.95);
    //display_menu();
    //display_param();
    //display_environment();
    //interpreter();
    printf("inizio il qlearning\n");
    qlearn();
    //ql_print_Qmatrix();
    terminate_graphics();
    return 0;
}