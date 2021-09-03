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
static int crawler_dl;

static pthread_mutex_t mux = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_reward = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_CR_dl = PTHREAD_MUTEX_INITIALIZER;
// Struttura per lo stato reale del robot
typedef struct {
    float q1;
    float q2;
    float q3;
    float q4;
    float q5;
    float q6;
    float energy;   //energia spesa utile per il reward
    float dt3;    //variazione dell'angolo della ruota per il reward  
} state;

// Struttura per lo stato desiderato del robot
typedef struct {
    float t1d;
    float t2d;
    int flag;
} target;

typedef struct {
    int state;
    int reward;
    int flag;
} rs_for_plot;

static target qd;
static rs_for_plot rw;

//Funzioni dall'interprete
extern int get_stop();
extern int get_pause();
extern int get_play();
//Funzioni dal modello
extern void get_state(state* s);

//funzione che incrementa la variabile delle deadline miss
//l'ho fatta così per non modificare la funzione di pthask
void inc_crawler_dl()
{
    pthread_mutex_lock(&mux_CR_dl);
    crawler_dl++;
    pthread_mutex_unlock(&mux_CR_dl);
}

//funzione che ottiene il valore delle deadline miss
void get_crawler_dl(int * dl_miss)
{
    pthread_mutex_lock(&mux_CR_dl);
    * dl_miss = crawler_dl;
    pthread_mutex_unlock(&mux_CR_dl);
}
void init_global_variables(){
    t1min = TH1MIN;
    t1max = TH1MAX;
    t2min = TH2MIN;
    t2max = TH2MAX;
    dt1 = DTH1;
    dt2 = DTH2;
    n2 = (t2max - t2min)/dt2 + 1;
}
/*
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
*/
void get_desired_joint(target* t){
  pthread_mutex_lock(&mux);
  t->t1d = qd.t1d;
  t->t2d = qd.t2d;
  t->flag = qd.flag;
  qd.flag = 0;
  pthread_mutex_unlock(&mux);
}
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

// Angles to state (non serve la matrice T di transizione dello stato)
int angles2state(float t1, float t2){

    int i, j;

    i = (t1 - t1min)/dt1;
    j = (t2 - t2min)/dt2;
    printf("i=%d, j=%d, n2=%d, ret = %d\n", i, j, n2, i*n2+j);
    return i*n2 + j;
}

int get_reward(int s, int snew, state robot){

    int r;

    r = (int)( ( (robot.dt3 > 0) ? 50 : -10) * RSCALE - 1);
    if (snew == s) 
        r += RHIT;        // hit the limit angle

    return r;
}

int next_desired_state(int a){
    pthread_mutex_lock(&mux);
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
    pthread_mutex_unlock(&mux);
    return angles2state(qd.t1d, qd.t2d);
}


// Learning loop 
void* qlearning(void* arg){

    printf("QLEARN: task started\n");
    int i;      // thread index
    int play;  // Serve per saltare entrambe le sezioni in 
                // cui è diviso il corpo del while per 
                // evitare di chiamare due volte get_pause()
    int s, a, r=0, snew;
    long step = 0;
    float newerr, err = 0;
    state robot;

    i = pt_get_index(arg);
    pt_set_activation(i);
    ql_init(NSTATES, NACTIONS);
    init_global_variables();
    printf("QLEARN: inizio ciclo while\n");
    while (!get_stop()){

        //Controllo se l'applicazione è in pausa
        play = get_play();
        if(play){
            //printf("QLEARN: sono dentro all'if\n");
            step++;

            get_state(&robot);
            //printf("Ottenuto stato attuale variabili di giunto del robot\n");
            s = angles2state(robot.q4, robot.q5);
            set_rs_for_plot(r,s);
            //printf("Quantizzato lo stato: s = %d\n", s);
            a = ql_egreedy_policy(s);
            //printf("Ottenuta l'azione\n");
            snew = next_desired_state(a);   //Questa funzione aggiorna 
                                            //anche le variabili di giunto desiderate "qd"
            //printf("Ottenuto il nuovo stato\n");
        }
            if(pt_deadline_miss(i))
                inc_crawler_dl();
            pt_wait_for_period(i);

        if(play){  
            r = get_reward(s, snew, robot);
            //printf("Ottenuto il reward r = %d\n", r);
            newerr = ql_updateQ(s, a, r, snew);
            ql_copy_Q();
            //printf("Aggioranta matrice Q\n");
            //err +=  (newerr - err)/step;
            if (step % 100 == 0)
                ql_print_Qmatrix();
            if (step % 1000 == 0)
                ql_reduce_exploration();
              
        }
    }
    printf("QLEARN: task finshed\n");
    return NULL;
}

//Funzioni dal model.c
extern void init_state();
extern void* dynamics(void*arg);
extern void* interface(void* arg);
extern void* update_graphic(void* arg);
extern void* update_graphic_DC(void* arg);



int main(){
       
    int i,ris;

    init_state();
    printf("MAIN: creo il task di interfaccia::\n");
    ris = pt_task_create( interface, 1, PER, DL, PRI);
    //printf("con il risultato %d\n",ris);
    printf("MAIN: creo il task di gestione della grafica\n");
    ris = pt_task_create( update_graphic, 2, PER, DL, PRI);
    //printf("con il risultato %d\n",ris);
    printf("MAIN: creo il task di qlearning\n");
    ris = pt_task_create( qlearning, 3, 100, DL, PRI);
    //printf("con il risultato %d\n",ris);
    printf("MAIN: creo il task per la risoluzione della dinamica\n");
    ris = pt_task_create( dynamics, 4, 1, DL, PRI);
    //printf("con il risultato %d\n",ris);
    for(i = 1; i <= 4; i++){
          pt_wait_for_end(i);
          printf("MAIN: fine ciclo %d\n", i);
    }
}