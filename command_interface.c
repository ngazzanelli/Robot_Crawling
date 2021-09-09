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

// Qlearning Parameters Change Constants
#define NPARAM  5               // Total Number of Possible Learning Parameters
#define STEP 0.1                // Increase/Decrease Step of Learning Parameters


// Functions from other modules
extern void init_state();
extern int next_desired_state(int a);
extern void set_dyn_dt(float dt);
extern float get_dyn_dt();


// Mutex 
static pthread_mutex_t mux_sys_state = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_parameter_values = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_parameter_selected = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_int_dl = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_pause_graphic = PTHREAD_MUTEX_INITIALIZER;


// Static Variables
static int sys_state;
static int pause_graphic;
static int interface_dl;        // Interface Deadline Miss
static int parameter_selected;  // Selected Parameter
static float values[5];         // Values of QLearning Parameters
//  value[0] -> alpha
//  value[1] -> gamma
//  value[2] -> decay
//  value[3] -> epsilon iniziale
//  value[4] -> epsilon finale



//-----------------------------------------------------
// The following Functions manage the Qlearning 
// Parameters Selection from user
//-----------------------------------------------------
void inc_parameter_selected(){
  pthread_mutex_lock(&mux_parameter_selected);
  parameter_selected = (parameter_selected + 1) % NPARAM;
  pthread_mutex_unlock(&mux_parameter_selected);
}

void dec_parameter_selected(){
  pthread_mutex_lock(&mux_parameter_selected);
  parameter_selected = (parameter_selected + 4) % NPARAM;
  pthread_mutex_unlock(&mux_parameter_selected);
}

int get_parameter_selected(){
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
void inc_parameter_value(int p){
    pthread_mutex_lock(&mux_parameter_values);
    if(values[p] < 1)
        values[p] += STEP;
    pthread_mutex_unlock(&mux_parameter_values);
}

void dec_parameter_value(int p){
    pthread_mutex_lock(&mux_parameter_values);
    if(values[p]>0.01)
       values[p] -= STEP;
    pthread_mutex_unlock(&mux_parameter_values);
}

void get_parameter_values(float *buff){
    pthread_mutex_lock(&mux_parameter_values);
    vector_copy(values, buff, 5);
    pthread_mutex_unlock(&mux_parameter_values);
}

void init_parameter_values(){
    values[0] = ql_get_learning_rate();
    values[1] = ql_get_discount_factor();
    values[2] = ql_get_expl_decay();
    values[3] = ql_get_epsini();
    values[4] = ql_get_epsfin();
}

void set_qlearning_values(){
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
// ????????????????????????????????????????????????????
//-----------------------------------------------------
void change_pause_graphic(){
    pthread_mutex_lock(&mux_pause_graphic);
    pause_graphic = (pause_graphic)?0:1;
    pthread_mutex_unlock(&mux_pause_graphic);
}
int get_pause_graphic(){
    int temp;
    pthread_mutex_lock(&mux_pause_graphic);
    temp = pause_graphic;
    pthread_mutex_unlock(&mux_pause_graphic);
    return temp;
}


//-----------------------------------------------------
// The following Functions manage Interface Task 
// Deadline Misses
//-----------------------------------------------------
void get_interface_dl(int * dl_miss)
{
    pthread_mutex_lock(&mux_int_dl);
    * dl_miss = interface_dl;
    pthread_mutex_unlock(&mux_int_dl);
}


void inc_interface_dl()
{
    pthread_mutex_lock(&mux_int_dl);
    interface_dl++;
    pthread_mutex_unlock(&mux_int_dl);
}


//-----------------------------------------------------
// The following Functions manage the System State
//-----------------------------------------------------
void set_sys_state(int i)
{
    if(i>=0 && i<4)
    {
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


// Funzione per ottenere il codice del tasto premuto da tastiera 
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
    int p;          //Serve per gestire il cambio di parametro del qlearning
    char cm;
    int pg;         //Serve per sapere il valore di pause_graphic

    cm = get_scancode();
    switch(cm){

        case KEY_R:
            printf("INTERFACE: hai premuto il tasto R\n");
            set_sys_state(RESET);
            init_state();
            break;

        case KEY_S:
            printf("INTERFACE: hai premuto il tasto S\n");
            if(exec == RESET)
                set_sys_state(PLAY);
            break;

        case KEY_P:
            printf("INTERFACE: hai premuto il tasto P\n");
            if(exec == PLAY)
                set_sys_state(PAUSE);
            else if(exec == PAUSE)
                set_sys_state(PLAY);
            break;

        case KEY_E:
            printf("INTERFACE: hai premuto il tasto E\n");
            set_sys_state(STOP); 
            break;

        case KEY_UP:
            if(sys_state == RESET){
                printf("INTERFACE: hai premuto il tasto UP\n");
                dec_parameter_selected();
            }
            break;

        case KEY_DOWN:
            if(sys_state == RESET){
                printf("INTERFACE: hai premuto il tasto DOWN\n");
                inc_parameter_selected();
            }
            break;

        case KEY_RIGHT:
            if(sys_state == RESET){
                printf("INTERFACE: hai premuto il tasto RIGHT\n");
                p = get_parameter_selected();
                inc_parameter_value(p);
            }
            break;

        case KEY_LEFT:
            if(sys_state == RESET){
                printf("INTERFACE: hai premuto il tasto LEFT\n");
                p = get_parameter_selected();
                dec_parameter_value(p);
            }
            break;
        case KEY_B:
            printf("hai premuto il tasto G\n");
            change_pause_graphic();
            pg = get_pause_graphic(); //pg = "pause graphic"
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
void* interface(void * arg)
{
    printf("INTERPRETER: task started\n");
    int i,  exec;
    i = pt_get_index(arg);
    pt_set_activation(i);
    set_sys_state(RESET);

    init_parameter_values();

    while(get_sys_state(&exec) != STOP)
    {
        key_manager(exec);
        if(pt_deadline_miss(i))
            inc_interface_dl();
        pt_wait_for_period(i);
    }
    printf("INTERPRETER: task end\n");
    return NULL;
}






// Manual Mode

void key_manager_manual(int *exec){
    int p;          //Serve per gestire il cambio di parametro del qlearning
    float step;     //Dice di quanto incrementare/decrementare value
    char cm;
    int s;

/*DA RIFARE TOGLIENDO GET_STOP ETC
    cm = get_scancode();
    switch(cm){

        case KEY_R:
            printf("INTERFACE: hai premuto il tasto R\n");
            set_sys_state(0);
            //TODO: gestire il reset dei vari task
            break;

        case KEY_S:
            printf("INTERFACE: hai premuto il tasto S\n");
            if(get_reset()){
                init_state();
                //init_graphics();
                //init_qlearning();
                set_sys_state(1);
            }
            break;

        case KEY_P:
            printf("INTERFACE: hai premuto il tasto P\n");
            if(get_play())
                set_sys_state(2);
            else if(get_pause())
                set_sys_state(1);
            break;

        case KEY_E:
            printf("INTERFACE: hai premuto il tasto E\n");
            set_sys_state(3);  
            *exec=0;
            break;

        case KEY_UP:
            printf("INTERFACE: hai premuto il tasto UP\n");
            s = next_desired_state(0);
            printf("sei andato nello stato %d\n", s);
            break;

        case KEY_DOWN:
            printf("INTERFACE: hai premuto il tasto DOWN\n");
            s = next_desired_state(1);
            printf("sei andato nello stato %d\n", s);
            break;

        case KEY_RIGHT:
            printf("INTERFACE: hai premuto il tasto RIGHT\n");
            s = next_desired_state(2);
            printf("sei andato nello stato %d\n", s);
            break;

        case KEY_LEFT:
            printf("INTERFACE: hai premuto il tasto LEFT\n");
            s = next_desired_state(3);
            printf("sei andato nello stato %d\n", s);
            break;

        default: break;
    }*/
}



void* manual_interface(void * arg)
{
    printf("INTERPRETER: task started\n");
    int i,  exec;
    i = pt_get_index(arg);
    pt_set_activation(i);
    set_sys_state(RESET);

    while(exec)
    {
        key_manager_manual(&exec);
        if(pt_deadline_miss(i))
            inc_interface_dl();
        pt_wait_for_period(i);
        if(exec==0)
        {
            printf("ho finito \n");
        }
    }
    printf("INTERPRETER: task end\n");
    return NULL;
}
