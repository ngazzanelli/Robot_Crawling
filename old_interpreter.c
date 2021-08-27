#include <pthread.h>
#include <allegro.h>
#include "ptask.h"
#include "qlearn.h"

#define NPARAM  5  //total number of possible learning parameters
#define STEP_ALPHA  1
#define STEP_GAMMA  0.1
#define STEP_DECAY   0.2
#define STEP_EPS  0.3

static int stop;  
static int pause;    
static int no_graphics;
static int parameter_selected;  //dice qual è il parametro attualmente selezionato
//Possibili valori di parameter_selected:
//  0 -> alpha
//  1 -> gamma
//  2 -> decay
//  3 -> epsilon iniziale
//  4 -> epsilon finale
static pthread_mutex_t mux_stop = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_pause = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_no_graphics = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mux_parameter_selected = PTHREAD_MUTEX_INITIALIZER;

void set_stop(){
  pthread_mutex_lock(&mux_stop);
  stop = 1;
  pthread_mutex_unlock(&mux_stop);
}
void change_pause(){
  pthread_mutex_lock(&mux_pause);
  pause = 1 - pause;
  pthread_mutex_unlock(&mux_pause);
}
void inc_parameter_selected(){
  pthread_mutex_lock(&mux_parameter_selected);
  parameter_selected = (parameter_selected+1)%NPARAM;
  pthread_mutex_unlock(&mux_parameter_selected);
}
void dec_parameter_selected(){
  pthread_mutex_lock(&mux_parameter_selected);
  parameter_selected = (parameter_selected+4)%NPARAM;
  pthread_mutex_unlock(&mux_parameter_selected);
}
void change_no_graphics(){
  pthread_mutex_lock(&mux_no_graphics);
  no_graphics = 1 - no_graphics;
  pthread_mutex_unlock(&mux_no_graphics);
}
int get_stop(){
  int ret;
  pthread_mutex_lock(&mux_stop);
  ret = stop;
  pthread_mutex_unlock(&mux_stop);
  return ret;
}
int get_pause(){
  int ret;
  pthread_mutex_lock(&mux_pause);
  ret = pause;
  pthread_mutex_unlock(&mux_pause);
  return ret;
}
int get_parameter_selected(){
  int ret;
  pthread_mutex_lock(&mux_parameter_selected);
  ret = parameter_selected;
  pthread_mutex_unlock(&mux_parameter_selected);
  return ret;
}
int get_no_graphics(){
  int ret;
  pthread_mutex_lock(&mux_no_graphics);
  ret = no_graphics;
  pthread_mutex_unlock(&mux_no_graphics);
  return ret;
}

void* interpreter (void *arg){
  int i;
  char scan;
  int resetting = 0;  //dice se siamo in modalità di reset e quindi possono essere cambiati i parametri del learning
  int p;
  float value;  //dice qual è il valore del parametro attualmente selezionato
  float step;    //dice di quanto incrementare/decrementare value

  i = get_task_index(arg);
  set_activation(i);  
  printf("Interprete: inizio (task %d)\n", i);

    do{
        scan = get_scancode_nb();
        switch (scan){
        case KEY_P:  // pause/play
            change_pause();
            break;
        case KEY_G:  // active/deactive graphic mode
            change_no_graphics();
            break;
        case KEY_R:  // reset the application
            resetting = 1;
            //resettare tutto
            break;
        case KEY_S:  // start the application (after reset)
            resetting = 0;
            //far ripartire tutto 
            break;
        case KEY_UP:
            if(resetting == 1){
            // Per prima cosa salviamo il valore precedentemente selezionato
            p = get_parameter_selected(); //si salva in locale per problemi di mutua esclusione
            // Passiamo quindi al successivo parametro
            inc_parameter_selected();
            //Possibili valori di parameter_selected:
            //  0 -> alpha  1 -> gamma  2 -> decay  3 -> eps ini  4 -> eps fin
            switch(p){
                case 0:  
                ql_set_learning_rate(value);
                step = STEP_GAMMA;
                //value = ...
                break;
                case 1:
                ql_set_discount_factor(value);
                step = STEP_DECAY;
                break;
                case 2:
                ql_set_expl_decay(value);
                step = STEP_EPS;
                break;
                case 3:
                ql_set_epsini(value);
                break; 
                case 4:
                ql_set_epsfin(value);
                step = STEP_ALPHA;
                break;
                default: break;
            }
            }
            break;
        case KEY_DOWN:
            if(resetting == 1){
            // Per prima cosa salviamo il valore precedentemente selezionato
            p = get_parameter_selected();
            // Passiamo quindi al successivo parametro
            dec_parameter_selected();
            //Possibili valori di parameter_selected:
            //  0 -> alpha  1 -> gamma  2 -> decay  3 -> eps ini  4 -> eps fin
            switch(p){
                case 0:  
                ql_set_learning_rate(value);
                step = STEP_EPS;
                break;
                case 1:
                ql_set_discount_factor(value);
                step = STEP_ALPHA;
                break;
                case 2:
                ql_set_expl_decay(value);
                step = STEP_GAMMA;
                break;
                case 3:
                ql_set_epsini(value);
                step = STEP_DECAY;
                break; 
                case 4:
                ql_set_epsfin(value);
                break;
                default: break;
            }
            }
            break;
        case KEY_RIGHT:
            if(resetting == 1)
            value += step;
            break;
        case KEY_LEFT:
            if(resetting == 1)
            value -= step;
            break;
        default: break;
        }
        deadline_miss(i);
        wait_for_period(i);
    } while (scan != KEY_ESC);
  
  set_stop();

  printf("Interprete: ho finito\n");
  return NULL;
}