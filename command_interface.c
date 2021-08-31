#include <allegro.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include "ptask.h"
#include "qlearn.h"

#define PI 3.14
#define PASSO 0.01

// Costanti utili per la modifica dei parametri del qlearning
#define NPARAM  5  //total number of possible learning parameters
#define STEP_ALPHA  1
#define STEP_GAMMA  0.1
#define STEP_DECAY   0.2
#define STEP_EPS  0.3

static int parameter_selected;  //Dice qual è il parametro attualmente selezionato
//Possibili valori di parameter_selected:
//  0 -> alpha
//  1 -> gamma
//  2 -> decay
//  3 -> epsilon iniziale
//  4 -> epsilon finale
static pthread_mutex_t mux_parameter_selected = PTHREAD_MUTEX_INITIALIZER;

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
int get_parameter_selected(){
  int ret;
  pthread_mutex_lock(&mux_parameter_selected);
  ret = parameter_selected;
  pthread_mutex_unlock(&mux_parameter_selected);
  return ret;
}

//per provare sia l'interprete che la capacità di scrivere 
//i igusti plot mi sono definito uno stato "falso" in cui va 
//a lavorare il generatore di onde; nella grafica basta inserire
//la dichiarazione esterna della get_state reale, inserire 
//nel corpo del thread la dichiarazione della struttura e 
//dovrebbe fungere lo stesso
static pthread_mutex_t mux_false_state = PTHREAD_MUTEX_INITIALIZER;

float FALSE_STATE[6];

static int com_state;
static pthread_mutex_t mux_int = PTHREAD_MUTEX_INITIALIZER;
// Variabile globale  che controlla il flusso delle operazioni:
   //case 0: sistema appena acceso, unico modo per settare
   //   le variabili del qlearning
   //case 1: sistema in play, i task sono attivi ed eseguono 
   //   tutte le istruzioni 
   //case 2: sistema in pausa, i task rimangono attivi 
   //   ma non eseguono il flusso principale
   //case 3: arresto del sistema 

int get_pause(){
    int temp ;
    pthread_mutex_lock(&mux_int);
    temp = com_state;
    pthread_mutex_unlock(&mux_int);
    if(temp == 2)
        return 1;
    else
        return 0;
}


int get_play(){
    int temp ;
    pthread_mutex_lock(&mux_int);
    temp = com_state;
    pthread_mutex_unlock(&mux_int);
    if(temp == 1)
        return 1;
    else
        return 0;
}


int get_stop(){
    int temp ;
    pthread_mutex_lock(&mux_int);
    temp = com_state;
    pthread_mutex_unlock(&mux_int);
    if(temp == 3)
        return 1;
    else
        return 0;
}

int get_reset(){
    int temp ;
    pthread_mutex_lock(&mux_int);
    temp = com_state;
    pthread_mutex_unlock(&mux_int);
    if(temp == 0)
        return 1;
    else
        return 0;
}

//inizio get set per lo stato falso 
void set_FALSE_ST(float *i)
{
    int j;
    pthread_mutex_lock(&mux_false_state);
    for(j=0;j<6;j++)
        FALSE_STATE[j]=i[j];
    pthread_mutex_unlock(&mux_false_state);
}

void get_FALSE_ST(float *i)
{
    int j;
    pthread_mutex_lock(&mux_false_state);
    for(j=0;j<6;j++)
        i[j]=FALSE_STATE[j];
    pthread_mutex_unlock(&mux_false_state);
}

//funzioni di interfaccia per accedere allo stato dell'applicazione
void set_com_variable(int i)
{
    if(i>=0 && i<4)
    {
        pthread_mutex_lock(&mux_int);
        com_state=i;
        pthread_mutex_unlock(&mux_int);
    }
}
void get_com_variable(int *ic) //LA MODIFICHEREI METTENDO IL RITORNO DI TIPO INTERO E TOGLIENDO IL PARAMETRO PASSATO
{
    pthread_mutex_lock(&mux_int);
    *ic = com_state;
    pthread_mutex_unlock(&mux_int);
}

//Funzioni esterne da altri moduli
extern void init_state();
//extern voif init_graphics();
//extern void init_qlearning();

// Funzione per ottenere il codice del tasto premuto da tastiera 
char get_scancode()
{
    if(keypressed())
        return readkey()>>8;
    else    
        return 0;
}

void init_com_inter()
{
    //inizializzazione dell'interprete
    set_com_variable(0);
}

float value; //Dice qual è il valore del parametro attualmente selezionato
void key_manager(int *exec)
{
    int p;          //Serve per gestire il cambio di parametro del qlearning
    float step;     //Dice di quanto incrementare/decrementare value
    char cm;


    cm = get_scancode();
    switch(cm){

        case KEY_R:
            printf("INTERFACE: hai premuto il tasto R\n");
            set_com_variable(0);
            //TODO: gestire il reset dei vari task
            break;

        case KEY_S:
            printf("INTERFACE: hai premuto il tasto S\n");
            if(get_reset()){
                init_state();
                //init_graphics();
                //init_qlearning();
                set_com_variable(1);
            }
            break;

        case KEY_P:
            printf("INTERFACE: hai premuto il tasto P\n");
            if(get_play())
                set_com_variable(2);
            else if(get_pause())
                set_com_variable(1);
            break;

        case KEY_E:
            printf("INTERFACE: hai premuto il tasto E\n");
            set_com_variable(3);  
            *exec=0;
            break;

        case KEY_UP:
            if(com_state == 0){
                printf("INTERFACE: hai premuto il tasto UP\n");
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
                        value = ql_get_discount_factor();
                        break;
                    case 1:
                        ql_set_discount_factor(value);
                        step = STEP_DECAY;
                        value = ql_get_expl_decay();
                        break;
                    case 2:
                        ql_set_expl_decay(value);
                        step = STEP_EPS;
                        value = ql_get_epsini();
                        break;
                    case 3:
                        ql_set_epsini(value);
                        value = ql_get_epsfin();
                        break; 
                    case 4:
                        ql_set_epsfin(value);
                        step = STEP_ALPHA;
                        value = ql_get_learning_rate();
                        break;
                    default: break;
                }
            }
            break;

        case KEY_DOWN:
            if(com_state == 0){
                printf("INTERFACE: hai premuto il tasto DOWN\n");
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
                        value = ql_get_epsfin();
                        break;
                    case 1:
                        ql_set_discount_factor(value);
                        step = STEP_ALPHA;
                        value = ql_get_learning_rate();
                        break;
                    case 2:
                        ql_set_expl_decay(value);
                        step = STEP_GAMMA;
                        value = ql_get_discount_factor();
                        break;
                    case 3:
                        ql_set_epsini(value);
                        step = STEP_DECAY;
                        value = ql_get_expl_decay();
                        break; 
                    case 4:
                        ql_set_epsfin(value);
                        value = ql_get_epsini();
                        break;
                    default: break;
                }
            }
            break;

        case KEY_RIGHT:
            if(com_state == 0){
                printf("INTERFACE: hai premuto il tasto RIGHT\n");
                value += step;
            }
            break;

        case KEY_LEFT:
            if(com_state == 0){
                printf("INTERFACE: hai premuto il tasto LEFT\n");
                value -= step;
            }
            break;

        default: break;
    }
}

void* interface(void * arg)
{
    printf("INTERPRETER: task started\n");
    int i,  exec = 1;
    i = pt_get_index(arg);
    pt_set_activation(i);
    set_com_variable(0);
    value = ql_get_learning_rate();

    while(exec)
    {
        key_manager(&exec);
        pt_deadline_miss(i);
        pt_wait_for_period(i);
        if(exec==0)
        {
            printf("ho finito \n");
        }
    }
    printf("INTERPRETER: task end\n");
    return NULL;
}

void * wave_gener(void *arg)
{
    printf("wave generator task started\n");
    int i,exec=1,ex_stat,j;
    static float counter =-3.14;
    float F_S[6];
    i = pt_get_index(arg);
    pt_set_activation(i);
    

    while(exec)
    {
        get_com_variable(&ex_stat);
        if(ex_stat==1)
        {
            counter=counter+PASSO;
            for(j=0;j<6;j++)
                F_S[j]=cos(counter+PI*(j/10));
            set_FALSE_ST(F_S);

        }
        else if(ex_stat==3)
            exec=0;
        
        pt_deadline_miss(i);
        pt_wait_for_period(i);
    }
    printf("wave generator task end\n");
    return NULL;
}