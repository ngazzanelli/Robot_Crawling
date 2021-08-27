#include <allegro.h>
#include <math.h>
#include <stdio.h>
#include "ptask.h"

#include <string.h>
#include <pthread.h>

#define PI 3.14
#define PASSO 0.01

static pthread_mutex_t mux_int = PTHREAD_MUTEX_INITIALIZER;
//per provare sia l'interprete che la capacità di scrivere 
//i igusti plot mi sono definito uno stato "falso" in cui va 
//a lavorare il generatore di onde; nella grafica basta inserire
//la dichiarazione esterna della get_state reale, inserire 
//nel corpo del thread la dichiarazione della struttura e 
//dovrebbe fungere lo stesso
static pthread_mutex_t mux_false_state = PTHREAD_MUTEX_INITIALIZER;

float FALSE_STATE[6];

int com_state;
// variabile globale  che controlla il flusso delle operazioni :
   //case 0: sistema appena acceso, unico modo per settare
   //le variabili del qlearning
   //case 1: sistema in play, eseguetutte le istruzioni 
   //case 2: sistema in pausa, i thread rimangono attivi 
   //ma non eseguono il flusso originale
   //case 3: arresto del sistema 

int get_pause(){
    return 0;
}
int get_stop(){
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
void set_com_variable(int i)
{
    if(i>=0 && i<4)
    {
        pthread_mutex_lock(&mux_int);
        com_state=i;
        pthread_mutex_unlock(&mux_int);
    }
}
void get_com_variable(int *ic)
{
    pthread_mutex_lock(&mux_int);
    *ic = com_state;
    pthread_mutex_unlock(&mux_int);
}
//funzione per ottenere il codice del tasto premuto da tastiera 
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

void key_mamager(int * exec)
{
    char cm;
    cm=get_scancode();
    switch(cm)
    {
        case KEY_S:
        printf("hai premuto il tasto S\n");
        set_com_variable(1);
        break;
        case KEY_P:
        printf("hai premuto il tasto P\n");
        set_com_variable(2);
        break;
        case KEY_E:
        printf("hai premuto il tasto E\n");
        set_com_variable(3);  
        *exec=0;
        break;   
    }
}
void* interface(void * arg)
{
    printf("interpreter task started\n");
    int i,exec=1;
    i = pt_get_index(arg);
    pt_set_activation(i);
    

    while(exec)
    {
        key_mamager(&exec);
        pt_deadline_miss(i);
        pt_wait_for_period(i);
        if(exec==0)
        {
            printf("ho finito \n");
        }
    }
    printf("interpreter task end\n");
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