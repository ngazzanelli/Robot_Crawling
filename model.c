#include <math.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include "ptask.h"
#include "matrices.h"

#define T   0.1       // Intervallo di generazione della traiettoria      [s]

// COSTANTI DEL CONTROLLO
#define KC  10 
#define KD  10

// COSTANTI PER LO STATO DEL SISTEMA
#define RESET   0
#define PLAY    1
#define PAUSE   2
#define STOP    3

//COSTANTI PER IL CALCOLO DELLA CINEMATICA
#define TRUE_ALPHA  1
#define FALSE_ALPHA 0

//COSTANTI DEL ROBOT
#define ROBOT_LENGHT    10.0
#define WHEEL_RADIUS    1.5

static float DT = 0.001;          // Intervallo di integrazione della dinamica        [s]
static pthread_mutex_t mux_dt = PTHREAD_MUTEX_INITIALIZER;
void set_dyn_dt(float dt){
    pthread_mutex_lock(&mux_dt);
    DT = dt;
    pthread_mutex_unlock(&mux_dt);
}
float get_dyn_dt(){
    float ret;
    pthread_mutex_lock(&mux_dt);
    ret = DT;
    pthread_mutex_unlock(&mux_dt);
    return ret;
}

//Struttura per lo stato desiderato del robot
typedef struct {
    float q1d;
    float q2d;
    int flag;
} target;
static target qd;
static int model_dl;

//Funzione dal crawler per l'aggiornamento dello stato desiderato
extern void get_desired_joint(target* t);

// Strutture per lo stato reale del robot: dichiarate static perchè possono essere accedute solo tramite
// funzioni di interfaccia
static state robot_globale;       
static dot_state dot_robot;

//Semaforo usato per accedere allo stato reale "robot"
static pthread_mutex_t mux = PTHREAD_MUTEX_INITIALIZER;
//semaforo usato per accedre alle deadline miss
static pthread_mutex_t mux_model_dl = PTHREAD_MUTEX_INITIALIZER;

void inc_model_dl()
{
    pthread_mutex_lock(&mux_model_dl);
    model_dl++;
    pthread_mutex_unlock(&mux_model_dl);
}

//funzione che ottiene il valore della deadline miss
void get_model_dl(int * dl_miss)
{
    pthread_mutex_lock(&mux_model_dl);
    * dl_miss = model_dl;
    pthread_mutex_unlock(&mux_model_dl);
}

void init_state(){

    pthread_mutex_lock(&mux);
    robot_globale.q1 = 0;
    robot_globale.q2 = 0;
    robot_globale.q3 = 0;
    robot_globale.q4 = 0;
    robot_globale.q5 = 0;
    robot_globale.q6 = 0;
    robot_globale.energy = 0;
    robot_globale.dt3 = 0;
    pthread_mutex_unlock(&mux);

    dot_robot.dq1 = 0;
    dot_robot.dq2 = 0;
    dot_robot.dq3 = 0;
    dot_robot.dq4 = 0;
    dot_robot.dq5 = 0;
    dot_robot.dq6 = 0;
}

void get_state(state* s){
    pthread_mutex_lock(&mux);
    s->q1 = robot_globale.q1;
    s->q2 = robot_globale.q2;
    s->q3 = robot_globale.q3;
    s->q4 = robot_globale.q4;
    s->q5 = robot_globale.q5;
    s->q6 = robot_globale.q6;
    s->energy = robot_globale.energy;
    s->dt3 = robot_globale.dt3;
    pthread_mutex_unlock(&mux);
}

void set_state(state s){
    pthread_mutex_lock(&mux);
    robot_globale.q1 = s.q1;
    robot_globale.q2 = s.q2;
    robot_globale.q3 = s.q3;
    robot_globale.q4 = s.q4;
    robot_globale.q5 = s.q5;
    robot_globale.q6 = s.q6 ;
    robot_globale.energy = s.energy;
    robot_globale.dt3 = s.dt3;
    pthread_mutex_unlock(&mux);
}

//Funzione dall'interprete
extern int get_sys_state(int *s);

void update_coefficients(float coef1[4], float coef2[4], state robot){

    //Prelevo lo stato attuale in modo da generare la traiettoria a partire da dove si trova il robot
    float qi_1 = robot.q4;
    float qi_2 = robot.q5;

    //Prelevo la posizione finale desiderata dalla variabile globale qd che è stata aggiornata in generate_tau
    float qf_1 = qd.q1d;
    float qf_2 = qd.q2d;

    //Update coefficienti per il primo giunto
    coef1[0] = qi_1;
    coef1[1] =  10*(qf_1 - qi_1)/pow(T, 3);
    coef1[2] = -15*(qf_1 - qi_1)/pow(T, 4);
    coef1[3] =   6*(qf_1 - qi_1)/pow(T, 5);

    //Update coefficienti per il secondo giunto
    coef2[0] = qi_2;
    coef2[1] =  10*(qf_2 - qi_2)/pow(T, 3);
    coef2[2] = -15*(qf_2 - qi_2)/pow(T, 4);
    coef2[3] =   6*(qf_2 - qi_2)/pow(T, 5);

    return;

}

void compute_qdt(float qdt[2], float dot_qdt[2], float dotdot_qdt[2], float coef1[4], float coef2[4], float t){
    if(t <= T){
        //Posizioni di giunto desiderate
        qdt[0] = coef1[0] + coef1[1]*pow(t, 3) + coef1[2]*pow(t, 4) + coef1[3]*pow(t, 5);
        qdt[1] = coef2[0] + coef2[1]*pow(t, 3) + coef2[2]*pow(t, 4) + coef2[3]*pow(t, 5);
        //Velocità di giunto desiderate
        dot_qdt[0] = 3*coef1[1]*pow(t, 2) + 4*coef1[2]*pow(t, 3) + 5*coef1[3]*pow(t, 4);
        dot_qdt[1] = 3*coef2[1]*pow(t, 2) + 4*coef2[2]*pow(t, 3) + 5*coef2[3]*pow(t, 4);
        //Accelerazioni di giunto desiderate
        dotdot_qdt[0] = 6*coef1[1]*t + 12*coef1[2]*pow(t, 2) + 20*coef1[3]*pow(t, 3);
        dotdot_qdt[1] = 6*coef2[1]*t + 12*coef2[2]*pow(t, 2) + 20*coef2[3]*pow(t, 3);
    }
    else{
        qdt[0] = qd.q1d;
        qdt[1] = qd.q2d;
        dot_qdt[0] = dot_qdt[1] = 0;
        dotdot_qdt[0] = dotdot_qdt[1] = 0;
    }

    return;
}

void generate_tau(float tau[2], state robot, float M[2][2], float C[2][2], float G[2]){
    static int step;
    static float coefficients1[4];
    static float coefficients2[4];

    float e[2];
    float dot_e[2];
    float s[2];
    float dot_qr_t[2];      //Velocità di giunto di riferimento all'istante t
    float dotdot_qr_t[2];   //Accelerazioni di giunto di riferimento all'istannte t
    float q_t[2];           //Variabili di giunto all'istante t
    float dot_q_t[2];       //Velocità di giunto all'istante t
    float dotdot_q_t[2];    //Accelerazioni di giunto all'istante t
    float qd_t[2];          //Variabili di giunto desiderate all'istante t
    float dot_qd_t[2];      //Velocità di giunto desiderate all'istante t
    float dotdot_qd_t[2];   //Accelerazioni di giunto desiderate all'istante t
    float t;                //Variabile per l'istante temporale
    float ris1[2], ris2[2], ris3[2], ris4[2], ris5[2], ris6[2], ris7[2]; // Variabili di appoggio per le operazioni

    get_desired_joint(&qd);

    /*******DEBUG****
    if(step == 0){
        qd.flag = 1;
        qd.q1d  = M_PI/5;
        qd.q2d = 0;
    }
    else
        qd.flag = 0;
    ****************/

    if(qd.flag == 1){
        step = 0;
        update_coefficients(coefficients1, coefficients2, robot);
    }

    t = step*get_dyn_dt();
    q_t[0] = robot.q4;
    q_t[1] = robot.q5;
    dot_q_t[0] = dot_robot.dq4;
    dot_q_t[1] = dot_robot.dq5;
    compute_qdt(qd_t, dot_qd_t, dotdot_qd_t, coefficients1, coefficients2, t);

    // e = q_d - q;
    vector_sub(qd_t, q_t, e, 2);

    // q_r_dot = q_d_dot + k_c * e
    vector_scal(e, KC, ris1, 2);
    vector_sum(dot_qd_t, ris1, dot_qr_t, 2);

    // s = q_r_dot - q_dot
    vector_sub(dot_qr_t, dot_q_t, s, 2);

    // e_dot = qd_dot - q_dot
    vector_sub(dot_qd_t, dot_q_t, dot_e, 2);

    // qr_dotdot = qd_dotdot + k_c*e_dot
    vector_scal(dot_e, KC, ris1, 2);
    vector_sum(dotdot_qd_t, ris1, dotdot_qr_t, 2);

    // tau = M * q_r_dotdot + C * q_r_dot + G + e + k_d * s
    vector_scal(s, KD, ris1, 2);
    vector_sum(e, ris1, ris2, 2);
    vector_sum(G, ris2, ris3, 2);
    matvec_mul(dot_qr_t, ris5, 2, 2, C);
    vector_sum(ris5, ris3, ris6, 2);
    matvec_mul(dotdot_qr_t, ris7, 2, 2, M);
    vector_sum(ris7, ris6, tau, 2);

    step++;
    return;
}

float adjuste_alpha(float x_p, float y_p){
    float theta1, theta2;
    y_p -= WHEEL_RADIUS;
    x_p += ROBOT_LENGHT/2; 
    theta1 = atan2(y_p, x_p);
    theta2 = atan2(WHEEL_RADIUS, x_p);
    return (theta1 + theta2);
}

void* dynamics(void* arg){

    printf("DYNAMIC: task started\n");
    init_state();
    int i,exec;            // thread index and system state
    float y_ee;
    float dt;
    float theta;        
    state robot;
    get_state(&robot);
    //Vettori per lo stato a un passo e al successivo
    float q_ind1[2], q_ind2[2], q_dip1[4], q_dip2[4], ris1[2], ris2[2], ris3[2], ris4[4];
    float qdotdot_ind[2];
    float qdot_ind1[2], qdot_ind2[2], qdot_dip1[4], qdot_dip2[4];
    //inizializziamo a 0 questi vettori
    vector_set_zero(q_ind1, 2);
    vector_set_zero(q_dip1, 4);
    vector_set_zero(qdot_ind1, 2);
    vector_set_zero(qdot_dip1, 4);

    //Matrici per la dinamica
    float M[2][2], C[2][2], M_inv[2][2], Tsee[4][4], S[4][2];
    float G[2], tau[2];

    i = pt_get_index(arg);
    pt_set_activation(i);
    //printf("DYN: il mio periodo è %d microsecondi\n", pt_get_period(i));
    
    while(get_sys_state(&exec) != STOP){

        //controllo se l'applicazione è in pausa o in reset
        if(exec==PLAY){
            //printf("DYNAMIC: il valore di q è: [%f %f %f %f %f %f]\n",robot.q1, robot.q2, robot.q3, robot.q4, robot.q5, robot.q6);

            dt = get_dyn_dt();
            //printf("DYN: dt vale %f\n", dt);
            update_kyn(Tsee, robot, FALSE_ALPHA);
            y_ee = Tsee[1][3];
            //printf("Y = %f\n", y_ee);
            if(y_ee > 0){
                //printf("DYN: Sono dentro y > 0 con y= %f\n",y_ee);
                update_M1(M, robot);
                update_C1(C, robot, dot_robot);
                update_G1(G, robot);
            }else{
                //printf("DYN: Sono dentro y < 0 con y= %f\n",y_ee);
                update_kyn(Tsee, robot, TRUE_ALPHA);
                y_ee = Tsee[1][3];
                if(y_ee > 0){
                    printf("DYN: Sono dentro y < 0 ma il secondo check dice y=%f\n", y_ee);
                    theta = adjuste_alpha(Tsee[0][3], Tsee[1][3]); 
                    robot.q3 -= theta;  
                    q_dip1[2] = robot.q3;
                }
                update_M2(M, robot);
                update_C2(C, robot, dot_robot);
                update_G2(G, robot);
            }

            generate_tau(tau, robot, M, C, G);
            //printf("G = [%f; %f]\n", G[0], G[1]);
            //printf("TAU = [%f; %f]\n", tau[0], tau[1]);
            matrix_inverse(M, M_inv);
            
            //qdotdot_ind = M_inv*(tau-G-C*qdot_ind1);    //accelerazione attuale variabili indipendenti 
            vector_sub(tau, G, ris1, 2);
            matvec_mul(qdot_ind1, ris2, 2, 2, C);
            vector_sub(ris1, ris2, ris3, 2);
            matvec_mul(ris3, qdotdot_ind, 2, 2, M_inv);    

            //qdot_ind2 = qdot_ind1 + dt*qdotdot_ind;     //velocità attuale variabili indipendenti 
            vector_scal(qdotdot_ind, dt, ris1, 2);
            vector_sum(qdot_ind1, ris1, qdot_ind2, 2);

            //q_ind2 = q_ind1 + dt*qdot_ind1;             //posizione attuale variabili indipendenti 
            vector_scal(qdot_ind1, dt, ris1, 2);
            vector_sum(q_ind1, ris1, q_ind2, 2);

            //q_dip2 = q_dip1 + dt*qdot_dip1;             //posizione attuale variabili dipendenti 
            vector_scal(qdot_dip1, dt, ris4, 4);
            vector_sum(q_dip1, ris4, q_dip2, 4);

            if(y_ee > 0)
                matrix_set_zero(4, 2, S);
            else{
                /* TODO: aggiornare lo stato prima di calcolare S2 */
                robot.q3 = q_dip2[2];
                robot.q4 = q_ind2[0];
                robot.q5 = q_ind2[1];
                update_S2(S, robot);
            }
                    
            //qdot_dip2 = S*qdot_ind2;                    //velocità attuale variabili dipendenti
            matvec_mul(qdot_ind2, qdot_dip2, 4, 2, S);

            //Aggiorniamo lo stato al passo precedente 
            vector_copy(q_dip2, q_dip1, 4);
            vector_copy(q_ind2, q_ind1, 2);
            vector_copy(qdot_dip2, qdot_dip1, 4);
            vector_copy(qdot_ind2, qdot_ind1, 2);


            //Aggiorniamo lo stato locale utilizzato soltanto per l'integrazione
            //robot.dt3 = q_dip2[0] - robot.q1;
            robot.q1 = q_dip2[0];
            robot.q2 = q_dip2[1];
            robot.q3 = q_dip2[2];

            /* CHECK SU ALPHA e Y */
            if(robot.q3 < 0){
                q_dip1[2] = 0;
                robot.q3 = 0;
            }

            if(robot.q2 != 0){
                q_dip1[1] = 0;
                robot.q2 = 0;   
            }
            /*********************/

            robot.q4 = q_ind2[0];
            robot.q5 = q_ind2[1];
            robot.q6 = q_dip2[3];
            //robot.energy++;


            dot_robot.dq1 = qdot_dip2[0]; 
            dot_robot.dq2 = qdot_dip2[1];
            dot_robot.dq3 = qdot_dip2[2];
            dot_robot.dq4 = qdot_ind2[0];
            dot_robot.dq5 = qdot_ind2[1];
            dot_robot.dq6 = qdot_dip2[3];

            //Aggiorniamo lo stato globale condiviso con gli altri task
            set_state(robot);
            //printf("DYN: state variables q = [%f; %f; %f; %f; %f; %f]\n", robot.q1, robot.q2, robot.q3, robot.q4, robot.q5, robot.q6);
        }

        //Riprendiamo il valore dello stato globale se siamo in reset
        if(exec==RESET){
            get_state(&robot);
            vector_set_zero(q_ind1, 2);
            vector_set_zero(q_dip1, 4);
            vector_set_zero(qdot_ind1, 2);
            vector_set_zero(qdot_dip1, 4);
        }

        if(pt_deadline_miss(i)){
            inc_model_dl();
            /*printf("DYN: ho missato una deadline\n");
            printf("DYN: il mio periodo vale %d microsecondi\n", pt_get_period(i));
            printf("DYN: la mia deadline relativa vale %d microsecondi\n", pt_get_deadline(i));*/
        }
        
        pt_wait_for_period(i);
    }

    printf("DYNAMIC: task finished\n");
    return NULL;
}
