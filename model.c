#include <math.h>
#include <stdio.h>
#include <pthread.h>
#include "ptask.h"
#include "matrices.h"

#define DT  10   //[ms]
// Queste probabilmente non serviranno se si suppone che i parametri non possano essere cambiati 
// (sono impliciti nella risoluzione della dinamica)
#define L1      0.06     // link1 length [m]
#define L2      0.06     // link2 length [m]
#define RH      0.03     // robot height [m]
#define RL      0.2     // robot length [m]

//Struttura per lo stato desiderato del robot
typedef struct {
    float q1d;
    float q2d;
    int flag;
} target;
static target qd;
//Funzione dal crawler per l'aggiornamento dello stato desiderato
extern void get_desired_joint(target* t);

// Strutture per lo stato reale del robot
state robot_globale;       
dot_state dot_robot;
//Semaforo usato per accedere allo stato reale "robot"
static pthread_mutex_t mux = PTHREAD_MUTEX_INITIALIZER;

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
extern int get_stop();
extern int get_pause();

void update_coefficients(float coef1[4], float coef2[4], target qd){
    
}

void compute_qdt(float qdt[2], float coef1[4], float coef2[4], float t){

}

void generate_tau(float tau[2], state robot){
    static int step;
    static float coefficients1[4];
    static float coefficients2[4];
    float qd_t[2];  //variabili di giunto desiderate all'istante t
    float t;        //variabile per l'istante temporale
    get_desired_joint(&qd);
    if(qd.flag == 1){
        step = 0;
        update_coefficients(coefficients1, coefficients2, qd, robot);
    }
    t = step*DT*0.001;
    compute_qdt(qd_t, coefficients1, coefficients2, t);

    step++;
}

//La funzione genera la traiettoria da inseguire e restituisce in res:
/*
    res[0] => variabile di giunto desiderata
    res[1] => velocità di giunto desiderata
    res[2] => accelerazione di giunto desiderata
*/
void generate_trajectory(float t, int joint, int direction, float res[3]){

    static float T = 1;
    static float old_q1 = 0;
    static float old_q2 = 0;

    float qi, qf;
    float a0, a3, a4, a5;

    if(joint == 0){
        qi = old_q1;
        qf = direction ? (old_q1 + M_PI/9) : (old_q1 - M_PI/9);
    } else {
        qi = old_q2;
        qf = direction ? (old_q2 + M_PI/9) : (old_q2 - M_PI/9);
    }

    a0 = qi;
    a3 =  10*(qf - qi)/pow(T, 3);
    a4 = -15*(qf - qi)/pow(T, 4);
    a5 =   6*(qf - qi)/pow(T, 5);

    if(joint == 0)
        old_q1 = qf;
    else 
        old_q2 = qf;

    if(t < T){
        res[0] = a0 + a3*pow(t, 3) + a4*pow(t, 4) + a5*pow(t, 5);
        res[1] = 3*a3*pow(t, 2) + 4*a4*pow(t, 3) + 5*a5*pow(t, 4);
        res[2] = 6*a3*t + 12*a4*pow(t, 2) + 20*a5*pow(t, 3);
    }
    else{
        res[0] = qf;
        res[1] = res[2] = 0;
    }
    
/*
    if(joint == 0){
        dq1 = qf;
        if(t < T){
            tq1 = a0 + a3*pow(t, 3) + a4*pow(t, 4) + a5*pow(t, 5);
            dtq1 = 3*a3*pow(t, 2) + 4*a4*pow(t, 3) + 5*a5*pow(t, 4);
            ddtq1 = 6*a3*t + 12*a4*pow(t, 2) + 20*a5*pow(t, 3);
        }
        else{
            tq1 = qf;
            dtq1 = ddtq1 = 0;
        }
        res[0] = tq1;
        res[1] = dtq1;
        res[2] = ddtq1;
        tq2 = dq2;

    }else{
        dq2 = qf;
    if(joint == 0)
        dq1 = qf;
    else 
        dq2 = qf;

    if(t < T){
        res[0] = a0 + a3*pow(t, 3) + a4*pow(t, 4) + a5*pow(t, 5);
        res[1] = 3*a3*pow(t, 2) + 4*a4*pow(t, 3) + 5*a5*pow(t, 4);
        res[2] = 6*a3*t + 12*a4*pow(t, 2) + 20*a5*pow(t, 3);
    }
    else{
        res[0] = qf;
        res[1] = res[2] = 0;
    }
                tq2 =  a0 + a3*pow(t, 3) + a4*pow(t, 4) + a5*pow(t, 5);
            dtq2 = 3*a3*pow(t, 2) + 4*a4*pow(t, 3) + 5*a5*pow(t, 4);
            ddtq2 = 6*a3*t + 12*a4*pow(t, 2) + 20*a5*pow(t, 3);
        }
        else{
            tq2 = qf;
            dtq2 = ddtq2 = 0;
        }
        res[0] = tq2;
        res[1] = dtq2;
        res[2] = ddtq2;
        tq1 = dq1;
    }
*/

}

void* dynamics(void* arg){

    printf("dynamic task started\n");
    int i;            // thread index
    float y_ee;
    float dt = 0.001; // 1 ms
    state robot;
    get_state(&robot);
    //Vettori per lo stato a un passo e al successivo
    float q_ind1[2], q_ind2[2], q_dip1[4], q_dip2[4], ris1[2], ris2[2], ris3[2], ris4[4];
    float qdotdot_ind[2];
    float qdot_ind1[2], qdot_ind2[2], qdot_dip1[4], qdot_dip2[4];

    //Matrici per la dinamica
    float M[2][2], C[2][2], M_inv[2][2], Tsee[4][4], S[4][2];
    float G[2], tau[2];

    i = pt_get_index(arg);
    pt_set_activation(i);

    while(!get_stop()){

        //controllo se l'applicazione è in pausa
        if(!get_pause()){

            update_kyn(Tsee, robot);
            y_ee = Tsee[1][3];
            if(y_ee >0){
                update_M1(M, robot);
                update_C1(C, robot, dot_robot);
                update_G1(G, robot);
            }else{
                update_M2(M, robot);
                update_C2(C, robot, dot_robot);
                update_G2(G, robot);
            }

            generate_tau(tau);
            matrix_inverse(M, M_inv);
            
            //qdotdot_ind = M_inv*(tau-G-C*qdot_ind1);    //accelerazione attuale variabili indipendenti 
            vector_sub(tau, G, ris1, 2);
            matvec_mul(qdot_ind1, ris2, 2, 2, C);
            vector_sum(ris1, ris2, ris3, 2);
            matvec_mul(ris3, qdotdot_ind, 2, 2, M_inv);    

            //qdot_ind2 = qdot_ind1 + dt*qdotdot_ind;     //velocità attuale variabili indipendenti 
            vector_scal(qdotdot_ind, dt, ris1, 2);
            vector_sum(qdot_ind1, ris1, qdot_ind2, 2);

            //q_ind2 = q_ind1 + dt*qdot_ind1;             //posizione attuale variabili indipendenti 
            vector_scal(qdot_ind1, dt, ris1, 2);
            vector_sum(q_ind1, ris1, q_ind2, 2);

            //q_dip2 = q_dip1 + dt*qdot_dip1;             //posizione attuale variabili dipendenti 
            vector_scal(qdot_dip1, dt, ris4, 4);
            vector_sum(q_dip1, ris1, q_dip2, 4);

            if(y_ee >0)
                matrix_set_zero(4, 2, S);
            else
                update_S2(S, robot);

            //qdot_dip2 = S*qdot_ind2;                    //velocità attuale variabili dipendenti
            matvec_mul(qdot_ind2, qdot_dip2, 4, 2, S);

            //Aggiorniamo lo stato globale condiviso con gli altri task
            set_state(robot);
        }
        
        pt_deadline_miss(i);
        pt_wait_for_period(i);
    }

    printf("dynamic task finished\n");
    return NULL;
}
