#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include "ptask.h"
#include <pthread.h>
#include "matrices.h"

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

extern void get_desired_joint(target* t);

// Struttura per lo stato reale del robot
state robot;

dot_state dot_robot;


static pthread_mutex_t mux = PTHREAD_MUTEX_INITIALIZER;


void init_state(){

    pthread_mutex_lock(&mux);
    robot.q1 = 0;
    robot.q2 = 0;
    robot.q3 = 0;
    robot.q4 = 0;
    robot.q5 = 0;
    robot.q6 = 0;
    robot.energy = 0;
    robot.dt3 = 0;
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
    s->q1 = robot.q1;
    s->q2 = robot.q2;
    s->q3 = robot.q3;
    s->q4 = robot.q4;
    s->q5 = robot.q5;
    s->q6 = robot.q6;
    s->energy = robot.energy;
    s->dt3 = robot.dt3;
    pthread_mutex_unlock(&mux);
}

extern int get_stop();
/*
void generate_tau(gsl_vector *tau){

}

gsl_matrix* compute_inverse(gsl_matrix *m){

    gsl_permutation *p;
    gsl_matrix *inv;
    int s;
    int size;

    size = m->size1;
    p = gsl_permutation_alloc(size);
    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(m, p, &s);

    // Compute the  inverse of the LU decomposition
    inv = gsl_matrix_alloc(size, size);
    gsl_linalg_LU_invert(m, p, inv);

    gsl_permutation_free(p);

    return inv;

}
*/
void* dynamics(void* arg){
    printf("dynamic task started\n");
    int i;      // thread index
    float y_ee;
    float dt = 0.001; // 1 ms
    //Vettori per lo stato a un passo e al successivo
    float q_ind1[2], q_ind2[2], q_dip1[4], q_dip2[4], ris1[2], ris2[2], ris3[2], ris4[4];
    float qdotdot_ind[2];
    float qdot_ind1[2], qdot_ind2[2], qdot_dip1[4], qdot_dip2[4];

    //Matrici per la dinamica
    float M[2][2], C[2][2], M_inv[2][2], Tsee[4][4], S[4][2];
    float G[2], tau[2];

    i = pt_get_index(arg);
    pt_set_activation(i);

    init_state();

    while(/*!get_stop()*/1){

        update_kyn(Tsee, robot);
        //y_ee = gsl_matrix_get(Tsee, 1, 3);
        y_ee = Tsee[1][3];

        if(y_ee <=0){
            update_M1(M, robot);
            update_C1(C, robot, dot_robot);
            update_G1(G, robot);
        }else{
            update_M2(M, robot);
            update_C2(C, robot, dot_robot);
            update_G2(G, robot);
        }

        //generate_tau(tau);
        //M_inv = compute_inverse(M);
        
        //qdotdot_ind = M_inv*(tau-G-C*qdot_ind1);    //accelerazione attuale variabili indipendenti 
        sub(tau, G, ris1, 2);
        mul(qdot_ind1, ris2, 2, 2, C);
        sum(ris1, ris2, ris3, 2);
        mul(ris3, qdotdot_ind, 2, 2, M_inv);    

        //qdot_ind2 = qdot_ind1 + dt*qdotdot_ind;     //velocità attuale variabili indipendenti 
        scal(qdotdot_ind, dt, ris1, 2);
        sum(qdot_ind1, ris1, qdot_ind2, 2);

        //q_ind2 = q_ind1 + dt*qdot_ind1;             //posizione attuale variabili indipendenti 
        scal(qdot_ind1, dt, ris1, 2);
        sum(q_ind1, ris1, q_ind2, 2);

        //q_dip2 = q_dip1 + dt*qdot_dip1;             //posizione attuale variabili dipendenti 
        scal(qdot_dip1, dt, ris4, 4);
        sum(q_dip1, ris1, q_dip2, 4);

        if(y_ee <=0)
            matrix_set_zero(4, 2, S);
        else
            update_S2(S, robot);

        //qdot_dip2 = S*qdot_ind2;                    //velocità attuale variabili dipendenti
        mul(qdot_ind2, qdot_dip2, 4, 2, S);

        pt_deadline_miss(i);
        pt_wait_for_period(i);
    }

    printf("dynamic task finished\n");
    return NULL;
}

/*int main(){
    gsl_matrix *M, *M_inv;
    M = gsl_matrix_alloc(3, 3);
    M_inv = gsl_matrix_alloc(3, 3);

    gsl_matrix_set(M, 0, 0, 5);
    gsl_matrix_set(M, 1, 1, 5);
    gsl_matrix_set(M, 2, 2, 5);

    M_inv = compute_inverse(M);
    print_matrix(M_inv);
}*/