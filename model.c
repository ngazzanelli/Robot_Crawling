#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include "ptask.h"
#include <pthread.h>

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
state robot;

typedef struct {
    float dq1;
    float dq2;
    float dq3;
    float dq4;
    float dq5;
    float dq6;
} dot_state;
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

extern void update_kyn(gsl_matrix *Tsee, state robot);
extern void update_M1(gsl_matrix *M1, state robot);
extern void update_M2(gsl_matrix *M2, state robot);
extern void update_C1(gsl_matrix *C1, state robot, dot_state dot_robot);
extern void update_C2(gsl_matrix *C2, state robot, dot_state dot_robot);
extern void update_G1(gsl_vector *G1, state robot);
extern void update_G2(gsl_vector *G2, state robot);
extern void update_S2(gsl_matrix *S, state robot);
extern void print_matrix(gsl_matrix *m);
extern int get_stop();

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

void* dynamics(void* arg){
    printf("dynamic task started\n");
    int i;      // thread index
    float y_ee;
    float dt = 0.001; // 1 ms
    //Vettori per lo stato a un passo e al successivo
    gsl_vector *q_ind1, *q_ind2, *q_dip1, *q_dip2;
    gsl_vector *qdotdot_ind;
    gsl_vector *qdot_ind1, *qdot_ind2, *qdot_dip1, *qdot_dip2;

    //Matrici per la dinamica
    gsl_matrix *M, *C, *M_inv, *Tsee, *S;
    gsl_vector *G, *tau;

    //Allocazione vettori e matrici
    q_ind1 = gsl_vector_alloc(2); 
    q_ind2 = gsl_vector_alloc(2); 
    qdot_ind1 = gsl_vector_alloc(2); 
    qdot_ind2 = gsl_vector_alloc(2); 
    q_dip1 = gsl_vector_alloc(4); 
    q_dip2 = gsl_vector_alloc(4); 
    qdot_dip1 = gsl_vector_alloc(4); 
    qdot_dip2 = gsl_vector_alloc(4);
    qdotdot_ind = gsl_vector_alloc(2);

    M = gsl_matrix_alloc(2, 2);
    M_inv = gsl_matrix_alloc(2, 2);
    C = gsl_matrix_alloc(2, 2);
    G = gsl_vector_alloc(2);
    S = gsl_matrix_alloc(4, 2);
    tau = gsl_vector_alloc(2);

    robot.q1 = gsl_vector_get(q_ind2, 0);

    i = pt_get_index(arg);
    pt_set_activation(i);

    Tsee = gsl_matrix_alloc(4, 4);
    init_state();

    while(/*!get_stop()*/1){

        update_kyn(Tsee, robot);
        y_ee = gsl_matrix_get(Tsee, 1, 3);

        if(y_ee <=0){
            update_M1(M, robot);
            update_C1(C, robot, dot_robot);
            update_G1(G, robot);
            gsl_matrix_set_zero(S);

        }else{
            update_M2(M, robot);
            update_C2(C, robot, dot_robot);
            update_G2(G, robot);
            update_S2(S, robot);
        }

        generate_tau(tau);
        M_inv = compute_inverse(M);
        //qdotdot_ind = prod(M_inv, sum(diff(tau, G), prod(C,qdot_ind1));    //accelerazione attuale variabili indipendenti 
        qdotdot_ind = M_inv*(tau-G-C*qdot_ind1);    //accelerazione attuale variabili indipendenti 
        qdot_ind2 = qdot_ind1 + dt*qdotdot_ind;     //velocità attuale variabili indipendenti 
        q_ind2 = q_ind1 + dt*qdot_ind1;             //posizione attuale variabili indipendenti 
        q_dip2 = q_dip1 + dt*ddot_dip1;             //posizione attuale variabili dipendenti 
        if(y_ee <=0)
            gsl_matrix_set_zero(S);
        else
            update_S2(S, robot);
        qdot_dip2 = S*qdot_ind2;                    //velocità attuale variabili dipendenti



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