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
extern void update_C1(gsl_matrix *C1, state robot);
extern void update_C2(gsl_matrix *C2, state robot);
extern void update_G1(gsl_vector *G1, state robot);
extern void update_G2(gsl_vector *G2, state robot);
extern void print_matrix(gsl_matrix *m);
extern int get_stop();

void generate_tau(gsl_vecctor *tau){

}

gsl_matrix* compute_inverse(gsl_matrix *m){

    gsl_permutation *p;
    gsl_matrix *inv;
    int s;
    int size;

    size = m->size1;
    gsl_permutation *p = gsl_permutation_alloc(size);
    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    inv = gsl_matrix_alloc(size, size);
    gsl_linalg_LU_invert(matrix, p, inv);

    gsl_permutation_free(p);

    return inv;

}

void* dynamics(void* arg){

    int i;      // thread index
    gsl_matrix *Tsee;
    float y_ee;
    //Matrici per la dinamica
    gsl_matrix *M, *C, *M_inv;
    gsl_vector *G, *tau;
    M = gsl_matrix_alloc(2, 2);
    M_inv = gsl_matrix_alloc(2, 2);
    C = gsl_matrix_alloc(2, 2);
    G = gsl_vector_alloc(2);
    tau = gsl_vector_alloc(2);


    i = pt_get_index(arg);
    pt_set_activation(i);

    Tsee = gsl_matrix_alloc(4, 4);
    init_state();

    while(!get_stop()){

        update_kyn(Tsee, robot);
        y_ee = gsl_matrx_get(Tsee, 1, 3);
        if(y_ee <=0){
            update_M1(M, robot);
            update_C1(C, robot, dot_robot);
            update_G1(G, robot, dot_robot);
            gsl_matrix_set_zero(S);

        }else{
            update_M2(M, robot);
            update_C2(C, robot, dot_robot);
            update_G2(G, robot, dot_robot);
            update_S2(S, robot, dot_robot);
        }

        generate_tau(tau);
        compute_inverse(M, M_inv);

        pt_deadline_miss(i);
        pt_wait_for_period(i);
    }
    return NULL;
}