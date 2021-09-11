#include <math.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include "ptask.h"
#include "matrices.h"


#define ROBOT_LENGHT    10.0
#define WHEEL_RADIUS    1.5
#define T               0.1  // Trajectory Generation Interval [s]

// Control Constants
#define KC  10 
#define KD  10

// System State Constants
#define RESET   0
#define PLAY    1
#define PAUSE   2
#define STOP    3

// Kynematic Computation Constants
#define TRUE_ALPHA  1
#define FALSE_ALPHA 0


// Desired Robot State (from crawler task)
typedef struct {
    float q1d;
    float q2d;
    int flag;
} target;


// Functions from other modules
extern int get_sys_state(int *s);           //(from command_interface task)
extern void get_desired_joint(target* t);   //(from crawler task)

// Mutexes
static pthread_mutex_t mux_dt = PTHREAD_MUTEX_INITIALIZER;          // Mutex for Dynamic Integration Interval dt
static pthread_mutex_t mux = PTHREAD_MUTEX_INITIALIZER;             // Mutex for Actual Robot State

// Static Variables
static float    delta_t = 0.001;    // Dynamic Integration Interval dt [s]
static target   qd;

// Actual Robot State: can be accessed only through Interface Function
static state global_robot;       
static dot_state dot_robot;

//-----------------------------------------------------
// The following Functions manage Robot State for 
// initialization and interfacing with other tasks
//-----------------------------------------------------
void init_state()
{
    pthread_mutex_lock(&mux);
    global_robot.q1 = 0;
    global_robot.q2 = 0;
    global_robot.q3 = 0;
    global_robot.q4 = 0;
    global_robot.q5 = 0;
    global_robot.q6 = 0;
    global_robot.energy = 0;
    global_robot.dt3 = 0;
    pthread_mutex_unlock(&mux);

    dot_robot.dq1 = 0;
    dot_robot.dq2 = 0;
    dot_robot.dq3 = 0;
    dot_robot.dq4 = 0;
    dot_robot.dq5 = 0;
    dot_robot.dq6 = 0;
}

void get_state(state* s)
{
    pthread_mutex_lock(&mux);
    s->q1 = global_robot.q1;
    s->q2 = global_robot.q2;
    s->q3 = global_robot.q3;
    s->q4 = global_robot.q4;
    s->q5 = global_robot.q5;
    s->q6 = global_robot.q6;
    s->energy = global_robot.energy;
    s->dt3 = global_robot.dt3;
    pthread_mutex_unlock(&mux);
}

void set_state(state s)
{
    pthread_mutex_lock(&mux);
    global_robot.q1 = s.q1;
    global_robot.q2 = s.q2;
    global_robot.q3 = s.q3;
    global_robot.q4 = s.q4;
    global_robot.q5 = s.q5;
    global_robot.q6 = s.q6 ;
    global_robot.energy = s.energy;
    global_robot.dt3 = s.dt3;
    pthread_mutex_unlock(&mux);
}


//-----------------------------------------------------
// The following Functions manage Dynamic Integration
// Interval
//-----------------------------------------------------
void set_dyn_dt(float t)
{
    pthread_mutex_lock(&mux_dt);
    delta_t = t;
    pthread_mutex_unlock(&mux_dt);
}

float get_dyn_dt()
{
    float ret;

    pthread_mutex_lock(&mux_dt);
    ret = delta_t;
    pthread_mutex_unlock(&mux_dt);
    return ret;
}


//-----------------------------------------------------
// The following Function updates the coefficients for
// Trajectory Generation.
// The Trajectory is a fifth grade polynomial
//-----------------------------------------------------
void update_coefficients(float coef1[4], float coef2[4], state robot)
{
    // Actual Initial Position 
    float qi_1 = robot.q4, qi_2 = robot.q5;
    // Desired Final Position controlled by Crawler Task
    float qf_1 = qd.q1d, qf_2 = qd.q2d;

    // First Joint coefficients update
    coef1[0] = qi_1;
    coef1[1] =  10*(qf_1 - qi_1)/pow(T, 3);
    coef1[2] = -15*(qf_1 - qi_1)/pow(T, 4);
    coef1[3] =   6*(qf_1 - qi_1)/pow(T, 5);

    // Second Joint coefficients update
    coef2[0] = qi_2;
    coef2[1] =  10*(qf_2 - qi_2)/pow(T, 3);
    coef2[2] = -15*(qf_2 - qi_2)/pow(T, 4);
    coef2[3] =   6*(qf_2 - qi_2)/pow(T, 5);

}


//-----------------------------------------------------
// The following Function computes the desired joint
// position, velocity and acceleration at time instant t 
//------------------------------------------------------
void compute_qdt(float qdt[2], float dot_qdt[2], float dotdot_qdt[2], float coef1[4], float coef2[4], float t)
{
    if(t <= T){

        // Desired Joint Position
        qdt[0] = coef1[0] + coef1[1]*pow(t, 3) + coef1[2]*pow(t, 4) + coef1[3]*pow(t, 5);
        qdt[1] = coef2[0] + coef2[1]*pow(t, 3) + coef2[2]*pow(t, 4) + coef2[3]*pow(t, 5);

        // Desired Joint Velocity
        dot_qdt[0] = 3*coef1[1]*pow(t, 2) + 4*coef1[2]*pow(t, 3) + 5*coef1[3]*pow(t, 4);
        dot_qdt[1] = 3*coef2[1]*pow(t, 2) + 4*coef2[2]*pow(t, 3) + 5*coef2[3]*pow(t, 4);

        // Desired Joint Acceleration
        dotdot_qdt[0] = 6*coef1[1]*t + 12*coef1[2]*pow(t, 2) + 20*coef1[3]*pow(t, 3);
        dotdot_qdt[1] = 6*coef2[1]*t + 12*coef2[2]*pow(t, 2) + 20*coef2[3]*pow(t, 3);
    } else {
        qdt[0] = qd.q1d;
        qdt[1] = qd.q2d;
        dot_qdt[0] = dot_qdt[1] = 0;
        dotdot_qdt[0] = dotdot_qdt[1] = 0;
    }

    return;
}


//-----------------------------------------------------
// The following Function computes the control torques  
// tau to be used in Dynamic Intergration.
// For more details please see the report  
//-----------------------------------------------------
void generate_tau(float tau[2], state robot, float M[2][2], float C[2][2], float G[2])
{
    static int step;
    static float coefficients1[4], coefficients2[4];
    float e[2], dot_e[2], s[2];     //  Error Variables
    float q_t[2];                   //  Current Joints Position(t)
    float dot_q_t[2];               //  Current Joints Velocity(t)
    float dot_qr_t[2];              //  Target Joints Velocity(t) 
    float dotdot_qr_t[2];           //  Target Joints Acceleration(t)
    float qd_t[2];                  //  Desired Joints Position(t)
    float dot_qd_t[2];              //  Desired Joints Velocity(t) 
    float dotdot_qd_t[2];           //  Desired Joints Acceleration(t) 
    float t;                        //  Time Instant

    // Support Variables for Linear Algebra Operation
    float ris1[2], ris2[2], ris3[2], ris4[2], ris5[2], ris6[2];

    get_desired_joint(&qd);

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
    matvec_mul(dot_qr_t, ris4, 2, 2, C);
    vector_sum(ris4, ris3, ris5, 2);
    matvec_mul(dotdot_qr_t, ris6, 2, 2, M);
    vector_sum(ris6, ris5, tau, 2);

    step++;
}


//------------------------------------------------------
// The following Function manages the numerical error
// due to Forward Euler Integration
//------------------------------------------------------
float adjust_alpha(float x_p, float y_p)
{
    float theta1, theta2;

    y_p -= WHEEL_RADIUS;
    x_p += ROBOT_LENGHT/2; 
    theta1 = atan2(y_p, x_p);
    theta2 = atan2(WHEEL_RADIUS, x_p);
    return (theta1 + theta2);
}


//------------------------------------------------------
// Dynamic Task
//------------------------------------------------------
void* dynamics(void* arg)
{
    printf("DYNAMIC: task started\n");
    int i, exec;            // thread index and system state
    float y_ee;
    float dt;
    float theta;

    init_state();
    state robot;
    get_state(&robot);

    // Dynamic Matrices and Vectors
    float M[2][2], C[2][2], M_inv[2][2], Tsee[4][4], S[4][2];
    float G[2], tau[2];

    // Vectors for dynamic integration: 
    // Subscript 1 => time instant k 
    // Subscript 2 => time instant k+1
    float q_ind1[2], q_ind2[2], qdot_ind1[2], qdot_ind2[2];     // Independent Configuration Variables
    float q_dip1[4], q_dip2[4], qdot_dip1[4], qdot_dip2[4];     // Dependent COnfiguration Variables
    float ris1[2], ris2[2], ris3[2], ris4[4];                   // Support Variables for Linear Algebra Operation
    float qdotdot_ind[2];                                       // Dynamic Integration Result

    // Vector Initialization
    vector_set_zero(q_ind1, 2);
    vector_set_zero(q_dip1, 4);
    vector_set_zero(qdot_ind1, 2);
    vector_set_zero(qdot_dip1, 4);

    i = pt_get_index(arg);
    pt_set_activation(i);
    
    while(get_sys_state(&exec) != STOP){
        if(exec == PLAY){

            dt = get_dyn_dt();
            update_kyn(Tsee, robot, FALSE_ALPHA);
            y_ee = Tsee[1][3];  // End-Effector Position
            

            if(y_ee > 0){
                update_M1(M, robot);
                update_C1(C, robot, dot_robot);
                update_G1(G, robot);

            }else{
                update_kyn(Tsee, robot, TRUE_ALPHA);
                y_ee = Tsee[1][3];

                if(y_ee > 0){
                    // Numerical Error Fix
                    //printf("DYN: Sono dentro y < 0 ma il secondo check dice y=%f\n", y_ee);
                    theta = adjust_alpha(Tsee[0][3], Tsee[1][3]); 
                    robot.q3 -= theta;  
                    q_dip1[2] = robot.q3;
                }

                update_M2(M, robot);
                update_C2(C, robot, dot_robot);
                update_G2(G, robot);

            }

            generate_tau(tau, robot, M, C, G);
            matrix_inverse(M, M_inv);
            
            //qdotdot_ind = M_inv*(tau-G-C*qdot_ind1);    // Current Independent Variables Acceleration
            vector_sub(tau, G, ris1, 2);
            matvec_mul(qdot_ind1, ris2, 2, 2, C);
            vector_sub(ris1, ris2, ris3, 2);
            matvec_mul(ris3, qdotdot_ind, 2, 2, M_inv);    

            //qdot_ind2 = qdot_ind1 + dt*qdotdot_ind;     // Current Independent Variables Velocity
            vector_scal(qdotdot_ind, dt, ris1, 2);
            vector_sum(qdot_ind1, ris1, qdot_ind2, 2);

            //q_ind2 = q_ind1 + dt*qdot_ind1;             // Current Independent Variables Position
            vector_scal(qdot_ind1, dt, ris1, 2);
            vector_sum(q_ind1, ris1, q_ind2, 2);

            //q_dip2 = q_dip1 + dt*qdot_dip1;             // Current Dependent Variables Position
            vector_scal(qdot_dip1, dt, ris4, 4);
            vector_sum(q_dip1, ris4, q_dip2, 4);

            if(y_ee > 0)
                matrix_set_zero(4, 2, S);
            else{
                robot.q3 = q_dip2[2];
                robot.q4 = q_ind2[0];
                robot.q5 = q_ind2[1];
                update_S2(S, robot);
            }
                    
            //qdot_dip2 = S*qdot_ind2;                    // Current Dependent Variables Velocity
            matvec_mul(qdot_ind2, qdot_dip2, 4, 2, S);

            // Update variables: k = k + 1 
            vector_copy(q_dip2, q_dip1, 4);
            vector_copy(q_ind2, q_ind1, 2);
            vector_copy(qdot_dip2, qdot_dip1, 4);
            vector_copy(qdot_ind2, qdot_ind1, 2);


            // Update of Local Robot State
            robot.q1 = q_dip2[0];
            robot.q2 = q_dip2[1];
            robot.q3 = q_dip2[2];

            // Check on ALPHA and Y */
            if(robot.q3 < 0){
                q_dip1[2] = 0;
                robot.q3 = 0;
            }

            if(robot.q2 != 0){
                q_dip1[1] = 0;
                robot.q2 = 0;   
            }

            // Update for next integration step
            robot.q4 = q_ind2[0];
            robot.q5 = q_ind2[1];
            robot.q6 = q_dip2[3];
            dot_robot.dq1 = qdot_dip2[0]; 
            dot_robot.dq2 = qdot_dip2[1];
            dot_robot.dq3 = qdot_dip2[2];
            dot_robot.dq4 = qdot_ind2[0];
            dot_robot.dq5 = qdot_ind2[1];
            dot_robot.dq6 = qdot_dip2[3];

            // Update the global Robot State accessed from the other tasks
            set_state(robot);
        }

        if(exec == RESET){
            // Update the local state with global state after a RESET action from user
            get_state(&robot);
            vector_set_zero(q_ind1, 2);
            vector_set_zero(q_dip1, 4);
            vector_set_zero(qdot_ind1, 2);
            vector_set_zero(qdot_dip1, 4);
        }

        pt_deadline_miss(i);
        pt_wait_for_period(i);
    }

    printf("DYNAMIC: task finished\n");
    return NULL;
}
