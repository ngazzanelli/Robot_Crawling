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


void dynamics(void* arg){






}