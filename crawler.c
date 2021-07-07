#include "qlearn.h"
#include "ptask.h"

#define NSTATES     36
#define NACTIONS    4

typedef struct {
    float q1d;
    float q2d;
    int flag;
} target;
static target qd;
static pthread_mutex_t mux = PTHREAD_MUTEX_INITIALIZER;

extern float get_energy();

void get_desired_joint(target* t){
  pthread_mutex_lock(&mux);
  t->q1d = qd.q1d;
  t->q2d = qd.q2d;
  t->flag = qd.flag;
  qd.flag = 0;
  pthread_mutex_unlock(&mux);
}

// Learning loop 
float qlearn()
{
  int s, a, r, snew;
  long step = 0;
  float newerr, err = 0;
      s = angles2state(rob.th1, rob.th2);
      //printf("Inizio il ciclo while dallo stato s = %d\n", s);
      while (!stop){
          step++;
          a = ql_egreedy_policy(s);
          //printf("Ottenuta l'azione\n");
          snew = next_state(a);
          //printf("Ottenuto il nuovo stato\n");
          r = get_reward(s, snew);
          //printf("Ottenuto il reward\n");
          newerr = ql_updateQ(s, a, r, snew);
          //printf("Aggioranta matrice Q\n");
          err +=  (newerr - err)/step;
          s = snew;
          if (step%1000 == 0){
              ql_reduce_exploration();
              update_info(step, rob.space);
          }
          // Monitoriamo la presenza del ciclo sugli stati
          //if (step > MAXSTEP*0.99)
          //    printf("Stato s = %d\n", s);
          //printf("Finito ciclo %ld\n", step);
          if(view){
              display_links(rob);
              //usleep(1000);
          }
          read_key();
      }
      return err;
}


int main(){
    srand(time(NULL));
    init_state();
    init_graphics();
    ql_init(NSTATES, NACTIONS);  //36 states, 4 actions
    task_create(dynamics, 1, PER, DL, PRI);
    task_create(graphics, 2, PER, DL, PRI);
    task_create(interpreter, 3, PER, DL, PRI);
    task_create(learning, 4, PER, DL, PRI);
    for(i=0; i<=MAX_BALLS; i++){
		  wait_for_task_end(i);
		  printf("fine ciclo %d\n", i);
	  }
}