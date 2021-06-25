#include "qlearn.h"
#include "ptask.h"

#define NSTATES     36
#define NACTIONS    4

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