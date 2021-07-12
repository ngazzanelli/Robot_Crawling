#include <stdio.h>

#define NSTATES     36
#define NACTIONS    4

// Global definitions
#define TH1UP   0       // action move up link 1
#define TH1DW   1       // action move down link 1
#define TH2UP   2       // action move up link 2
#define TH2DW   3       // action move down link 2
#define TH1MAX  1.55    // max theta1 angle [rad]
#define TH1MIN  -0.55   // min theta1 angle [rad]
#define TH2MAX  1.05    // max theta2 angle [rad]
#define TH2MIN  -1.05   // min theta2 angle [rad]
#define DTH1    0.35    // theta1 quantization step [rad]
#define DTH2    0.35    // theta2 quantization step [rad]
#define PRI     10
#define DL      20
#define PER     20

// Global variables (inutili se tanto ho le costanti, no?)
static float t1min, t2min;
static float t1max, t2max;
static float dt1, dt2;
static int n2;          // # of theta2 quantiations 

void init_global_variables()
{
    t1min = TH1MIN;
    t1max = TH1MAX;
    t2min = TH2MIN;
    t2max = TH2MAX;
    dt1 = DTH1;
    dt2 = DTH2;
    n2 = (t2max - t2min)/dt2 + 1;
}
//-------------------------------------------

// Angles to state (non serve la matrice T di transizione dello stato)
int angles2state(float t1, float t2)
{
int i, j;

    i = (t1 - t1min)/dt1;
    j = (t2 - t2min)/dt2;
    //printf("i=%d, j=%d, n2=%d, ret = %d\n", i, j, n2, i*n2+j);
    return i*n2 + j;
}

int main(){
	int s, a, b;
	a = 1;
	b = 2;

	init_global_variables();
    s = angles2state(TH1MAX-a*DTH1, TH2MAX-b*DTH2);
    printf("lo stato è %d\n", s);
}