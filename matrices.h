#ifndef	MATRICES_H
#define	MATRICES_H

typedef struct {
    float q1;
    float q2;
    float q3;
    float q4;
    float q5;
    float q6;
    /*float dotq1;
    float dotq2;
    float dotq3;
    float dotq4;
    float dotq5;
    float dotq6;*/
    float energy;
    float dt3;
} state;

typedef struct {
    float dq1;
    float dq2;
    float dq3;
    float dq4;
    float dq5;
    float dq6;
} dot_state;

extern void update_kyn(float Tsee[4][4], state robot);
extern void update_S2(float S2[4][2], state robot);
extern void update_M1(float M1[2][2], state robot);
extern void update_G1(float G1[2], state robot);
extern void update_C1(float C1[2][2], state robot, dot_state dot_robot);
extern void update_G2(float G2[2], state robot);
extern void update_M2(float M2[2][2], state robot);
extern void update_C2(float C2[2][2], state robot, dot_state dot_robot);
extern float* sum(float *a, float *b, float *c, int dim);
extern float* sub(float *a, float *b, float *c, int dim);
extern float* scal(float *a, float b, float *c, int dim);
extern float* mul(float *x, float *y, int d1, int d2, float A[d1][d2]);
extern void print_matrix(int row, int column, float m[row][column]);
extern void matrix_set_zero(int row, int column, float m[row][column]);

#endif