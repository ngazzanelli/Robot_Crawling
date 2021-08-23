#ifndef	MATRICES_H
#define	MATRICES_H

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

extern float* vector_sum(float *a, float *b, float *c, int dim);
extern float* vector_sub(float *a, float *b, float *c, int dim);
extern float* vector_scal(float *a, float b, float *c, int dim);
extern float* matvec_mul(float *x, float *y, int d1, int d2, float A[d1][d2]);
extern void matrix_print(int row, int column, float m[row][column]);
extern void vector_print(int column, float v[column]);
extern void matrix_set_zero(int row, int column, float m[row][column]);
extern void matrix_inverse(float A[2][2], float res[2][2]);

#endif