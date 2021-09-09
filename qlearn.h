#ifndef QLEARN_H
#define QLEARN_H

extern float frand(float xmin, float xmax);

// Functions used by Crawler Task
extern void ql_copy_Q();
extern void ql_init(int ns, int na);
extern void ql_reduce_exploration();
extern float ql_maxQ(int s);
extern float ql_best_action(int s);
extern int  ql_egreedy_policy (int s);
extern float  ql_updateQ(int s, int a, int r, int snew);
extern void ql_print_Qmatrix();

// Functions used by Interpreter Task
extern void ql_inc_learning_rate(float alpha);
extern void ql_inc_discount_factor(float gam);
extern void ql_inc_expl_decay(float decay);
extern void ql_inc_epsini(float eps_ini);
extern void ql_inc_epsfin(float eps_fin);
extern void ql_dec_learning_rate(float alpha);
extern void ql_dec_discount_factor(float gam);
extern void ql_dec_expl_decay(float decay);
extern void ql_dec_epsini(float eps_ini);
extern void ql_dec_epsfin(float eps_fin);

// Functions used by Graphic Task
extern float ql_get_learning_rate();
extern float ql_get_discount_factor();
extern float ql_get_epsini();
extern float ql_get_epsfin();
extern float ql_get_expl_decay();
extern float ql_get_epsilon();
extern void ql_get_Q(float *dest);

extern void ql_set_learning_rate(float alpha);
extern void ql_set_discount_factor(float gam);
extern void ql_set_expl_range(float eps_ini, float eps_fin);
extern void ql_set_expl_decay(float decay);
extern void ql_set_epsini(float eps_ini);
extern void ql_set_epsfin(float eps_fin);
extern void ql_set_epsilon(float e);

#endif