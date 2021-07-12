#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

typedef struct {
    float q1;
    float q2;
    float q3;
    float q4;
    float q5;
    float q6;
    float dotq1;
    float dotq2;
    float dotq3;
    float dotq4;
    float dotq5;
    float dotq6;
    float energy;
    float dt3;
} state;

void update_kyn(gsl_matrix *Tsee, state robot){
	float value;
	float q1, q2, q3, q4, q5, q6;
	q1 = robot.q1;
	q2 = robot.q2;
	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;
	q6 = robot.q6;


	value = cos(q4 + q5)*sin(q3) + cos(q3)*sin(q4 + q5);
	gsl_matrix_set(Tsee, 0, 0, value);

  	value = cos(q3)*cos(q4 + q5) - sin(q3)*sin(q4 + q5);
  	gsl_matrix_set(Tsee, 0, 1, value);

  	value = 0.0;
  	gsl_matrix_set(Tsee,0,2,value);


  	value = 1.0*(-0.1 + q1 + 0.1*cos(q3) - 0.03*sin(q3)) +
			cos(q3)*(0.1 + 0.06*cos(q5)*sin(q4) +
			cos(q4)*(0.06 + 0.06*sin(q5))) -
			sin(q3)*(0.015 - 0.06*cos(q4)*cos(q5) +
			sin(q4)*(0.06 + 0.06*sin(q5)));
	gsl_matrix_set(Tsee, 0, 3, value);


	value = 0.0 - cos(q3)*cos(q4 + q5) + sin(q3)*sin(q4 + q5);
	gsl_matrix_set(Tsee, 1, 0, value);

	value = cos(q4 + q5)*sin(q3) + cos(q3)*sin(q4 + q5);
	gsl_matrix_set(Tsee, 1, 1, value);

	value = 0.0;
	gsl_matrix_set(Tsee, 1, 2, value);

	value = 1.0*(0.015 + q2 + 0.03*cos(q3) + 0.1*sin(q3)) + 
		sin(q3)*(0.1 + 0.06*cos(q5)*sin(q4) +
		cos(q4)*(0.06 + 0.06*sin(q5))) +
		cos(q3)*(0.015 - 0.06*cos(q4)*cos(q5) +
		sin(q4)*(0.06 + 0.06*sin(q5)));
	gsl_matrix_set(Tsee, 1, 3, value);

	//Terza e quarta riga
	gsl_matrix_set(Tsee, 2, 0, 0.0);
	gsl_matrix_set(Tsee, 2, 1, 0.0);
	gsl_matrix_set(Tsee, 2, 2, 1.0);
	gsl_matrix_set(Tsee, 2, 3, 0.0);

	gsl_matrix_set(Tsee, 3, 0, 0.0);
	gsl_matrix_set(Tsee, 3, 1, 0.0);
	gsl_matrix_set(Tsee, 3, 2, 0.0);
	gsl_matrix_set(Tsee, 3, 3, 1.0);
}

void update_S2(gsl_matrix *S2, state robot){
	float value;
	float q3, q4, q5;
	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;

	//PRIMA RIGA
	value = -0.06*cos(q3 + q4 + q5) +
	0.06*sin(q3 + q4) + ((0.000225 + 0.00045*cos(q3) -
	0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) + 0.0009*sin(q4) +
	0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
	0.06*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
	cos(q4)*(0.0009 + 0.0009*sin(q5)));
	gsl_matrix_set(S2, 0, 0, value);

	value =  -0.06*cos(q3 + q4 + q5) - 
	(0.03*(0.5 + 1.0*cos(q3) - 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) +
	2.0*sin(q4) + 2.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5))/(3.33333 + 1.0*cos(q5)*sin(q4) +
	cos(q4)*(1.0 + 1.0*sin(q5)));
	gsl_matrix_set(S2, 0, 1, value);

	//SECONDA RIGA
	value = ((0.0015 - 0.0015*cos(q3) +
	0.00045*sin(q3))*(-0.06*cos(q3 + q4) - 0.06*sin(q3 + q4 + q5)))/(
	0.003 + 0.0009*cos(q5)*sin(q4) +
	cos(q4)*(0.0009 + 0.0009*sin(q5)));
	gsl_matrix_set(S2, 1, 0, value);

	value = (0.1*(-1.0 + 1.0*cos(q3) - 0.3*sin(q3))*sin(q3 + q4 + q5))/(
	3.33333 + 1.0*cos(q5)*sin(q4) +
	cos(q4)*(1.0 + 1.0*sin(q5)));
	gsl_matrix_set(S2, 1, 1, value);

	//TERZA RIGA
	value = (-1.0*cos(q3 + q4) -
	1.0*sin(q3 + q4 + q5))/(
	3.33333 + 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)));
	gsl_matrix_set(S2, 2, 0, value);

	value =  -((1.0*sin(q3 + q4 + q5))/(
	3.33333 + 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5))));
	gsl_matrix_set(S2, 2, 1, value);

	//QUARTA RIGA
	value = (3.0*cos(q3 + q4) - 4.0*cos(q3 - q5) +
	4.0*cos(q3 + q5) + (13.3333)*cos(q3 + q4 + q5) -
	8.0*sin(q3) - (13.3333)*sin(q3 + q4) +
	3.0*sin(q3 + q4 + q5))/(
	3.33333 + 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)));
	gsl_matrix_set(S2, 3, 0, value);


	value = (4.0*cos(q3 + q5) + (13.3333)*cos(q3 + q4 + q5) -
	4.0*sin(q3) + 3.0*sin(q3 + q4 + q5))/(3.33333 + 1.0*cos(q5)*sin(q4) +
	cos(q4)*(1.0 + 1.0*sin(q5)));
	gsl_matrix_set(S2, 3, 1, value);
}


void update_M1(gsl_matrix *M1, state robot){
	float value;
	float q5;
	q5 = robot.q5;

	//PRIMA RIGA
	value = 0.0000822813 + 0.000018*sin(q5);
	gsl_matrix_set(M1, 0, 0, value);

	value = 0.0000231406 + 9.*pow(10, -6)*sin(q5);
	gsl_matrix_set(M1, 0, 1, value);

	//SECONDA RIGA
	value = 0.0000231406 + 9.*pow(10, -6)*sin(q5);
	gsl_matrix_set(M1, 1, 0, value);

	value = 0.0000231406;
	gsl_matrix_set(M1, 1, 1, value);
}

void update_G1(gsl_vector *G1, state robot){
	float value;
	float q3, q4, q5;
	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;
	
	//PRIMA RIGA
	value = 0.0 + 0.00441*cos(q3)*cos(q4) + 0.00294*cos(q3 + q4) -
	0.00441*sin(q3)*sin(q4) + 0.00147*sin(q3 + q4 + q5);
	gsl_vector_set(G1, 0, value);

	//SECONDA RIGA
	value = 0.0 + 0.00147*sin(q3 + q4 + q5);
	gsl_vector_set(G1, 1, value);
}

void update_C1(gsl_matrix *C1, state robot, dot_state dot_robot){
	float value;
	float q3, q4, q5, dotq1, dotq2, dotq3, dotq4, dotq5;
	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;
	dotq1 = dot_robot.dotq1;
	dotq2 = dot_robot.dotq2;
	dotq3 = dot_robot.dotq3;
	dotq4 = dot_robot.dotq4;
	dotq5 = dot_robot.dotq5;

	//PRIMA RIGA
	value = (0.000375*dotq1 - 5.625*pow(10, -6)*dotq3)*cos(q3 + q4) +
	0.000018*dotq5*cos(q5) - 7.5*pow(10, -6)*dotq3*cos(q4 + q5) -
	0.000075*dotq2*cos(q3 + q4 + q5) -
	7.5*pow(10, -6)*dotq3*cos(q3 + q4 + q5) +
	cos(q4)*(-0.00001125*dotq3 + 0.000225*dotq2*sin(q3)) +
	0.0000375*dotq3*sin(q4) + 0.000225*dotq2*cos(q3)*sin(q4) +
	0.00015*dotq2*sin(q3 + q4) + 0.0000375*dotq3*sin(q3 + q4) -
	2.25*pow(10, -6)*dotq3*sin(q4 + q5) + 0.000075*dotq1*sin(q3 + q4 + q5) -
	1.125*pow(10, -6)*dotq3*sin(q3 + q4 + q5);
	gsl_matrix_set(C1, 0, 0, value);

	value = 9.*pow(10, -6)*dotq5*cos(q5) - 7.5*pow(10, -6)*dotq3*cos(q4 + q5) -
	0.000075*dotq2*cos(q3 + q4 + q5) -
	7.5*pow(10, -6)*dotq3*cos(q3 + q4 + q5) -
	2.25*pow(10, -6)*dotq3*sin(q4 + q5) + 0.000075*dotq1*sin(q3 + q4 + q5) -
	1.125*pow(10, -6)*dotq3*sin(q3 + q4 + q5);
	gsl_matrix_set(C1, 0, 1, value);

	//SECONDA RIGA
	value = -9.*pow(10, -6)*dotq4*cos(q5) +
	4.5*pow(10, -6)*dotq5*cos(q5) - 0.000075*dotq2*cos(q3 + q4 + q5) -
	4.5*pow(10, -6)*dotq3*(1.0*cos(q3 - q5) + 1.0*cos(q3 + q5) +
	1.66667*cos(q4 + q5) + 1.66667*cos(q3 + q4 + q5) +
	0.5*sin(q4 + q5) + 0.25*sin(q3 + q4 + q5)) +
	0.000075*dotq1*sin(q3 + q4 + q5);
	gsl_matrix_set(C1, 1, 0, value);

	value = -4.5*pow(10, -6)*dotq4*cos(q5) -
	4.5*pow(10, -6)*dotq3*cos(q3 + q5) - 7.5*pow(10, -6)*dotq3*cos(q4 + q5) -
	0.000075*dotq2*cos(q3 + q4 + q5) -
	7.5*pow(10, -6)*dotq3*cos(q3 + q4 + q5) -
	2.25*pow(10, -6)*dotq3*sin(q4 + q5) + 0.000075*dotq1*sin(q3 + q4 + q5) -
	1.125*pow(10, -6)*dotq3*sin(q3 + q4 + q5);
	gsl_matrix_set(C1, 1, 1, value);
}

void update_G2(gsl_vector *G2, state robot){
	float value;
	float q3, q4, q5;
	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;

	//PRIMA RIGA
	value = 0.0 + 0.00441*cos(q3)*cos(q4) + 0.00294*cos(q3 + q4) -
	0.00441*sin(q3)*sin(q4) + ((-0.11074 + 0.4792*cos(q3) - 
	0.143766*sin(q3) + 0.00147*cos(q5)*sin(q4) +
	cos(q4)*(0.00735 + 0.00147*sin(q5)))*(-1.0*cos(q3 + q4) -
	1.0*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) + 
	cos(q4)*(1.0 + 1.0*sin(q5))) + (4.7922*(0.0015 - 0.0015*cos(q3) +
	0.00045*sin(q3))*(-0.06*cos(q3 + q4) -
	0.06*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
	cos(q4)*(0.0009 + 0.0009*sin(q5))) +
	0.00147*sin(q3 + q4 + q5);
	gsl_vector_set(G2, 0, value);

	//SECONDA RIGA
	value = 0.0 + 0.00147*sin(q3 + q4 + q5) + 
	 (0.47922*(-1.0 + 1.0*cos(q3) - 0.3*sin(q3))*sin(q3 + q4 + q5))/(
	3.33333 + 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5))) - (
	1.0*(-0.11074 + 0.4792*cos(q3) - 0.143766*sin(q3) +
	0.00147*cos(q5)*sin(q4) +
	cos(q4)*(0.00735 + 0.00147*sin(q5)))*sin(q3 + q4 + q5))/(
	3.33333 + 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)));
	gsl_vector_set(G2, 1, value);
}


void update_M2(gsl_matrix *M2, state robot){
 	float value;
 	float q1, q2, q3, q4, q5, q6;
 	q1 = robot.q1;
 	q2 = robot.q2;
 	q3 = robot.q3;
 	q4 = robot.q4;
 	q5 = robot.q5;
 	q6 = robot.q6;

//M2(0,0)
 	 value = 0.0000822813 +
	 0.000018*sin(q5) + (0.00015*cos(q3 + q4 + q5) - 
	 0.00075*sin(q3 + q4))*(-0.06*cos(q3 + q4 + q5) +
	 0.06*sin(q3 + q4) + ((0.000225 + 0.00045*cos(q3) -
	 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) + 0.0009*sin(q4) +
	 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
	 0.06*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
	 cos(q4)*(0.0009 + 0.0009*sin(q5)))) + ((-0.06*cos(q3 + q4 + q5) +
	 0.06*sin(q3 + q4) + ((0.000225 + 0.00045*cos(q3) -
	 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) + 0.0009*sin(q4) +
	 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
	 0.06*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
	 cos(q4)*(0.0009 + 0.0009*sin(q5))))*(-0.011025*cos(q3 + q4) + 0.02889*cos(q3 - q5) - 0.02889*cos(q3 + q5) -
	 0.0973*cos(q3 + q4 + q5) + 0.05778*sin(q3) +
	 0.0953*sin(q3 + q4) - 0.011025*sin(q3 + q4 + q5)))/(3.33333 +
	 1.0*cos(q5)*sin(q4) +
	 cos(q4)*(1.0 + 1.0*sin(q5))) + ((-1.0*cos(q3 + q4) -
	 1.0*sin(q3 + q4 + q5))*(0.0000372812 + 0.000045*cos(q3) +
	 0.000075*cos(q4) + 0.000075*cos(q3 + q4) -
	 4.5*pow(10, -6)*cos(q4 + q5) - 2.25*pow(10, -6)*cos(q3 + q4 + q5) +
	 0.0000225*sin(q4) + 0.00001125*sin(q3 + q4) -
	 9.*pow(10, -6)*sin(q3 - q5) + 9.*pow(10, -6)*sin(q3 + q5) +
	 0.000015*sin(q4 + q5) +
	 0.000015*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) +
	 cos(q4)*(1.0 + 1.0*sin(q5))) + ((-1.0*cos(q3 + q4) -
	 1.0*sin(q3 + q4 + q5))*(0.0000372812 + 0.000045*cos(q3) +
	 0.000075*cos(q4) + 0.000075*cos(q3 + q4) -
	 4.5*pow(10, -6)*cos(q4 + q5) - 2.25*pow(10, -6)*cos(q3 + q4 + q5) +
	 0.0000225*sin(q4) + 0.00001125*sin(q3 + q4) -
	 9.*pow(10, -6)*sin(q3 - q5) + 9.*pow(10, -6)*sin(q3 + q5) +
	 0.000015*sin(q4 + q5) + (0.00369 - 0.01467*cos(q3) +
	 0.00015*cos(q4)*cos(q5) - 0.0489*sin(q3) -
	 0.00075*sin(q4) -
	 0.00015*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4 + q5) +
	 0.06*sin(q3 + q4) + ((0.000225 + 0.00045*cos(q3) -
	 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
	 0.0009*sin(q4) +
	 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
	 0.06*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
	 cos(q4)*(0.0009 + 0.0009*sin(q5)))) + ((0.00728604 -
	 0.0024814*cos(q3) + 0.00015*cos(q3 - q4) +
	 0.00015*cos(q4) - 9.*pow(10, -6)*cos(q3 - q4 - q5) -
	 4.5*pow(10, -6)*cos(q4 + q5) - 0.00006*sin(q3) -
	 0.000045*sin(q3 - q4) + 0.0000225*sin(q4) -
	 0.00003*sin(q3 - q4 - q5) + 0.000018*sin(q5) +
	 0.00003*sin(q4 + q5))*(-1.0*cos(q3 + q4) -
	 1.0*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) +
	 cos(q4)*(1.0 + 1.0*sin(q5))) + ((0.0015 - 0.0015*cos(q3) +
	 0.00045*sin(q3))*(-0.0113 + 0.0489*cos(q3) -
	 0.01467*sin(q3) + 0.00015*cos(q5)*sin(q4) +
	 cos(q4)*(0.00075 + 0.00015*sin(q5)))*(-0.06*cos(q3 + q4) - 0.06*sin(q3 + q4 + q5)))/(0.003 +
	 0.0009*cos(q5)*sin(q4) + 
	 cos(q4)*(0.0009 + 0.0009*sin(q5))) + (0.0000448875*cos(q3 + q4) - 0.00005985*cos(q3 - q5) +
	 0.00005985*cos(q3 + q5) + 0.0001995*cos(q3 + q4 + q5) -
	 0.0001197*sin(q3) - 0.0001995*sin(q3 + q4) +
	 0.0000448875*sin(q3 + q4 + q5))/(3.33333 +
	 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5))) +
	 0.000015*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) +
	 cos(q4)*(1.0 + 1.0*sin(q5))) + ((0.0015 - 0.0015*cos(q3) +
	 0.00045*sin(q3))*(-0.06*cos(q3 + q4) -
	 0.06*sin(q3 + q4 + q5))*(0.00045*cos(q3)*cos(q4) +
	 0.0003*cos(q3 + q4) - 0.00045*sin(q3)*sin(q4) +
	 0.00015*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
	 cos(q4)*(0.0009 + 0.0009*sin(q5))) + (0.000045*(-1.0 +
	 1.0*cos(q3) - 0.3*sin(q3))*(-78.0*cos(q3 + q4) -
	 0.666667*sin(q3 - q5) - 0.666667*sin(q3 + q5) -
	 82.4444*sin(q3 + q4 + q5))*(1.0*cos(q3 + q4) +
	 1.0*sin(q3 + q4 + q5)))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
	 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + (0.000089775*(1.0*cos(q3 + q4) -
	 2.0*cos(q3 - q5) + 2.0*cos(q3 + q5) +
	 6.66667*cos(q3 + q4 + q5) - 4.0*sin(q3) -
	 6.66667*sin(q3 + q4) +
	 1.0*sin(q3 + q4 + q5))*(1.0*cos(q3 + q4) -
	 1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
	 4.44444*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
	 4.44444*sin(q3 + q4) + 1.0*sin(q3 + q4 + q5)))/pow((3.33333 +
	 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5))), 2);
	 gsl_matrix_set(M2, 0, 0, value);

//M2(0,1)
	value = 1/pow((3.33333 + 1.0*cos(q4) + 1.0*sin(q4 + q5)),2)*(0.0165437 -
	 0.0000372813*cos(q3) - 0.0029547*cos(2*q3) + 0.0100558*cos(q4) +
	 0.0000171953*cos(2*q4) - 0.0000621354*cos(q3 + q4) -
	 7.875*pow(10, -6)*cos(2*(q3 + q4)) - 0.0065685*cos(2*q3 + q4) -
	 9.32031*pow(10, -6)*cos(q3 + 2*q4) + 0.0003906*cos(q4 - q5) -
	 0.0009774*cos(2*q5) + 0.0009774*cos(2*(q3 + q5)) -
	 0.0011718*cos(q4 + q5) - 0.0000104453*cos(2*(q4 + q5)) +
	 0.00870785*cos(2*(q3 + q4 + q5)) +
	 0.0015624*cos(2*q3 + q4 + q5) + 0.002604*cos(2*q3 + 2*q4 + q5) -
	 0.003258*cos(q4 + 2*q5) + 0.006516*cos(2*q3 + q4 + 2*q5) +
	 0.0000279609*cos(q3 + 2*(q4 + q5)) + 0.0011718*sin(q4) -
	 0.0007812*sin(2*q3 + q4) + 9.32031*pow(10, -6)*sin(q3 - q5) +
	 0.000985275*sin(2*q3 - q5) - 0.003258*sin(q4 - q5) +
	 0.016964*sin(q5) - 0.0000279609*sin(q3 + q5) -
	 0.00294683*sin(2*q3 + q5) + 0.0100108*sin(q4 + q5) -
	 0.000186406*sin(q3 + q4 + q5) + 0.002604*sin(2*(q3 + q4 + q5)) -
	 0.0131295*sin(2*q3 + q4 + q5) + 0.0000276406*sin(2*q4 + q5) -
	 0.0000372813*sin(q3 + 2*q4 + q5) -
	 0.00871572*sin(2*q3 + 2*q4 + q5) - 0.0003906*sin(q4 + 2*q5) +
	 0.0007812*sin(2*q3 + q4 + 2*q5));
	 gsl_matrix_set(M2, 0, 1, value);

//M2(1,0)
	value = 1/pow((3.33333 + 1.0*cos(q4) +
 1.0*sin(q4 + q5)), 2)*(0.0165437 - 0.0000372812*cos(q3) -
 0.0029547*cos(2*q3) + 0.0100558*cos(q4) +
 0.0000171953*cos(2*q4) - 0.0000621354*cos(q3 + q4) -
 7.875*pow(10, -6)*cos(2*(q3 + q4)) - 0.0065685*cos(2*q3 + q4) -
 9.32031*pow(10, -6)*cos(q3 + 2*q4) + 0.0003906*cos(q4 - q5) -
 0.0009774*cos(2*q5) + 0.0009774*cos(2*(q3 + q5)) -
 0.0011718*cos(q4 + q5) - 0.0000104453*cos(2*(q4 + q5)) +
 0.00870785*cos(2*(q3 + q4 + q5)) +
 0.0015624*cos(2*q3 + q4 + q5) + 0.002604*cos(2*q3 + 2*q4 + q5) -
 0.003258*cos(q4 + 2*q5) + 0.006516*cos(2*q3 + q4 + 2*q5) +
 0.0000279609*cos(q3 + 2*(q4 + q5)) + 0.0011718*sin(q4) -
 0.0007812*sin(2*q3 + q4) + 9.32031*pow(10, -6)*sin(q3 - q5) +
 0.000985275*sin(2*q3 - q5) - 0.003258*sin(q4 - q5) +
 0.016964*sin(q5) - 0.0000279609*sin(q3 + q5) -
 0.00294682*sin(2*q3 + q5) + 0.0100108*sin(q4 + q5) -
 0.000186406*sin(q3 + q4 + q5) + 0.002604*sin(2*(q3 + q4 + q5)) -
 0.0131295*sin(2*q3 + q4 + q5) + 0.0000276406*sin(2*q4 + q5) -
 0.0000372812*sin(q3 + 2*q4 + q5) -
 0.00871572*sin(2*q3 + 2*q4 + q5) - 0.0003906*sin(q4 + 2*q5) +
 0.0007812*sin(2*q3 + q4 + 2*q5));
 	gsl_matrix_set(M2, 1, 0, value);

//M2(1,1)
	value = (0.0155731 -
 0.0000186406*cos(q3) - 0.00099315*cos(2*q3) +
 0.00685027*cos(q4) + 0.0000250703*cos(2*q4) +
 0.0009774*cos(2*(q3 + q5)) - 0.0007812*cos(q4 + q5) -
 9.32031*pow(10, -6)*cos(2*(q4 + q5)) +
 0.00870897*cos(2*(q3 + q4 + q5)) +
 0.0007812*cos(2*q3 + q4 + q5) + 0.006516*cos(2*q3 + q4 + 2*q5) +
 0.0000186406*cos(q3 + 2*(q4 + q5)) + 0.0007812*sin(q4) +
 0.00198919*sin(q5) - 0.0000186406*sin(q3 + q5) -
 0.00197055*sin(2*q3 + q5) + 0.00674527*sin(q4 + q5) -
 0.000124271*sin(q3 + q4 + q5) + 0.002604*sin(2*(q3 + q4 + q5)) -
 0.006621*sin(2*q3 + q4 + q5) + 0.0000343906*sin(2*q4 + q5) -
 0.0000186406*sin(q3 + 2*q4 + q5) -
 0.00001575*sin(2*q3 + 2*q4 + q5) +
 0.0007812*sin(2*q3 + q4 + 2*q5))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2);
 	gsl_matrix_set(M2, 1, 1, value);
}

void update_C2(gsl_matrix *C2, state robot){

	float value;
	float q1, q2, q3, q4, q5, q6;
    float dotq1, dotq2, dotq3, dotq4, dotq5, dotq6;

    q1 = robot.q1;
	q2 = robot.q2;
	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;
	q6 = robot.q6;
    dotq1 = robot.dotq1;
    dotq2 = robot.dotq2;
    dotq3 = robot.dotq3;
    dotq4 = robot.dotq4;
    dotq5 = robot.dotq5;
    dotq6 = robot.dotq6;

//C2(0,0)
	value = 0.000375*dotq1*cos(q3 + q4) + 0.000018*dotq5*cos(q5) -
 0.000075*dotq2*cos(q3 + q4 + q5) + 0.000375*dotq2*sin(q3 + q4) -
 0.00001125*dotq3*(1.0*cos(q4) + 0.5*cos(q3 + q4) +
 0.666667*cos(q4 + q5) + 0.666667*cos(q3 + q4 + q5) -
 3.33333*sin(q4) - 3.33333*sin(q3 + q4) + 0.2*sin(q4 + q5) +
 0.1*sin(q3 + q4 + q5)) +
 0.000075*dotq1*sin(q3 + q4 +q5) - (0.000045*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(0.166667*dotq3*cos(q4 + q5) + (-0.333333*dotq3 - 0.166667*dotq4 -
 0.166667*dotq5)*cos(q3 + q4 + q5) -
 0.833333*dotq3*sin(q4) + 1.66667*dotq3*sin(q3 + q4) +
 0.833333*dotq4*sin(q3 + q4))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + (0.000075*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5))*(0.3*dotq3*cos(q3 - q4) -
 0.15*dotq3*cos(q3 + q4) - 0.075*dotq4*cos(q3 + q4) +
 0.12*dotq3*cos(q3 - q5) - 0.12*dotq5*cos(q3 - q5) +
 0.2*dotq3*cos(q3 - q4 - q5) - 0.12*dotq3*cos(q3 + q5) -
 0.12*dotq5*cos(q3 + q5) + 0.2*dotq3*cos(q4 + q5) -
 0.1*dotq4*cos(q4 + q5) - 0.1*dotq5*cos(q4 + q5) -
 0.2*dotq3*cos(q3 + q4 + q5) - 0.1*dotq4*cos(q3 + q4 + q5) -
 0.1*dotq5*cos(q3 + q4 + q5) + 0.6*dotq3*sin(q3) +
 1.0*dotq3*sin(q3 - q4) - 5.0*dotq2*sin(q4) - 1.0*dotq3*sin(q4) +
 0.5*dotq4*sin(q4) - 1.0*dotq1*cos(q5)*sin(q4) +
 1.0*dotq3*sin(q3 + q4) + 0.5*dotq4*sin(q3 + q4) -
 0.06*dotq3*sin(q3 - q4 - q5) - 1.0*dotq2*sin(q4)*sin(q5) +
 cos(q4)*(-5.0*dotq1 + 0.15*dotq3 - 0.15*dotq4 +
 1.0*dotq2*cos(q5) - 1.0*dotq1*sin(q5)) +
 0.03*dotq3*sin(q4 + q5) - 0.03*dotq4*sin(q4 + q5) -
 0.03*dotq5*sin(q4 + q5) - 0.03*dotq3*sin(q3 + q4 + q5) -
 0.015*dotq4*sin(q3 + q4 + q5) -
 0.015*dotq5*sin(q3 + q4 + q5)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + (-0.06*cos(q3 + q4 + q5) +
 0.06*sin(q3 + q4) + ((0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) + 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))))*((-0.00075*dotq3 -
 0.000375*dotq4)*cos(q3 + q4) +
 0.000075*dotq3*cos(q5)*sin(q4) +
 dotq3*cos(q4)*(0.000375 + 0.000075*sin(q5)) -
 0.00015*dotq3*sin(q3 + q4 + q5) -
 0.000075*dotq4*sin(q3 + q4 + q5) -
 0.000075*dotq5*sin(q3 + q4 + q5)) + ((3.0*cos(q3 + q4) -
 4.0*cos(q3 - q5) + 4.0*cos(q3 + q5) +
 13.3333*cos(q3 + q4 + q5) - 8.0*sin(q3) -
 13.3333*sin(q3 + q4) +
 3.0*sin(q3 + q4 + q5))*(0.0000149625*dotq4*cos(q3 + q4)*cos(q3 + q4 + q5) + 0.000029925*dotq4*sin(q3) +
 0.0000149625*dotq5*sin(q3) +
 0.000049875*dotq4*cos(q4)*sin(q3) +
 0.000049875*dotq4*cos(q5)*sin(q3)*sin(q4) +
 0.000049875*dotq5*cos(q5)*sin(q3)*sin(q4) -
 7.48125*pow(10, -6)*dotq4*sin(2*(q3 + q4)) +
 0.000029925*dotq4*sin(q3)*sin(q5) +
 0.0000149625*dotq5*sin(q3)*sin(q5) +
 0.000049875*dotq4*cos(q4)*sin(q3)*sin(q5) +
 0.000049875*dotq5*cos(q4)*sin(q3)*sin(q5) +
 0.0000149625*cos(q3)*((-3.33333*dotq4 - 3.33333*dotq5)*cos(q4)*cos(q5) +
 sin(q4)*(3.33333*dotq4 + (3.33333*dotq4 +
 3.33333*dotq5)*sin(q5))) -
 0.0000149625*dotq4*sin(q3 + q4)*sin(q3 + q4 + q5) -
 0.0000149625*dotq5*sin(q3 + q4)*sin(q3 + q4 + q5) +
 0.0001197*(1.0*cos(q3) + 1.66667*cos(q3 + q4) -
 0.375*cos(q3 + q4 + q5) + 0.375*sin(q3 + q4) -
 0.5*sin(q3 - q5) + 0.5*sin(q3 + q5) +
 1.66667*sin(q3 + q4 + q5))*(1.0*dotq4*cos(q3 + q4) + 
 (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5)) +
 0.0000149625*dotq4*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(-13.3333*cos(q3 + q4) +
 3.0*cos(q3 + q4 + q5) - 3.0*sin(q3 + q4) -
 13.3333*sin(q3 + q4 + q5)) -
 3.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(1.0*cos(q3 + q4) -
 1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
 4.44444*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
 4.44444*sin(q3 + q4) + 1.0*sin(q3 + q4 + q5))) +
 0.0000149625*dotq5*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 4.0*sin(q3 - q5) - 4.0*sin(q3 + q5) -
 13.3333*sin(q3 + q4 + q5)) -
 3.0*(1.0*cos(q4)*cos(q5) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) -
 1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
 4.44444*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
 4.44444*sin(q3 + q4) + 1.0*sin(q3 + q4 + q5))) +
 7.48125*pow(10, -6)*dotq4*sin(2*(q3 + q4 + q5)) +
 7.48125*pow(10, -6)*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 3) + ((0.0000372812 +
 0.000045*cos(q3) + 0.000075*cos(q4) + 0.000075*cos(q3 + q4) -
 4.5*pow(10, -6)*cos(q4 + q5) - 2.25*pow(10, -6)*cos(q3 + q4 + q5) +
 0.0000225*sin(q4) + 0.00001125*sin(q3 + q4) -
 9.*pow(10, -6)*sin(q3 - q5) + 9.*pow(10, -6)*sin(q3 + q5) +
 0.000015*sin(q4 + q5) +
 0.000015*sin(q3 + q4 + q5))*(1.0*dotq4*cos(q3 + q4)*cos(q3 + q4 + q5) + 
 2.0*dotq4*sin(q3) + 1.0*dotq5*sin(q3) +
 3.33333*dotq4*cos(q4)*sin(q3) +
 3.33333*dotq4*cos(q5)*sin(q3)*sin(q4) +
 3.33333*dotq5*cos(q5)*sin(q3)*sin(q4) -
 0.5*dotq4*sin(2*(q3 + q4)) + 2.0*dotq4*sin(q3)*sin(q5) +
 1.0*dotq5*sin(q3)*sin(q5) +
 3.33333*dotq4*cos(q4)*sin(q3)*sin(q5) +
 3.33333*dotq5*cos(q4)*sin(q3)*sin(q5) +
 cos(q3)*((-3.33333*dotq4 - 3.33333*dotq5)*cos(q4)*cos(q5) +
 sin(q4)*(3.33333*dotq4 + 3.33333*dotq4*sin(q5) +
 3.33333*dotq5*sin(q5))) -
 1.0*dotq4*sin(q3 + q4)*sin(q3 + q4 + q5) -
 1.0*dotq5*sin(q3 + q4)*sin(q3 + q4 + q5) +
 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((0.00045*cos(q3)*cos(q4) +
 0.0003*cos(q3 + q4) - 0.00045*sin(q3)*sin(q4) +
 0.00015*sin(q3 + q4 + q5))*(-0.1*dotq5*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(sin(q3)*(1.0 +
 3.33333*cos(q5)*sin(q4) + (1.0 + 3.33333*cos(q4))*sin(
 q5)) + cos(q3)*(-3.33333*cos(q4)*cos(q5) +
 3.33333*sin(q4)*sin(q5))) -
 0.1*dotq4*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(sin(q3)*(2.0 + 3.33333*cos(q5)*sin(q4) + 2.0*sin(q5) +
 cos(q4)*(3.33333 + 3.33333*sin(q5))) +
 cos(q3)*(-3.33333*cos(q4)*cos(q5) +
 sin(q4)*(3.33333 + 3.33333*sin(q5)))) +
 0.03*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))*(1.0*cos(2*q3 + q4) +
 3.33333*cos(q3 + q4 + q5) - 3.33333*cos(2*q3 + q4 + q5) -
 3.33333*sin(q3 + q4) + 3.33333*sin(2*q3 + q4) +
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + (0.00015*cos(q3 + q4 + q5) -
 0.00075*sin(q3 + q4))*(dotq5*(-((
 0.06*cos(q3 + q4 + q5)*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) + 
 cos(q4)*(0.0009 + 0.0009*sin(q5)))) + ((0.0009*cos(q5)*sin(q4) + 
 0.0009*cos(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) - ((0.0009*cos(q4)*cos(q5) -
 0.0009*sin(q4)*sin(q5))*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/pow((0.003 +
 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))), 2) +
 0.06*sin(q3 + q4 + q5)) +
 dotq4*(0.06*cos(q3 + q4) + ((-0.06*cos(q3 + q4 + q5) +
 0.06*sin(q3 + q4))*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) + ((0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5)))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))) +
 0.06*sin(q3 + q4 +q5) - (0.06*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(-0.25 - 0.5*cos(q3) +
 1.0*cos(q4)*cos(q5) - 1.66667*sin(q3) - 1.0*sin(q4) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)) + (0.1*(1.0*dotq4*cos(q3 + q4) + 
 (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))*(-1.2*cos(q3) - 
 2.0*cos(q3 + q4) +
 1.0*cos(2*q3 + q4) + 0.15*cos(q3 + q4 + q5) +
 0.3*cos(2*q3 + q4 + q5) - 0.15*sin(q3 + q4) -
 0.3*sin(2*q3 + q4) + 0.6*sin(q3 - q5) - 0.6*sin(q3 + q5) -
 2.0*sin(q3 + q4 + q5) + 1.0*sin(2*q3 + q4 + q5)))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)) + ((0.0015 - 0.0015*cos(q3) +
 0.00045*sin(q3))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5))*(0.00015*dotq3*cos(q3 + q4 + q5) +
 0.00015*dotq4*cos(q3 + q4 + q5) +
 0.00015*dotq5*cos(q3 + q4 + q5) -
 0.00075*dotq3*sin(q3 + q4) -
 0.00075*dotq4*sin(q3 + q4) + (0.01467*(1.0*dotq3*cos(q3) + 
 (-0.0102249*dotq4 - 0.0102249*dotq5)*cos(q4)*cos(q5) + 3.33333*dotq3*sin(q3) +
 0.0511247*dotq4*sin(q4) +
 0.0102249*dotq4*sin(q4)*sin(q5) +
 0.0102249*dotq5*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + ((-0.0113 + 0.0489*cos(q3) -
 0.01467*sin(q3) + 0.00015*cos(q5)*sin(q4) +
 cos(q4)*(0.00075 + 0.00015*sin(q5)))*(1.0*dotq4*cos(q3 + q4)*cos(q3 + q4 + q5) + 
 2.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 3.33333*dotq4*cos(q4)*sin(q3) +
 3.33333*dotq4*cos(q5)*sin(q3)*sin(q4) +
 3.33333*dotq5*cos(q5)*sin(q3)*sin(q4) -
 0.5*dotq4*sin(2*(q3 + q4)) + 2.0*dotq4*sin(q3)*sin(q5) +
 1.0*dotq5*sin(q3)*sin(q5) +
 3.33333*dotq4*cos(q4)*sin(q3)*sin(q5) +
 3.33333*dotq5*cos(q4)*sin(q3)*sin(q5) +
 cos(q3)*((-3.33333*dotq4 - 3.33333*dotq5)*cos(q4)*cos(q5) + sin(q4)*(3.33333*dotq4 + 3.33333*dotq4*sin(q5) +
 3.33333*dotq5*sin(q5))) -
 1.0*dotq4*sin(q3 + q4)*sin(q3 + q4 + q5) -
 1.0*dotq5*sin(q3 + q4)*sin(q3 + q4 + q5) +
 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))), 2) + (0.489*(-0.1*dotq5*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(sin(q3)*(1.0 +
 3.33333*cos(q5)*sin(q4) + (1.0 + 3.33333*cos(q4))*sin(q5)) +
 cos(q3)*(-3.33333*cos(q4)*cos(q5) +
 3.33333*sin(q4)*sin(q5))) -
 0.1*dotq4*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(sin(q3)*(2.0 + 3.33333*cos(q5)*sin(q4) + 2.0*sin(q5) +
 cos(q4)*(3.33333 + 3.33333*sin(q5))) +
 cos(q3)*(-3.33333*cos(q4)*cos(q5) +
 sin(q4)*(3.33333 + 3.33333*sin(q5)))) +
 0.03*(1.0*dotq4*cos(q3 + q4) +
 (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))*(1.0*cos(2*q3 + q4) +
 3.33333*cos(q3 + q4 + q5) -
 3.33333*cos(2*q3 + q4 + q5) - 3.33333*sin(q3 + q4) +
 3.33333*sin(2*q3 + q4) +
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)))/(0.003 +
 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))) + (-0.06*cos(q3 + q4 + q5) +
 0.06*sin(q3 + q4) + ((0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) + 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))))*(dotq3*(-0.00075*cos(q3 + q4) -
 0.00015*sin(q3 + q4 + q5)) +
 dotq4*(-0.00075*cos(q3 + q4) - 0.00015*sin(q3 + q4 + q5)) -
 0.00015*dotq5*sin(q3 + q4 + q5) + (0.0489*(1.0*dotq3*cos(q3) - 0.3*dotq3*sin(q3) +
 0.00306748*dotq4*cos(q5)*sin(q4) +
 0.00306748*dotq5*cos(q5)*sin(q4) +
 cos(q4)*(0.0153374*dotq4 + 0.00306748*dotq4*sin(q5) +
 0.00306748*dotq5*sin(q5)))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + ((0.00369 - 0.01467*cos(q3) +
 0.00015*cos(q4)*cos(q5) - 0.0489*sin(q3) -
 0.00075*sin(q4) -
 0.00015*sin(q4)*sin(q5))*(1.0*dotq4*cos(q3 + q4)*cos(q3 + q4 + q5) + 2.0*dotq4*sin(q3) + 1.0*dotq5*sin(q3) +
 3.33333*dotq4*cos(q4)*sin(q3) +
 3.33333*dotq4*cos(q5)*sin(q3)*sin(q4) +
 3.33333*dotq5*cos(q5)*sin(q3)*sin(q4) -
 0.5*dotq4*sin(2*(q3 + q4)) + 2.0*dotq4*sin(q3)*sin(q5) +
 1.0*dotq5*sin(q3)*sin(q5) +
 3.33333*dotq4*cos(q4)*sin(q3)*sin(q5) +
 3.33333*dotq5*cos(q4)*sin(q3)*sin(q5) +
 cos(q3)*((-3.33333*dotq4 - 3.33333*dotq5)*cos(q4)*cos(q5) +
 sin(q4)*(3.33333*dotq4 + 3.33333*dotq4*sin(q5) +
 3.33333*dotq5*sin(q5))) -
 1.0*dotq4*sin(q3 + q4)*sin(q3 + q4 + q5) -
 1.0*dotq5*sin(q3 + q4)*sin(q3 + q4 + q5) +
 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 + 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5))), 2) +
 0.489*(dotq5*(-((0.06*cos(q3 + q4 + q5)*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5)))) + ((0.0009*cos(q5)*sin(q4) +
 0.0009*cos(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) - ((0.0009*cos(q4)*cos(q5) -
 0.0009*sin(q4)*sin(q5))*(0.000225 +
 0.00045*cos(q3) - 0.0009*cos(q4)*cos(q5) +
 0.0015*sin(q3) + 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/pow((0.003 +
 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))), 2) +
 0.06*sin(q3 + q4 + q5)) +
 dotq4*(0.06*cos(q3 + q4) + ((-0.06*cos(q3 + q4 + q5) +
 0.06*sin(q3 + q4))*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) + ((0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5)))*(-0.06*cos(q3 + q4) - 0.06*sin(q3 + q4 + q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))) +
 0.06*sin(q3 + q4 +q5) - (0.06*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(-0.25 - 0.5*cos(q3) +
 1.0*cos(q4)*cos(q5) - 1.66667*sin(q3) - 1.0*sin(q4) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/pow((3.33333 + 
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)) + (0.1*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(-1.2*cos(q3) - 2.0*cos(q3 + q4) +
 1.0*cos(2*q3 + q4) + 0.15*cos(q3 + q4 + q5) +
 0.3*cos(2*q3 + q4 + q5) - 0.15*sin(q3 + q4) -
 0.3*sin(2*q3 + q4) + 0.6*sin(q3 - q5) -
 0.6*sin(q3 + q5) - 2.0*sin(q3 + q4 + q5) +
 1.0*sin(2*q3 + q4 + q5)))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2))) + ((-1.0*cos(q3 + q4) -
 1.0*sin(q3 + q4 + q5))*(0.000375*dotq1*cos(q3 + q4) -
 0.000075*dotq2*cos(q3 + q4 + q5) +
 0.000375*dotq2*sin(q3 + q4) +
 dotq5*(9.*pow(10, -6)*cos(q3 - q5) + 9.*pow(10, -6)*cos(q3 + q5) +
 0.000015*cos(q4 + q5) + 0.000015*cos(q3 + q4 + q5) +
 4.5*pow(10, -6)*sin(q4 + q5) + 2.25*pow(10, -6)*sin(q3 + q4 + q5)) +
 dotq4*(0.0000225*cos(q4) + 0.00001125*cos(q3 + q4) +
 0.000015*cos(q4 + q5) + 0.000015*cos(q3 + q4 + q5) -
 0.000075*sin(q4) - 0.000075*sin(q3 + q4) +
 4.5*pow(10, -6)*sin(q4 + q5) + 2.25*pow(10, -6)*sin(q3 + q4 + q5)) +
 5.625*pow(10, -6)*dotq3*(1.0*cos(q3 + q4) - 0.8*cos(q3 - q5) +
 0.8*cos(q3 + q5) + 1.33333*cos(q3 + q4 + q5) -
 4.0*sin(q3) - 6.66667*sin(q3 + q4) +
 0.2*sin(q3 + q4 + q5)) +
 0.000075*dotq1*sin(q3 + q4 +q5) - (0.0007335*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(1.0*dotq3*cos(q3) + (-0.0204499*dotq4 - 0.0204499*dotq5)*cos(
 q4 + q5) + 0.0102249*dotq4*cos(q3 + q4 + q5) +
 0.0102249*dotq5*cos(q3 + q4 + q5) +
 3.33333*dotq3*sin(q3) + 0.102249*dotq4*sin(q4) -
 0.0511247*dotq4*sin(q3 + q4))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + (-0.06*cos(q3 + q4 + q5) +
 0.06*sin(q3 + q4) + ((0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))))*(-0.02445*dotq3*cos(q3) + 0.000375*dotq4*cos(q3 + q4) +
 0.007335*dotq3*sin(q3) - 0.00015*dotq4*cos(q5)*sin(q4) -
 0.00015*dotq5*cos(q5)*sin(q4) +
 cos(q4)*(-0.00075*dotq4 - 0.00015*dotq4*sin(q5) -
 0.00015*dotq5*sin(q5)) +
 0.000075*dotq4*sin(q3 + q4 + q5) +
 0.000075*dotq5*sin(q3 + q4 + q5)) + (5.625*pow(10, -6)*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5))*((-4346.67*dotq1 - 1304.0*dotq2 +
 5.33333*dotq3)*cos(q3) + (4.0*dotq3 - 8.0*dotq4)*cos(q3 - q4) - 4.0*dotq4*cos(q4) + 1.0*dotq4*cos(q3 + q4) -
 0.8*dotq4*cos(q3 - q5) +
 2.66667*dotq3*cos(q3 - q4 - q5) -
 5.33333*dotq4*cos(q3 - q4 - q5) -
 5.33333*dotq5*cos(q3 - q4 - q5) - 3.2*dotq5*cos(q5) +
 0.8*dotq4*cos(q3 + q5) + 0.8*dotq5*cos(q3 + q5) -
 5.33333*dotq4*cos(q4 + q5) - 
 5.33333*dotq5*cos(q4 + q5) +
 1.33333*dotq4*cos(q3 + q4 + q5) +
 1.33333*dotq5*cos(q3 + q4 + q5) + 1304.0*dotq1*sin(q3) -
 4346.67*dotq2*sin(q3) - 220.569*dotq3*sin(q3) -
 4.0*dotq4*sin(q3) - 0.4*dotq5*sin(q3) +
 13.3333*dotq3*sin(q3 - q4) -
 26.6667*dotq4*sin(q3 - q4) + 26.6667*dotq4*sin(q4) -
 6.66667*dotq4*sin(q3 + q4) -
 0.8*dotq3*sin(q3 - q4 - q5) +
 1.6*dotq4*sin(q3 - q4 - q5) +
 1.6*dotq5*sin(q3 - q4 - q5) - 0.8*dotq4*sin(q4 + q5) -
 0.8*dotq5*sin(q4 + q5) + 0.2*dotq4*sin(q3 + q4 + q5) +
 0.2*dotq5*sin(q3 + q4 + q5)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))) + (0.0000149625*(8.0*(1.0*cos(q3) +
 1.66667*cos(q3 + q4) - 0.375*cos(q3 + q4 + q5) +
 0.375*sin(q3 + q4) - 0.5*sin(q3 - q5) +
 0.5*sin(q3 + q5) +
 1.66667*sin(q3 + q4 + q5))*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5)) +
 dotq4*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(-13.3333*cos(
 q3 + q4) + 3.0*cos(q3 + q4 + q5) -
 3.0*sin(q3 + q4) - 13.3333*sin(q3 + q4 + q5)) -
 3.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(1.0*cos(q3 + q4) -
 1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
 4.44444*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
 4.44444*sin(q3 + q4) + 1.0*sin(q3 + q4 + q5))) +
 dotq5*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 4.0*sin(q3 - q5) - 4.0*sin(q3 + q5) -
 13.3333*sin(q3 + q4 + q5)) -
 3.0*(1.0*cos(q4)*cos(q5) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) -
 1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
 4.44444*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
 4.44444*sin(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((0.00728604 -
 0.0024814*cos(q3) + 0.00015*cos(q3 - q4) +
 0.00015*cos(q4) - 9.*pow(10, -6)*cos(q3 - q4 - q5) -
 4.5*pow(10, -6)*cos(q4 + q5) - 0.00006*sin(q3) -
 0.000045*sin(q3 - q4) + 0.0000225*sin(q4) -
 0.00003*sin(q3 - q4 - q5) + 0.000018*sin(q5) +
 0.00003*sin(q4 + q5))*(1.0*dotq4*cos(q3 + q4)*cos(
 q3 + q4 + q5) + 2.0*dotq4*sin(q3) + 1.0*dotq5*sin(q3) +
 3.33333*dotq4*cos(q4)*sin(q3) +
 3.33333*dotq4*cos(q5)*sin(q3)*sin(q4) +
 3.33333*dotq5*cos(q5)*sin(q3)*sin(q4) -
 0.5*dotq4*sin(2*(q3 + q4)) + 2.0*dotq4*sin(q3)*sin(q5) +
 1.0*dotq5*sin(q3)*sin(q5) +
 3.33333*dotq4*cos(q4)*sin(q3)*sin(q5) +
 3.33333*dotq5*cos(q4)*sin(q3)*sin(q5) +
 cos(q3)*((-3.33333*dotq4 - 3.33333*dotq5)*cos(q4)*cos(
 q5) + sin(q4)*(3.33333*dotq4 + 3.33333*dotq4*sin(q5) +
 3.33333*dotq5*sin(q5))) -
 1.0*dotq4*sin(q3 + q4)*sin(q3 + q4 + q5) -
 1.0*dotq5*sin(q3 + q4)*sin(q3 + q4 + q5) +
 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((-0.0113 +
 0.0489*cos(q3) - 0.01467*sin(q3) +
 0.00015*cos(q5)*sin(q4) +
 cos(q4)*(0.00075 + 0.00015*sin(q5)))*(-0.1*dotq5*(-1.0 +
 1.0*cos(q3) -
 0.3*sin(q3))*(sin(q3)*(1.0 +
 3.33333*cos(q5)*sin(q4) + (1.0 + 3.33333*cos(q4))*sin(q5)) +
 cos(q3)*(-3.33333*cos(q4)*cos(q5) +
 3.33333*sin(q4)*sin(q5))) -
 0.1*dotq4*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(sin(
 q3)*(2.0 + 3.33333*cos(q5)*sin(q4) + 2.0*sin(q5) +
 cos(q4)*(3.33333 + 3.33333*sin(q5))) +
 cos(q3)*(-3.33333*cos(q4)*cos(q5) +
 sin(q4)*(3.33333 + 3.33333*sin(q5)))) +
 0.03*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(1.0*cos(2*q3 + q4) +
 3.33333*cos(q3 + q4 + q5) -
 3.33333*cos(2*q3 + q4 + q5) - 3.33333*sin(q3 + q4) +
 3.33333*sin(2*q3 + q4) +
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + (0.00369 -
 0.01467*cos(q3) + 0.00015*cos(q4)*cos(q5) -
 0.0489*sin(q3) - 0.00075*sin(q4) -
 0.00015*sin(q4)*sin(
 q5))*(dotq5*(-((
 0.06*cos(q3 + q4 + q5)*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5)))) + ((0.0009*cos(
 q5)*sin(q4) +
 0.0009*cos(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) - ((0.0009*cos(q4)*cos(q5) -
 0.0009*sin(q4)*sin(q5))*(0.000225 +
 0.00045*cos(q3) - 0.0009*cos(q4)*cos(q5) +
 0.0015*sin(q3) + 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/pow((0.003 +
 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))), 2) +
 0.06*sin(q3 + q4 + q5)) +
 dotq4*(0.06*cos(
 q3 + q4) + ((-0.06*cos(q3 + q4 + q5) +
 0.06*sin(q3 + q4))*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) + 
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) + ((0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5)))*(-0.06*cos(
 q3 + q4) - 0.06*sin(q3 + q4 + q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))) +
 0.06*sin(q3 + q4 + q5) - (0.06*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(-0.25 - 0.5*cos(q3) +
 1.0*cos(q4)*cos(q5) - 1.66667*sin(q3) -
 1.0*sin(q4) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)) + (0.1*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(-1.2*cos(q3) - 2.0*cos(q3 + q4) +
 1.0*cos(2*q3 + q4) + 0.15*cos(q3 + q4 + q5) +
 0.3*cos(2*q3 + q4 + q5) - 0.15*sin(q3 + q4) -
 0.3*sin(2*q3 + q4) + 0.6*sin(q3 - q5) -
 0.6*sin(q3 + q5) - 2.0*sin(q3 + q4 + q5) +
 1.0*sin(2*q3 + q4 + q5)))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2))))/(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)));
	gsl_matrix_set(C2, 0, 0, value);

//C2(0,1)
	value = 9.*pow(10, -6)*dotq5*cos(q5) - 0.000075*dotq2*cos(q3 + q4 + q5) +
 dotq3*(-7.5*pow(10, -6)*cos(q4 + q5) - 7.5*pow(10, -6)*cos(q3 + q4 + q5) -
 2.25*pow(10, -6)*sin(q4 + q5) - 1.125*pow(10, -6)*sin(q3 + q4 + q5)) +
 0.000075*dotq1*sin(
 q3 + q4 +
 q5) - (0.000045*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(0.166667*dotq3*cos(
 q4 + q5) + (-0.333333*dotq3 - 0.166667*dotq4 -
 0.166667*dotq5)*cos(q3 + q4 + q5) -
 0.833333*dotq3*sin(q4) + 1.66667*dotq3*sin(q3 + q4) +
 0.833333*dotq4*sin(q3 + q4))*sin(q3 + q4 + q5))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + (0.000075*sin(
 q3 + q4 + q5)*(0.3*dotq3*cos(q3 - q4) -
 0.15*dotq3*cos(q3 + q4) - 0.075*dotq4*cos(q3 + q4) +
 0.12*dotq3*cos(q3 - q5) - 0.12*dotq5*cos(q3 - q5) +
 0.2*dotq3*cos(q3 - q4 - q5) - 0.12*dotq3*cos(q3 + q5) -
 0.12*dotq5*cos(q3 + q5) + 0.2*dotq3*cos(q4 + q5) -
 0.1*dotq4*cos(q4 + q5) - 0.1*dotq5*cos(q4 + q5) -
 0.2*dotq3*cos(q3 + q4 + q5) - 0.1*dotq4*cos(q3 + q4 + q5) -
 0.1*dotq5*cos(q3 + q4 + q5) + 0.6*dotq3*sin(q3) +
 1.0*dotq3*sin(q3 - q4) - 5.0*dotq2*sin(q4) - 1.0*dotq3*sin(q4) +
 0.5*dotq4*sin(q4) - 1.0*dotq1*cos(q5)*sin(q4) +
 1.0*dotq3*sin(q3 + q4) + 0.5*dotq4*sin(q3 + q4) -
 0.06*dotq3*sin(q3 - q4 - q5) - 1.0*dotq2*sin(q4)*sin(q5) +
 cos(q4)*(-5.0*dotq1 + 0.15*dotq3 - 0.15*dotq4 +
 1.0*dotq2*cos(q5) - 1.0*dotq1*sin(q5)) +
 0.03*dotq3*sin(q4 + q5) - 0.03*dotq4*sin(q4 + q5) -
 0.03*dotq5*sin(q4 + q5) - 0.03*dotq3*sin(q3 + q4 + q5) -
 0.015*dotq4*sin(q3 + q4 + q5) -
 0.015*dotq5*sin(q3 + q4 + q5)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))) + (4.5*pow(10, -6)*((-10.0*dotq3 - 5.0*dotq4)*cos(
 q3 + q4) + 1.0*dotq3*cos(q5)*sin(q4) +
 dotq3*cos(q4)*(5.0 + 1.0*sin(q5)) -
 2.0*dotq3*sin(q3 + q4 + q5) - 1.0*dotq4*sin(q3 + q4 + q5) -
 1.0*dotq5*sin(q3 + q4 + q5))*(cos(
 q3 + q4 + q5)*(-3.33333 - 1.0*cos(q5)*sin(q4) +
 cos(q4)*(-1.0 - 1.0*sin(q5))) + (-0.25 - 0.5*cos(q3) +
 1.0*cos(q4)*cos(q5) - 1.66667*sin(q3) - 1.0*sin(q4) -
 1.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + ((0.0000372812 +
 0.000045*cos(q3) + 0.000075*cos(q4) + 0.000075*cos(q3 + q4) -
 4.5*pow(10, -6)*cos(q4 + q5) - 2.25*pow(10, -6)*cos(q3 + q4 + q5) +
 0.0000225*sin(q4) + 0.00001125*sin(q3 + q4) -
 9.*pow(10, -6)*sin(q3 - q5) + 9.*pow(10, -6)*sin(q3 + q5) +
 0.000015*sin(q4 + q5) +
 0.000015*sin(q3 + q4 + q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
 q3 + q5) + (-3.33333*dotq4 - 3.33333*dotq5 +
 1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
 0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((3.0*cos(q3 + q4) -
 4.0*cos(q3 - q5) + 4.0*cos(q3 + q5) +
 13.3333*cos(q3 + q4 + q5) - 8.0*sin(q3) -
 13.3333*sin(q3 + q4) +
 3.0*sin(q3 + q4 +
 q5))*(0.0000149625*(dotq4*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 13.3333*sin(q3 + q4 + q5)) -
 4.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3) +
 0.75*sin(q3 + q4 + q5))) +
 dotq5*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 4.0*sin(q3 + q5) - 13.3333*sin(q3 + q4 + q5)) -
 4.0*(1.0*cos(q4)*cos(q5) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3) +
 0.75*sin(q3 + q4 + q5))) +
 4.0*(1.0*cos(q3) - 0.75*cos(q3 + q4 + q5) +
 1.0*sin(q3 + q5) +
 3.33333*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))) +
 0.0000149625*((-1.0*dotq4 - 0.5*dotq5)*cos(
 q3 + q5) + (-3.33333*dotq4 - 3.33333*dotq5 +
 1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
 0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5)))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 3) + (-0.06*cos(q3 + q4 + q5) +
 0.06*sin(
 q3 + q4) + ((0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) + 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))))*(-0.00015*dotq3*sin(
 q3 + q4 + q5) - 0.00015*dotq4*sin(q3 + q4 + q5) -
 0.00015*dotq5*sin(
 q3 + q4 +
 q5) + (0.0489*(1.0*dotq3*cos(q3) - 0.3*dotq3*sin(q3) +
 0.00306748*dotq4*cos(q5)*sin(q4) +
 0.00306748*dotq5*cos(q5)*sin(q4) +
 cos(q4)*(0.0153374*dotq4 + 0.00306748*dotq4*sin(q5) +
 0.00306748*dotq5*sin(q5)))*sin(
 q3 + q4 + q5))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + ((0.00369 - 0.01467*cos(q3) +
 0.00015*cos(q4)*cos(q5) - 0.0489*sin(q3) -
 0.00075*sin(q4) -
 0.00015*sin(q4)*sin(q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
 q3 + q5) + (-3.33333*dotq4 - 3.33333*dotq5 +
 1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
 0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))), 2) - (0.01467*(1.0*dotq5*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) -
 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) + 2.0*sin(q4) +
 2.0*sin(q4)*sin(q5))*(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) +
 2.0*(1.0*cos(q4)*cos(q5) - 1.0*sin(q4)*sin(q5))*(-0.25 -
 0.5*cos(q3) + 1.0*cos(q4)*cos(q5) - 1.66667*sin(q3) -
 1.0*sin(q4) - 1.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5) +
 1.0*(2.0*cos(q5)*sin(q4) + 2.0*cos(q4)*sin(q5))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)))*sin(
 q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +
 1.0*dotq4*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) -
 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) + 2.0*sin(q4) +
 2.0*sin(q4)*sin(q5))*(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) -
 2.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(0.25 + 0.5*cos(q3) -
 1.0*cos(q4)*cos(q5) + 1.66667*sin(q3) + 1.0*sin(q4) +
 1.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5) +
 2.0*(1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)))*sin(
 q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +
 3.33333*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(0.6*cos(q3) -
 0.15*cos(q3 + q4 + q5) - 0.3*cos(2*q3 + q4 + q5) +
 0.6*sin(q3 + q5) + 2.0*sin(q3 + q4 + q5) -
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 + 
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)) + ((0.0015 - 0.0015*cos(q3) +
 0.00045*sin(q3))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5))*(0.00015*dotq3*cos(q3 + q4 + q5) +
 0.00015*dotq4*cos(q3 + q4 + q5) +
 0.00015*dotq5*cos(
 q3 + q4 +
 q5) + (0.01467*(1.0*dotq3*cos(
 q3) + (-0.0102249*dotq4 - 0.0102249*dotq5)*cos(q4)*cos(
 q5) + 3.33333*dotq3*sin(q3) +
 0.0511247*dotq4*sin(q4) +
 0.0102249*dotq4*sin(q4)*sin(q5) +
 0.0102249*dotq5*sin(q4)*sin(q5))*sin(
 q3 + q4 + q5))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + ((-0.0113 + 0.0489*cos(q3) -
 0.01467*sin(q3) + 0.00015*cos(q5)*sin(q4) +
 cos(q4)*(0.00075 + 0.00015*sin(q5)))*((-1.0*dotq4 -
 0.5*dotq5)*cos(
 q3 + q5) + (-3.33333*dotq4 - 3.33333*dotq5 +
 1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
 0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))), 2) + (0.489*(-0.015*dotq5*(1.0*cos(
 q3 + q5) + 6.66667*cos(q3 + q4 + q5) +
 1.0*cos(q3 + 2*q4 + q5) - 2.0*sin(q3))*(3.33333 -
 3.33333*cos(q3) + 1.0*sin(q3)) -
 0.03*dotq4*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3))*(3.33333 -
 3.33333*cos(q3) + 1.0*sin(q3)) +
 0.1*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(1.0*cos(q3 + q4 + q5) -
 1.0*cos(2*q3 + q4 + q5) +
 0.3*sin(2*q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)))/(0.003 +
 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))) + ((-1.0*cos(q3 + q4) -
 1.0*sin(q3 + q4 + q5))*(-0.000075*dotq2*cos(q3 + q4 + q5) +
 dotq4*(0.000015*cos(q4 + q5) + 0.000015*cos(q3 + q4 + q5) +
 4.5*pow(10, -6)*sin(q4 + q5) + 2.25*pow(10, -6)*sin(q3 + q4 + q5)) +
 dotq5*(9.*pow(10, -6)*cos(q3 + q5) + 0.000015*cos(q4 + q5) +
 0.000015*cos(q3 + q4 + q5) + 4.5*pow(10, -6)*sin(q4 + q5) +
 2.25*pow(10, -6)*sin(q3 + q4 + q5)) +
 4.5*pow(10, -6)*dotq3*(1.0*cos(q3 + q5) +
 1.66667*cos(q3 + q4 + q5) - 0.5*sin(q3) +
 0.25*sin(q3 + q4 + q5)) +
 0.000075*dotq1*sin(
 q3 + q4 +
 q5) - (0.0007335*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(1.0*dotq3*cos(
 q3) + (-0.0204499*dotq4 - 0.0204499*dotq5)*cos(
 q4 + q5) + 0.0102249*dotq4*cos(q3 + q4 + q5) +
 0.0102249*dotq5*cos(q3 + q4 + q5) +
 3.33333*dotq3*sin(q3) + 0.102249*dotq4*sin(q4) -
 0.0511247*dotq4*sin(q3 + q4))*sin(
 q3 + q4 + q5))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + (1.125*pow(10, -6)*sin(
 q3 + q4 +
 q5)*((-21733.3*dotq1 - 6520.0*dotq2 + 26.6667*dotq3)*cos(
 q3) + (20.0*dotq3 - 40.0*dotq4)*cos(q3 - q4) -
 20.0*dotq4*cos(q4) + 5.0*dotq4*cos(q3 + q4) -
 4.0*dotq4*cos(q3 - q5) +
 13.3333*dotq3*cos(q3 - q4 - q5) -
 26.6667*dotq4*cos(q3 - q4 - q5) -
 26.6667*dotq5*cos(q3 - q4 - q5) - 16.0*dotq5*cos(q5) +
 4.0*dotq4*cos(q3 + q5) + 4.0*dotq5*cos(q3 + q5) -
 26.6667*dotq4*cos(q4 + q5) -
 26.6667*dotq5*cos(q4 + q5) +
 6.66667*dotq4*cos(q3 + q4 + q5) +
 6.66667*dotq5*cos(q3 + q4 + q5) + 6520.0*dotq1*sin(q3) -
 21733.3*dotq2*sin(q3) - 1102.84*dotq3*sin(q3) -
 20.0*dotq4*sin(q3) - 2.0*dotq5*sin(q3) +
 66.6667*dotq3*sin(q3 - q4) -
 133.333*dotq4*sin(q3 - q4) + 133.333*dotq4*sin(q4) -
 33.3333*dotq4*sin(q3 + q4) -
 4.0*dotq3*sin(q3 - q4 - q5) +
 8.0*dotq4*sin(q3 - q4 - q5) +
 8.0*dotq5*sin(q3 - q4 - q5) - 4.0*dotq4*sin(q4 + q5) -
 4.0*dotq5*sin(q4 + q5) + 1.0*dotq4*sin(q3 + q4 + q5) +
 1.0*dotq5*sin(q3 + q4 + q5)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))) + (0.0007335*(1.0*dotq3*cos(q3) -
 0.0153374*dotq4*cos(q3 + q4) - 0.3*dotq3*sin(q3) +
 0.00613497*dotq4*cos(q5)*sin(q4) +
 0.00613497*dotq5*cos(q5)*sin(q4) +
 cos(q4)*(0.0306748*dotq4 + 0.00613497*dotq4*sin(q5) +
 0.00613497*dotq5*sin(q5)) -
 0.00306748*dotq4*sin(q3 + q4 + q5) -
 0.00306748*dotq5*sin(q3 + q4 + q5))*(cos(
 q3 + q4 + q5)*(6.66667 + 2.0*cos(q5)*sin(q4) +
 cos(q4)*(2.0 + 2.0*sin(q5))) + (0.5 + 1.0*cos(q3) -
 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) + 2.0*sin(q4) +
 2.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))) + (0.0000149625*(dotq4*((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 13.3333*sin(q3 + q4 + q5)) -
 4.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3) +
 0.75*sin(q3 + q4 + q5))) +
 dotq5*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 4.0*sin(q3 + q5) - 13.3333*sin(q3 + q4 + q5)) -
 4.0*(1.0*cos(q4)*cos(q5) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3) +
 0.75*sin(q3 + q4 + q5))) +
 4.0*(1.0*cos(q3) - 0.75*cos(q3 + q4 + q5) +
 1.0*sin(q3 + q5) +
 3.33333*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))))/pow((3.33333 + 1.0*cos(q5)*sin(q4) + 
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((0.00728604 -
 0.0024814*cos(q3) + 0.00015*cos(q3 - q4) +
 0.00015*cos(q4) - 9.*pow(10, -6)*cos(q3 - q4 - q5) -
 4.5*pow(10, -6)*cos(q4 + q5) - 0.00006*sin(q3) -
 0.000045*sin(q3 - q4) + 0.0000225*sin(q4) -
 0.00003*sin(q3 - q4 - q5) + 0.000018*sin(q5) +
 0.00003*sin(q4 + q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
 q3 + q5) + (-3.33333*dotq4 - 3.33333*dotq5 +
 1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
 0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) - (0.03*(0.00369 -
 0.01467*cos(q3) + 0.00015*cos(q4)*cos(q5) -
 0.0489*sin(q3) - 0.00075*sin(q4) -

 0.00015*sin(q4)*sin(
 q5))*(1.0*dotq5*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) -
 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) +
 2.0*sin(q4) + 2.0*sin(q4)*sin(q5))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5))) +
 2.0*(1.0*cos(q4)*cos(q5) - 1.0*sin(q4)*sin(q5))*(-0.25 -
 0.5*cos(q3) + 1.0*cos(q4)*cos(q5) -
 1.66667*sin(q3) - 1.0*sin(q4) -
 1.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5) +
 1.0*(2.0*cos(q5)*sin(q4) +
 2.0*cos(q4)*sin(q5))*(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*sin(q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +
 1.0*dotq4*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) -
 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) +
 2.0*sin(q4) + 2.0*sin(q4)*sin(q5))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5))) -
 2.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(0.25 + 0.5*cos(q3) -
 1.0*cos(q4)*cos(q5) + 1.66667*sin(q3) +
 1.0*sin(q4) + 1.0*sin(q4)*sin(q5))*sin(
 q3 + q4 + q5) +
 2.0*(1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*sin(q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +
 3.33333*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(0.6*cos(q3) -
 0.15*cos(q3 + q4 + q5) - 0.3*cos(2*q3 + q4 + q5) +
 0.6*sin(q3 + q5) + 2.0*sin(q3 + q4 + q5) -
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((-0.0113 +
 0.0489*cos(q3) - 0.01467*sin(q3) +
 0.00015*cos(q5)*sin(q4) +
 cos(q4)*(0.00075 +
 0.00015*sin(q5)))*(-0.015*dotq5*(1.0*cos(q3 + q5) + 
 6.66667*cos(q3 + q4 + q5) + 1.0*cos(q3 + 2*q4 + q5) -
 2.0*sin(q3))*(3.33333 - 3.33333*cos(q3) +
 1.0*sin(q3)) -
 0.03*dotq4*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3))*(3.33333 -
 3.33333*cos(q3) + 1.0*sin(q3)) +
 0.1*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(1.0*cos(q3 + q4 + q5) -
 1.0*cos(2*q3 + q4 + q5) +
 0.3*sin(2*q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) - (0.03*(0.00015*cos(q3 + q4 + q5) -
 0.00075*sin(
 q3 + q4))*(1.0*dotq5*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4)*cos(q5) +
 3.33333*sin(q3) + 2.0*sin(q4) +
 2.0*sin(q4)*sin(q5))*(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) +
 2.0*(1.0*cos(q4)*cos(q5) - 1.0*sin(q4)*sin(q5))*(-0.25 -
 0.5*cos(q3) + 1.0*cos(q4)*cos(q5) - 1.66667*sin(q3) -
 1.0*sin(q4) - 1.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5) +
 1.0*(2.0*cos(q5)*sin(q4) + 2.0*cos(q4)*sin(q5))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)))*sin(
 q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +
 1.0*dotq4*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4)*cos(q5) +
 3.33333*sin(q3) + 2.0*sin(q4) +
 2.0*sin(q4)*sin(q5))*(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) -
 2.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(0.25 + 0.5*cos(q3) -
 1.0*cos(q4)*cos(q5) + 1.66667*sin(q3) + 1.0*sin(q4) +
 1.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5) +
 2.0*(1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)))*sin(
 q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +
 3.33333*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(0.6*cos(q3) - 0.15*cos(q3 + q4 + q5) -
 0.3*cos(2*q3 + q4 + q5) + 0.6*sin(q3 + q5) +
 2.0*sin(q3 + q4 + q5) -
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((0.00045*cos(q3)*cos(q4) +
 0.0003*cos(q3 + q4) - 0.00045*sin(q3)*sin(q4) +
 0.00015*sin(
 q3 + q4 + q5))*(-0.015*dotq5*(1.0*cos(q3 + q5) +
 6.66667*cos(q3 + q4 + q5) + 1.0*cos(q3 + 2*q4 + q5) -
 2.0*sin(q3))*(3.33333 - 3.33333*cos(q3) + 1.0*sin(q3)) -
 0.03*dotq4*(1.0*cos(q3 + q5) + 3.33333*cos(q3 + q4 + q5) -
 1.0*sin(q3))*(3.33333 - 3.33333*cos(q3) + 1.0*sin(q3)) +
 0.1*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(1.0*cos(q3 + q4 + q5) -
 1.0*cos(2*q3 + q4 + q5) +
 0.3*sin(2*q3 + q4 + q5))))/pow((3.33333 + 1.0*cos(q5)*sin(q4) + 
 cos(q4)*(1.0 + 1.0*sin(q5))), 2);
 	gsl_matrix_set(C2, 0, 1, value);

//C2(1,0)
	value = -9.*pow(10, -6)*dotq4*cos(q5) +
 4.5*pow(10, -6)*dotq5*cos(q5) - 0.000075*dotq2*cos(q3 + q4 + q5) -
 4.5*pow(10, -6)*dotq3*(1.0*cos(q3 - q5) + 1.0*cos(q3 + q5) +
 1.66667*cos(q4 + q5) + 1.66667*cos(q3 + q4 + q5) +
 0.5*sin(q4 + q5) + 0.25*sin(q3 + q4 + q5)) +
 0.000075*dotq1*sin(
 q3 + q4 +
 q5) - (7.5*pow(10, -6)*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(1.0*dotq3*cos(q4)*cos(
 q5) + (-2.0*dotq3 - 1.0*dotq4 - 1.0*dotq5)*cos(q3 + q4 + q5) -
 1.0*dotq3*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + (-0.06*cos(q3 + q4 + q5) +
 0.06*sin(
 q3 + q4) + ((0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) + 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))))*(0.000075*dotq3*cos(
 q5)*sin(q4) +
 0.000075*dotq3*cos(q4)*sin(
 q5) + (-0.00015*dotq3 - 0.000075*dotq4 - 0.000075*dotq5)*sin(
 q3 + q4 + q5)) + (0.000075*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5))*(0.06*dotq4*cos(q3 - q5) +
 0.2*dotq3*cos(q3 - q4 - q5) + 0.12*dotq3*cos(q5) +
 1.0*dotq2*cos(q4)*cos(q5) - 0.12*dotq3*cos(q3 + q5) +
 0.06*dotq4*cos(q3 + q5) - 0.06*dotq5*cos(q3 + q5) +
 0.2*dotq3*cos(q4 + q5) - 0.1*dotq4*cos(q4 + q5) -
 0.1*dotq5*cos(q4 + q5) - 0.2*dotq3*cos(q3 + q4 + q5) -
 0.1*dotq4*cos(q3 + q4 + q5) - 0.1*dotq5*cos(q3 + q4 + q5) +
 0.06*dotq3*sin(q3) - 1.0*dotq1*cos(q5)*sin(q4) -
 0.06*dotq3*sin(q3 - q4 - q5) - 1.0*dotq1*cos(q4)*sin(q5) -
 1.0*dotq2*sin(q4)*sin(q5) + 0.03*dotq3*sin(q4 + q5) -
 0.03*dotq4*sin(q4 + q5) - 0.03*dotq5*sin(q4 + q5) -
 0.03*dotq3*sin(q3 + q4 + q5) -
 0.015*dotq4*sin(q3 + q4 + q5) -
 0.015*dotq5*sin(q3 + q4 + q5)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + ((4.0*cos(q3 + q5) +
 13.3333*cos(q3 + q4 + q5) - 4.0*sin(q3) +
 3.0*sin(q3 + q4 + q5))*(0.0000149625*dotq4*cos(q3 + q4)*cos(
 q3 + q4 + q5) + 0.000029925*dotq4*sin(q3) +
 0.0000149625*dotq5*sin(q3) +
 0.000049875*dotq4*cos(q4)*sin(q3) +
 0.000049875*dotq4*cos(q5)*sin(q3)*sin(q4) +
 0.000049875*dotq5*cos(q5)*sin(q3)*sin(q4) -
 7.48125*pow(10, -6)*dotq4*sin(2*(q3 + q4)) +
 0.000029925*dotq4*sin(q3)*sin(q5) +
 0.0000149625*dotq5*sin(q3)*sin(q5) +
 0.000049875*dotq4*cos(q4)*sin(q3)*sin(q5) +
 0.000049875*dotq5*cos(q4)*sin(q3)*sin(q5) +
 0.0000149625*cos(
 q3)*((-3.33333*dotq4 - 3.33333*dotq5)*cos(q4)*cos(q5) +
 sin(q4)*(3.33333*dotq4 + (3.33333*dotq4 +
 3.33333*dotq5)*sin(q5))) -
 0.0000149625*dotq4*sin(q3 + q4)*sin(q3 + q4 + q5) -
 0.0000149625*dotq5*sin(q3 + q4)*sin(q3 + q4 + q5) + 
 0.0001197*(1.0*cos(q3) + 1.66667*cos(q3 + q4) -
 0.375*cos(q3 + q4 + q5) + 0.375*sin(q3 + q4) -
 0.5*sin(q3 - q5) + 0.5*sin(q3 + q5) +
 1.66667*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5)) +
 0.0000149625*dotq4*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(-13.3333*cos(q3 + q4) +
 3.0*cos(q3 + q4 + q5) - 3.0*sin(q3 + q4) -
 13.3333*sin(q3 + q4 + q5)) -
 3.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(1.0*cos(q3 + q4) -
 1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
 4.44444*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
 4.44444*sin(q3 + q4) + 1.0*sin(q3 + q4 + q5))) +
 0.0000149625*dotq5*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 4.0*sin(q3 - q5) - 4.0*sin(q3 + q5) -
 13.3333*sin(q3 + q4 + q5)) -
 3.0*(1.0*cos(q4)*cos(q5) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) -
 1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
 4.44444*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
 4.44444*sin(q3 + q4) + 1.0*sin(q3 + q4 + q5))) +
 7.48125*pow(10, -6)*dotq4*sin(2*(q3 + q4 + q5)) +
 7.48125*pow(10, -6)*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 3) + ((0.0000186406 +
 4.5*pow(10, -6)*cos(q3) - 4.5*pow(10, -6)*cos(q4 + q5) -
 2.25*pow(10, -6)*cos(q3 + q4 + q5) + 9.*pow(10, -6)*sin(q3 + q5) +
 0.000015*sin(q4 + q5) +
 0.000015*sin(q3 + q4 + q5))*(1.0*dotq4*cos(q3 + q4)*cos(
 q3 + q4 + q5) + 2.0*dotq4*sin(q3) + 1.0*dotq5*sin(q3) +
 3.33333*dotq4*cos(q4)*sin(q3) +
 3.33333*dotq4*cos(q5)*sin(q3)*sin(q4) +
 3.33333*dotq5*cos(q5)*sin(q3)*sin(q4) -
 0.5*dotq4*sin(2*(q3 + q4)) + 2.0*dotq4*sin(q3)*sin(q5) +
 1.0*dotq5*sin(q3)*sin(q5) +
 3.33333*dotq4*cos(q4)*sin(q3)*sin(q5) +
 3.33333*dotq5*cos(q4)*sin(q3)*sin(q5) +
 cos(q3)*((-3.33333*dotq4 - 3.33333*dotq5)*cos(q4)*cos(q5) +

 sin(q4)*(3.33333*dotq4 + 3.33333*dotq4*sin(q5) +
 3.33333*dotq5*sin(q5))) -
 1.0*dotq4*sin(q3 + q4)*sin(q3 + q4 + q5) -
 1.0*dotq5*sin(q3 + q4)*sin(q3 + q4 + q5) +
 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + (0.00015*sin(
 q3 + q4 +
 q5)*(-0.1*dotq5*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(sin(
 q3)*(1.0 +
 3.33333*cos(q5)*sin(q4) + (1.0 + 3.33333*cos(q4))*sin(
 q5)) + cos(
 q3)*(-3.33333*cos(q4)*cos(q5) +
 3.33333*sin(q4)*sin(q5))) -
 0.1*dotq4*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(sin(
 q3)*(2.0 + 3.33333*cos(q5)*sin(q4) + 2.0*sin(q5) + 
 cos(q4)*(3.33333 + 3.33333*sin(q5))) +
 cos(q3)*(-3.33333*cos(q4)*cos(q5) +
 sin(q4)*(3.33333 + 3.33333*sin(q5)))) +
 0.03*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(1.0*cos(2*q3 + q4) +
 3.33333*cos(q3 + q4 + q5) - 3.33333*cos(2*q3 + q4 + q5) -
 3.33333*sin(q3 + q4) + 3.33333*sin(2*q3 + q4) +
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) +
 0.00015*cos(
 q3 + q4 +
 q5)*(dotq5*(-((
 0.06*cos(
 q3 + q4 + q5)*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5)))) + ((0.0009*cos(q5)*sin(
 q4) + 0.0009*cos(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) - ((0.0009*cos(q4)*cos(q5) -
 0.0009*sin(q4)*sin(q5))*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/pow((0.003 +
 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))), 2) +
 0.06*sin(q3 + q4 + q5)) +
 dotq4*(0.06*cos(
 q3 + q4) + ((-0.06*cos(q3 + q4 + q5) +
 0.06*sin(q3 + q4))*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) + ((0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5)))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))) +
 0.06*sin(
 q3 + q4 +
 q5) - (0.06*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(-0.25 - 0.5*cos(q3) +
 1.0*cos(q4)*cos(q5) - 1.66667*sin(q3) - 1.0*sin(q4) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)) + (0.1*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(-1.2*cos(q3) - 2.0*cos(q3 + q4) +
 1.0*cos(2*q3 + q4) + 0.15*cos(q3 + q4 + q5) +
 0.3*cos(2*q3 + q4 + q5) - 0.15*sin(q3 + q4) -
 0.3*sin(2*q3 + q4) + 0.6*sin(q3 - q5) - 0.6*sin(q3 + q5) -
 2.0*sin(q3 + q4 + q5) + 1.0*sin(2*q3 + q4 + q5)))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)) + (0.1*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*sin(
 q3 + q4 + q5)*(0.00015*dotq3*cos(q3 + q4 + q5) +
 0.00015*dotq4*cos(q3 + q4 + q5) +
 0.00015*dotq5*cos(q3 + q4 + q5) -
 0.00075*dotq3*sin(q3 + q4) -
 0.00075*dotq4*sin(
 q3 + q4) + (0.01467*(1.0*dotq3*cos(
 q3) + (-0.0102249*dotq4 - 0.0102249*dotq5)*cos(q4)*cos(
 q5) + 3.33333*dotq3*sin(q3) +
 0.0511247*dotq4*sin(q4) +
 0.0102249*dotq4*sin(q4)*sin(q5) +
 0.0102249*dotq5*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + ((-0.0113 + 0.0489*cos(q3) -
 0.01467*sin(q3) + 0.00015*cos(q5)*sin(q4) +
 cos(q4)*(0.00075 + 0.00015*sin(q5)))*(1.0*dotq4*cos(
 q3 + q4)*cos(q3 + q4 + q5) + 2.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 3.33333*dotq4*cos(q4)*sin(q3) +
 3.33333*dotq4*cos(q5)*sin(q3)*sin(q4) +
 3.33333*dotq5*cos(q5)*sin(q3)*sin(q4) -
 0.5*dotq4*sin(2*(q3 + q4)) + 2.0*dotq4*sin(q3)*sin(q5) +
 1.0*dotq5*sin(q3)*sin(q5) +
 3.33333*dotq4*cos(q4)*sin(q3)*sin(q5) +
 3.33333*dotq5*cos(q4)*sin(q3)*sin(q5) +
 cos(q3)*((-3.33333*dotq4 - 3.33333*dotq5)*cos(q4)*cos(
 q5) + sin(
 q4)*(3.33333*dotq4 + 3.33333*dotq4*sin(q5) +
 3.33333*dotq5*sin(q5))) -
 1.0*dotq4*sin(q3 + q4)*sin(q3 + q4 + q5) -
 1.0*dotq5*sin(q3 + q4)*sin(q3 + q4 + q5) +
 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))), 2) + (0.489*(-0.1*dotq5*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(sin(
 q3)*(1.0 +
 3.33333*cos(q5)*sin(
 q4) + (1.0 + 3.33333*cos(q4))*sin(q5)) +
 cos(q3)*(-3.33333*cos(q4)*cos(q5) +
 3.33333*sin(q4)*sin(q5))) -
 0.1*dotq4*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(sin(
 q3)*(2.0 + 3.33333*cos(q5)*sin(q4) + 2.0*sin(q5) +
 cos(q4)*(3.33333 + 3.33333*sin(q5))) +
 cos(q3)*(-3.33333*cos(q4)*cos(q5) +
 sin(q4)*(3.33333 + 3.33333*sin(q5)))) +
 0.03*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(1.0*cos(2*q3 + q4) +
 3.33333*cos(q3 + q4 + q5) -
 3.33333*cos(2*q3 + q4 + q5) - 3.33333*sin(q3 + q4) +
 3.33333*sin(2*q3 + q4) +
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + (-0.06*cos(q3 + q4 + q5) - (
 0.03*(0.5 + 1.0*cos(q3) - 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) +
 2.0*sin(q4) + 2.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5))/(
 3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))))*(dotq3*(-0.00075*cos(q3 + q4) - 
 0.00015*sin(q3 + q4 + q5)) +
 dotq4*(-0.00075*cos(q3 + q4) - 0.00015*sin(q3 + q4 + q5)) -
 0.00015*dotq5*sin(
 q3 + q4 +
 q5) + (0.0489*(1.0*dotq3*cos(q3) - 0.3*dotq3*sin(q3) +
 0.00306748*dotq4*cos(q5)*sin(q4) +
 0.00306748*dotq5*cos(q5)*sin(q4) +
 cos(q4)*(0.0153374*dotq4 + 0.00306748*dotq4*sin(q5) +
 0.00306748*dotq5*sin(q5)))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + ((0.00369 - 0.01467*cos(q3) +
 0.00015*cos(q4)*cos(q5) - 0.0489*sin(q3) -
 0.00075*sin(q4) -
 0.00015*sin(q4)*sin(q5))*(1.0*dotq4*cos(q3 + q4)*cos(
 q3 + q4 + q5) + 2.0*dotq4*sin(q3) + 1.0*dotq5*sin(q3) +
 3.33333*dotq4*cos(q4)*sin(q3) +
 3.33333*dotq4*cos(q5)*sin(q3)*sin(q4) +
 3.33333*dotq5*cos(q5)*sin(q3)*sin(q4) -
 0.5*dotq4*sin(2*(q3 + q4)) + 2.0*dotq4*sin(q3)*sin(q5) +
 1.0*dotq5*sin(q3)*sin(q5) +
 3.33333*dotq4*cos(q4)*sin(q3)*sin(q5) +
 3.33333*dotq5*cos(q4)*sin(q3)*sin(q5) +
 cos(q3)*((-3.33333*dotq4 - 3.33333*dotq5)*cos(q4)*cos(q5) +
 sin(q4)*(3.33333*dotq4 + 3.33333*dotq4*sin(q5) +
 3.33333*dotq5*sin(q5))) -
 1.0*dotq4*sin(q3 + q4)*sin(q3 + q4 + q5) -
 1.0*dotq5*sin(q3 + q4)*sin(q3 + q4 + q5) +
 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5))), 2) +
 0.489*(dotq5*(-((
 0.06*cos(
 q3 + q4 + q5)*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5)))) + ((0.0009*cos(
 q5)*sin(q4) +
 0.0009*cos(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) - ((0.0009*cos(q4)*cos(q5) -
 0.0009*sin(q4)*sin(q5))*(0.000225 +
 0.00045*cos(q3) - 0.0009*cos(q4)*cos(q5) +
 0.0015*sin(q3) + 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/pow((0.003 +
 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))), 2) +
 0.06*sin(q3 + q4 + q5)) +
 dotq4*(0.06*cos(
 q3 + q4) + ((-0.06*cos(q3 + q4 + q5) +
 0.06*sin(q3 + q4))*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) + ((0.0009*cos(q5)*sin(q4) +
 
 cos(q4)*(0.0009 + 0.0009*sin(q5)))*(-0.06*cos(
 q3 + q4) - 0.06*sin(q3 + q4 + q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))) +
 0.06*sin(
 q3 + q4 +
 q5) - (0.06*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(-0.25 - 0.5*cos(q3) +
 1.0*cos(q4)*cos(q5) - 1.66667*sin(q3) - 1.0*sin(q4) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)) + (0.1*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(-1.2*cos(q3) - 2.0*cos(q3 + q4) +
 1.0*cos(2*q3 + q4) + 0.15*cos(q3 + q4 + q5) +
 0.3*cos(2*q3 + q4 + q5) - 0.15*sin(q3 + q4) -
 0.3*sin(2*q3 + q4) + 0.6*sin(q3 - q5) -
 0.6*sin(q3 + q5) - 2.0*sin(q3 + q4 + q5) +
 1.0*sin(2*q3 + q4 + q5)))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2))) - (1.0*sin(
 q3 + q4 + q5)*(0.000375*dotq1*cos(q3 + q4) -
 0.000075*dotq2*cos(q3 + q4 + q5) +
 0.000375*dotq2*sin(q3 + q4) +
 dotq5*(9.*pow(10, -6)*cos(q3 - q5) + 9.*pow(10, -6)*cos(q3 + q5) +
 0.000015*cos(q4 + q5) + 0.000015*cos(q3 + q4 + q5) +
 4.5*pow(10, -6)*sin(q4 + q5) + 2.25*pow(10, -6)*sin(q3 + q4 + q5)) +
 dotq4*(0.0000225*cos(q4) + 0.00001125*cos(q3 + q4) +
 0.000015*cos(q4 + q5) + 0.000015*cos(q3 + q4 + q5) -
 0.000075*sin(q4) - 0.000075*sin(q3 + q4) +
 4.5*pow(10, -6)*sin(q4 + q5) + 2.25*pow(10, -6)*sin(q3 + q4 + q5)) +
 5.625*pow(10, -6)*dotq3*(1.0*cos(q3 + q4) - 0.8*cos(q3 - q5) +
 0.8*cos(q3 + q5) + 1.33333*cos(q3 + q4 + q5) -
 4.0*sin(q3) - 6.66667*sin(q3 + q4) +
 0.2*sin(q3 + q4 + q5)) +
 0.000075*dotq1*sin(
 q3 + q4 +
 q5) - (0.0007335*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(1.0*dotq3*cos(
 q3) + (-0.0204499*dotq4 - 0.0204499*dotq5)*cos(
 q4 + q5) + 0.0102249*dotq4*cos(q3 + q4 + q5) +
 0.0102249*dotq5*cos(q3 + q4 + q5) +
 3.33333*dotq3*sin(q3) + 0.102249*dotq4*sin(q4) -
 0.0511247*dotq4*sin(q3 + q4))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + (-0.06*cos(q3 + q4 + q5) +
 0.06*sin(q3 + q4) + ((0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))))*(-0.02445*dotq3*cos(q3) +
 0.000375*dotq4*cos(q3 + q4) +
 0.007335*dotq3*sin(q3) - 0.00015*dotq4*cos(q5)*sin(q4) -
 0.00015*dotq5*cos(q5)*sin(q4) +
 cos(q4)*(-0.00075*dotq4 - 0.00015*dotq4*sin(q5) -
 0.00015*dotq5*sin(q5)) +
 0.000075*dotq4*sin(q3 + q4 + q5) + 
 0.000075*dotq5*sin(
 q3 + q4 + q5)) + (5.625*pow(10, -6)*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5))*((-4346.67*dotq1 - 1304.0*dotq2 +
 5.33333*dotq3)*cos(q3) + (4.0*dotq3 - 8.0*dotq4)*cos(
 q3 - q4) - 4.0*dotq4*cos(q4) + 1.0*dotq4*cos(q3 + q4) -
 0.8*dotq4*cos(q3 - q5) +
 2.66667*dotq3*cos(q3 - q4 - q5) -
 5.33333*dotq4*cos(q3 - q4 - q5) -
 5.33333*dotq5*cos(q3 - q4 - q5) - 3.2*dotq5*cos(q5) +
 0.8*dotq4*cos(q3 + q5) + 0.8*dotq5*cos(q3 + q5) -
 5.33333*dotq4*cos(q4 + q5) -
 5.33333*dotq5*cos(q4 + q5) +
 1.33333*dotq4*cos(q3 + q4 + q5) +
 1.33333*dotq5*cos(q3 + q4 + q5) + 1304.0*dotq1*sin(q3) -
 4346.67*dotq2*sin(q3) - 220.569*dotq3*sin(q3) -
 4.0*dotq4*sin(q3) - 0.4*dotq5*sin(q3) +
 13.3333*dotq3*sin(q3 - q4) -
 26.6667*dotq4*sin(q3 - q4) + 26.6667*dotq4*sin(q4) -
 6.66667*dotq4*sin(q3 + q4) -
 0.8*dotq3*sin(q3 - q4 - q5) +
 1.6*dotq4*sin(q3 - q4 - q5) +
 1.6*dotq5*sin(q3 - q4 - q5) - 0.8*dotq4*sin(q4 + q5) -
 0.8*dotq5*sin(q4 + q5) + 0.2*dotq4*sin(q3 + q4 + q5) +
 0.2*dotq5*sin(q3 + q4 + q5)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))) + (0.0000149625*(8.0*(1.0*cos(q3) +
 1.66667*cos(q3 + q4) - 0.375*cos(q3 + q4 + q5) +
 0.375*sin(q3 + q4) - 0.5*sin(q3 - q5) +
 0.5*sin(q3 + q5) +
 1.66667*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5)) +
 dotq4*((3.33333 + 1.0*cos(q5)*sin(q4) +

 cos(q4)*(1.0 + 1.0*sin(q5)))*(-13.3333*cos(
 q3 + q4) + 3.0*cos(q3 + q4 + q5) -
 3.0*sin(q3 + q4) - 13.3333*sin(q3 + q4 + q5)) -
 3.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(1.0*cos(q3 + q4) -
 1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
 4.44444*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
 4.44444*sin(q3 + q4) + 1.0*sin(q3 + q4 + q5))) +
 dotq5*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 4.0*sin(q3 - q5) - 4.0*sin(q3 + q5) -
 13.3333*sin(q3 + q4 + q5)) -
 3.0*(1.0*cos(q4)*cos(q5) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) -
 1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
 4.44444*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
 4.44444*sin(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((0.00728604 -
 0.0024814*cos(q3) + 0.00015*cos(q3 - q4) +
 0.00015*cos(q4) - 9.*pow(10, -6)*cos(q3 - q4 - q5) -
 4.5*pow(10, -6)*cos(q4 + q5) - 0.00006*sin(q3) -
 0.000045*sin(q3 - q4) + 0.0000225*sin(q4) -
 0.00003*sin(q3 - q4 - q5) + 0.000018*sin(q5) + 
 0.00003*sin(q4 + q5))*(1.0*dotq4*cos(q3 + q4)*cos(
 q3 + q4 + q5) + 2.0*dotq4*sin(q3) + 1.0*dotq5*sin(q3) +
 3.33333*dotq4*cos(q4)*sin(q3) +
 3.33333*dotq4*cos(q5)*sin(q3)*sin(q4) +
 3.33333*dotq5*cos(q5)*sin(q3)*sin(q4) -
 0.5*dotq4*sin(2*(q3 + q4)) + 2.0*dotq4*sin(q3)*sin(q5) +
 1.0*dotq5*sin(q3)*sin(q5) +
 3.33333*dotq4*cos(q4)*sin(q3)*sin(q5) +
 3.33333*dotq5*cos(q4)*sin(q3)*sin(q5) +
 cos(q3)*((-3.33333*dotq4 - 3.33333*dotq5)*cos(q4)*cos(
 q5) + sin(
 q4)*(3.33333*dotq4 + 3.33333*dotq4*sin(q5) +
 3.33333*dotq5*sin(q5))) -
 1.0*dotq4*sin(q3 + q4)*sin(q3 + q4 + q5) -
 1.0*dotq5*sin(q3 + q4)*sin(q3 + q4 + q5) +
 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((-0.0113 +
 0.0489*cos(q3) - 0.01467*sin(q3) +
 0.00015*cos(q5)*sin(q4) +
 cos(q4)*(0.00075 + 0.00015*sin(q5)))*(-0.1*dotq5*(-1.0 +
 1.0*cos(q3) -
 0.3*sin(q3))*(sin(
 q3)*(1.0 +
 3.33333*cos(q5)*sin(
 q4) + (1.0 + 3.33333*cos(q4))*sin(q5)) +
 cos(q3)*(-3.33333*cos(q4)*cos(q5) +
 3.33333*sin(q4)*sin(q5))) -
 0.1*dotq4*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(sin(
 q3)*(2.0 + 3.33333*cos(q5)*sin(q4) + 2.0*sin(q5) +
 cos(q4)*(3.33333 + 3.33333*sin(q5))) +
 cos(q3)*(-3.33333*cos(q4)*cos(q5) +
 sin(q4)*(3.33333 + 3.33333*sin(q5)))) +
 0.03*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(1.0*cos(2*q3 + q4) +
 3.33333*cos(q3 + q4 + q5) -
 3.33333*cos(2*q3 + q4 + q5) - 3.33333*sin(q3 + q4) +
 3.33333*sin(2*q3 + q4) +
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + (0.00369 -
 0.01467*cos(q3) + 0.00015*cos(q4)*cos(q5) -
 0.0489*sin(q3) - 0.00075*sin(q4) -
 0.00015*sin(q4)*sin(
 q5))*(dotq5*(-((
 0.06*cos(
 q3 + q4 + q5)*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(

 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5)))) + ((0.0009*cos(
 q5)*sin(q4) +
 0.0009*cos(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) - ((0.0009*cos(q4)*cos(q5) - 
 0.0009*sin(q4)*sin(q5))*(0.000225 +
 0.00045*cos(q3) - 0.0009*cos(q4)*cos(q5) +
 0.0015*sin(q3) + 0.0009*sin(q4) +
 0.0009*sin(q4)*sin(q5))*(-0.06*cos(q3 + q4) -
 0.06*sin(q3 + q4 + q5)))/pow((0.003 +
 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))), 2) +
 0.06*sin(q3 + q4 + q5)) +
 dotq4*(0.06*cos(
 q3 + q4) + ((-0.06*cos(q3 + q4 + q5) +
 0.06*sin(q3 + q4))*(0.000225 + 0.00045*cos(q3) -
 0.0009*cos(q4)*cos(q5) + 0.0015*sin(q3) +
 0.0009*sin(q4) + 0.0009*sin(q4)*sin(q5)))/(
 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 +
 0.0009*sin(q5))) + ((0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5)))*(-0.06*cos(
 q3 + q4) - 0.06*sin(q3 + q4 + q5)))/(

 0.003 + 0.0009*cos(q5)*sin(q4) +
 cos(q4)*(0.0009 + 0.0009*sin(q5))) +
 0.06*sin(
 q3 + q4 +
 q5) - (0.06*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(-0.25 - 0.5*cos(q3) +
 1.0*cos(q4)*cos(q5) - 1.66667*sin(q3) -
 1.0*sin(q4) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q4) +
 1.0*sin(q3 + q4 + q5)))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)) + (0.1*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(-1.2*cos(q3) - 2.0*cos(q3 + q4) +
 1.0*cos(2*q3 + q4) + 0.15*cos(q3 + q4 + q5) +
 0.3*cos(2*q3 + q4 + q5) - 0.15*sin(q3 + q4) -
 0.3*sin(2*q3 + q4) + 0.6*sin(q3 - q5) -
 0.6*sin(q3 + q5) - 2.0*sin(q3 + q4 + q5) +
 1.0*sin(2*q3 + q4 + q5)))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2))))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)));
	gsl_matrix_set(C2, 1, 0, value);

//C2(1,1)

 value = -4.5*pow(10, -6)*dotq4*cos(q5) -
 0.000075*dotq2*cos(q3 + q4 + q5) -
 4.5*pow(10, -6)*dotq3*(1.0*cos(q3 + q5) + 1.66667*cos(q4 + q5) +
 1.66667*cos(q3 + q4 + q5) + 0.5*sin(q4 + q5) +
 0.25*sin(q3 + q4 + q5)) + 0.000075*dotq1*sin(q3 + q4 + q5) - (
 7.5*pow(10, -6)*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(1.0*dotq3*cos(q4)*cos(
 q5) + (-2.0*dotq3 - 1.0*dotq4 - 1.0*dotq5)*cos(q3 + q4 + q5) -
 1.0*dotq3*sin(q4)*sin(q5))*sin(q3 + q4 + q5))/(
 3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + (0.000075*sin(
 q3 + q4 + q5)*(0.06*dotq4*cos(q3 - q5) +
 0.2*dotq3*cos(q3 - q4 - q5) + 0.12*dotq3*cos(q5) +
 1.0*dotq2*cos(q4)*cos(q5) - 0.12*dotq3*cos(q3 + q5) +
 0.06*dotq4*cos(q3 + q5) - 0.06*dotq5*cos(q3 + q5) +
 0.2*dotq3*cos(q4 + q5) - 0.1*dotq4*cos(q4 + q5) -
 0.1*dotq5*cos(q4 + q5) - 0.2*dotq3*cos(q3 + q4 + q5) - 
 0.1*dotq4*cos(q3 + q4 + q5) - 0.1*dotq5*cos(q3 + q4 + q5) +
 0.06*dotq3*sin(q3) - 1.0*dotq1*cos(q5)*sin(q4) -
 0.06*dotq3*sin(q3 - q4 - q5) - 1.0*dotq1*cos(q4)*sin(q5) -
 1.0*dotq2*sin(q4)*sin(q5) + 0.03*dotq3*sin(q4 + q5) -
 0.03*dotq4*sin(q4 + q5) - 0.03*dotq5*sin(q4 + q5) -
 0.03*dotq3*sin(q3 + q4 + q5) -
 0.015*dotq4*sin(q3 + q4 + q5) -
 0.015*dotq5*sin(q3 + q4 + q5)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + (0.000075*dotq3*cos(q5)*sin(q4) +
 0.000075*dotq3*cos(q4)*sin(
 q5) + (-0.00015*dotq3 - 0.000075*dotq4 - 0.000075*dotq5)*sin(
 q3 + q4 + q5))*(-0.06*cos(q3 + q4 + q5) - (
 0.03*(0.5 + 1.0*cos(q3) - 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) +
 2.0*sin(q4) + 2.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5))/(
 3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))) + ((0.0000186406 +
 4.5*pow(10, -6)*cos(q3) - 4.5*pow(10, -6)*cos(q4 + q5) -
 2.25*pow(10, -6)*cos(q3 + q4 + q5) + 9.*pow(10, -6)*sin(q3 + q5) +
 0.000015*sin(q4 + q5) +
 0.000015*sin(q3 + q4 + q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
 q3 + q5) + (-3.33333*dotq4 - 3.33333*dotq5 +
 1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
 0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((4.0*cos(q3 + q5) +
 13.3333*cos(q3 + q4 + q5) - 4.0*sin(q3) +
 3.0*sin(q3 + q4 +
 q5))*(0.0000149625*(dotq4*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 13.3333*sin(q3 + q4 + q5)) -
 4.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3) +
 0.75*sin(q3 + q4 + q5))) +
 dotq5*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 4.0*sin(q3 + q5) - 13.3333*sin(q3 + q4 + q5)) -
 4.0*(1.0*cos(q4)*cos(q5) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3) +
 0.75*sin(q3 + q4 + q5))) +
 4.0*(1.0*cos(q3) - 0.75*cos(q3 + q4 + q5) +
 1.0*sin(q3 + q5) +
 3.33333*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))) +
 0.0000149625*((-1.0*dotq4 - 0.5*dotq5)*cos(
 q3 + q5) + (-3.33333*dotq4 - 3.33333*dotq5 +
 1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
 0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5)))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 3) + (-0.06*cos(q3 + q4 + q5) - (
 0.03*(0.5 + 1.0*cos(q3) - 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) +
 2.0*sin(q4) + 2.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5))/(
 3.33333 + 1.0*cos(q5)*sin(q4) + 
 cos(q4)*(1.0 + 1.0*sin(q5))))*(-0.00015*dotq3*sin(
 q3 + q4 + q5) - 0.00015*dotq4*sin(q3 + q4 + q5) -
 0.00015*dotq5*sin(
 q3 + q4 +
 q5) + (0.0489*(1.0*dotq3*cos(q3) - 0.3*dotq3*sin(q3) +
 0.00306748*dotq4*cos(q5)*sin(q4) +
 0.00306748*dotq5*cos(q5)*sin(q4) +
 cos(q4)*(0.0153374*dotq4 + 0.00306748*dotq4*sin(q5) +
 0.00306748*dotq5*sin(q5)))*sin(
 q3 + q4 + q5))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + ((0.00369 - 0.01467*cos(q3) +
 0.00015*cos(q4)*cos(q5) - 0.0489*sin(q3) -
 0.00075*sin(q4) -
 0.00015*sin(q4)*sin(q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
 q3 + q5) + (-3.33333*dotq4 - 3.33333*dotq5 +
 1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
 0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))), 2) - (0.01467*(1.0*dotq5*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) -
 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) + 2.0*sin(q4) +
 2.0*sin(q4)*sin(q5))*(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) +
 2.0*(1.0*cos(q4)*cos(q5) - 1.0*sin(q4)*sin(q5))*(-0.25 -
 0.5*cos(q3) + 1.0*cos(q4)*cos(q5) - 1.66667*sin(q3) -
 1.0*sin(q4) - 1.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5) +
 1.0*(2.0*cos(q5)*sin(q4) + 2.0*cos(q4)*sin(q5))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)))*sin(
 q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +

 1.0*dotq4*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) -
 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) + 2.0*sin(q4) +
 2.0*sin(q4)*sin(q5))*(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) -
 2.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(0.25 + 0.5*cos(q3) -
 1.0*cos(q4)*cos(q5) + 1.66667*sin(q3) + 1.0*sin(q4) +
 1.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5) +
 2.0*(1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)))*sin(
 q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +
 3.33333*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(0.6*cos(q3) -
 0.15*cos(q3 + q4 + q5) - 0.3*cos(2*q3 + q4 + q5) +
 0.6*sin(q3 + q5) + 2.0*sin(q3 + q4 + q5) -
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)) + (0.1*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*sin(
 q3 + q4 + q5)*(0.00015*dotq3*cos(q3 + q4 + q5) +
 0.00015*dotq4*cos(q3 + q4 + q5) + 
 0.00015*dotq5*cos(
 q3 + q4 +
 q5) + (0.01467*(1.0*dotq3*cos(
 q3) + (-0.0102249*dotq4 - 0.0102249*dotq5)*cos(q4)*cos(
 q5) + 3.33333*dotq3*sin(q3) +
 0.0511247*dotq4*sin(q4) +
 0.0102249*dotq4*sin(q4)*sin(q5) +
 0.0102249*dotq5*sin(q4)*sin(q5))*sin(
 q3 + q4 + q5))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + ((-0.0113 + 0.0489*cos(q3) -
 0.01467*sin(q3) + 0.00015*cos(q5)*sin(q4) +
 cos(q4)*(0.00075 + 0.00015*sin(q5)))*((-1.0*dotq4 -
 0.5*dotq5)*cos(
 q3 + q5) + (-3.33333*dotq4 - 3.33333*dotq5 +
 1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
 0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))), 2) + (0.489*(-0.015*dotq5*(1.0*cos(
 q3 + q5) + 6.66667*cos(q3 + q4 + q5) +
 1.0*cos(q3 + 2*q4 + q5) - 2.0*sin(q3))*(3.33333 -
 3.33333*cos(q3) + 1.0*sin(q3)) -

 0.03*dotq4*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3))*(3.33333 -
 3.33333*cos(q3) + 1.0*sin(q3)) +
 0.1*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(1.0*cos(q3 + q4 + q5) -
 1.0*cos(2*q3 + q4 + q5) +
 0.3*sin(2*q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) - (1.0*sin(
 q3 + q4 + q5)*(-0.000075*dotq2*cos(q3 + q4 + q5) +
 dotq4*(0.000015*cos(q4 + q5) + 0.000015*cos(q3 + q4 + q5) +
 4.5*pow(10, -6)*sin(q4 + q5) + 2.25*pow(10, -6)*sin(q3 + q4 + q5)) +
 dotq5*(9.*pow(10, -6)*cos(q3 + q5) + 0.000015*cos(q4 + q5) +
 0.000015*cos(q3 + q4 + q5) + 4.5*pow(10, -6)*sin(q4 + q5) +
 2.25*pow(10, -6)*sin(q3 + q4 + q5)) +
 4.5*pow(10, -6)*dotq3*(1.0*cos(q3 + q5) +
 1.66667*cos(q3 + q4 + q5) - 0.5*sin(q3) +
 0.25*sin(q3 + q4 + q5)) +
 0.000075*dotq1*sin(
 q3 + q4 +
 q5) - (0.0007335*(-1.0 + 1.0*cos(q3) -
 0.3*sin(q3))*(1.0*dotq3*cos(
 q3) + (-0.0204499*dotq4 - 0.0204499*dotq5)*cos(
 q4 + q5) + 0.0102249*dotq4*cos(q3 + q4 + q5) +
 0.0102249*dotq5*cos(q3 + q4 + q5) +
 3.33333*dotq3*sin(q3) + 0.102249*dotq4*sin(q4) -
 0.0511247*dotq4*sin(q3 + q4))*sin(
 q3 + q4 + q5))/(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) + (1.125*pow(10, -6)*sin(
 q3 + q4 +
 q5)*((-21733.3*dotq1 - 6520.0*dotq2 + 26.6667*dotq3)*cos(
 q3) + (20.0*dotq3 - 40.0*dotq4)*cos(q3 - q4) -
 20.0*dotq4*cos(q4) + 5.0*dotq4*cos(q3 + q4) - 
 4.0*dotq4*cos(q3 - q5) +
 13.3333*dotq3*cos(q3 - q4 - q5) -
 26.6667*dotq4*cos(q3 - q4 - q5) -
 26.6667*dotq5*cos(q3 - q4 - q5) - 16.0*dotq5*cos(q5) +
 4.0*dotq4*cos(q3 + q5) + 4.0*dotq5*cos(q3 + q5) -
 26.6667*dotq4*cos(q4 + q5) -
 26.6667*dotq5*cos(q4 + q5) +
 6.66667*dotq4*cos(q3 + q4 + q5) +
 6.66667*dotq5*cos(q3 + q4 + q5) + 6520.0*dotq1*sin(q3) -
 21733.3*dotq2*sin(q3) - 1102.84*dotq3*sin(q3) -
 20.0*dotq4*sin(q3) - 2.0*dotq5*sin(q3) +
 66.6667*dotq3*sin(q3 - q4) - 133.333*dotq4*sin(q3 - q4) +
 133.333*dotq4*sin(q4) - 33.3333*dotq4*sin(q3 + q4) -
 4.0*dotq3*sin(q3 - q4 - q5) +
 8.0*dotq4*sin(q3 - q4 - q5) +
 8.0*dotq5*sin(q3 - q4 - q5) - 4.0*dotq4*sin(q4 + q5) -
 4.0*dotq5*sin(q4 + q5) + 1.0*dotq4*sin(q3 + q4 + q5) +
 1.0*dotq5*sin(q3 + q4 + q5)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))) + (0.0007335*(1.0*dotq3*cos(q3) -
 0.0153374*dotq4*cos(q3 + q4) - 0.3*dotq3*sin(q3) +
 0.00613497*dotq4*cos(q5)*sin(q4) +
 0.00613497*dotq5*cos(q5)*sin(q4) +
 cos(q4)*(0.0306748*dotq4 + 0.00613497*dotq4*sin(q5) +
 0.00613497*dotq5*sin(q5)) -
 0.00306748*dotq4*sin(q3 + q4 + q5) -
 0.00306748*dotq5*sin(q3 + q4 + q5))*(cos(
 q3 + q4 + q5)*(6.66667 + 2.0*cos(q5)*sin(q4) +
 cos(q4)*(2.0 + 2.0*sin(q5))) + (0.5 + 1.0*cos(q3) -
 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) + 2.0*sin(q4) +
 2.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 +
 1.0*sin(q5))) + (0.0000149625*(dotq4*((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 13.3333*sin(q3 + q4 + q5)) -
 4.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3) +
 0.75*sin(q3 + q4 + q5))) +
 dotq5*((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.0*cos(q3 + q4 + q5) -
 4.0*sin(q3 + q5) - 13.3333*sin(q3 + q4 + q5)) -
 4.0*(1.0*cos(q4)*cos(q5) -
 1.0*sin(q4)*sin(q5))*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3) +
 0.75*sin(q3 + q4 + q5))) +
 4.0*(1.0*cos(q3) - 0.75*cos(q3 + q4 + q5) +
 1.0*sin(q3 + q5) +
 3.33333*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((0.00728604 -
 0.0024814*cos(q3) + 0.00015*cos(q3 - q4) +
 0.00015*cos(q4) - 9.*pow(10, -6)*cos(q3 - q4 - q5) -
 4.5*pow(10, -6)*cos(q4 + q5) - 0.00006*sin(q3) -
 0.000045*sin(q3 - q4) + 0.0000225*sin(q4) -
 0.00003*sin(q3 - q4 - q5) + 0.000018*sin(q5) + 
 0.00003*sin(q4 + q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
 q3 + q5) + (-3.33333*dotq4 - 3.33333*dotq5 +
 1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
 0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
 1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
 0.5*dotq5*sin(2*(q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) - (0.03*(0.00369 -
 0.01467*cos(q3) + 0.00015*cos(q4)*cos(q5) -
 0.0489*sin(q3) - 0.00075*sin(q4) -
 0.00015*sin(q4)*sin(
 q5))*(1.0*dotq5*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) -
 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) +
 2.0*sin(q4) + 2.0*sin(q4)*sin(q5))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5))) +
 2.0*(1.0*cos(q4)*cos(q5) - 1.0*sin(q4)*sin(q5))*(-0.25 -
 0.5*cos(q3) + 1.0*cos(q4)*cos(q5) -
 1.66667*sin(q3) - 1.0*sin(q4) -
 1.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5) +
 1.0*(2.0*cos(q5)*sin(q4) +
 2.0*cos(q4)*sin(q5))*(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*sin(q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +
 1.0*dotq4*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) -
 2.0*cos(q4)*cos(q5) + 3.33333*sin(q3) +
 2.0*sin(q4) + 2.0*sin(q4)*sin(q5))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5))) -
 2.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(0.25 + 0.5*cos(q3) -
 1.0*cos(q4)*cos(q5) + 1.66667*sin(q3) +
 1.0*sin(q4) + 1.0*sin(q4)*sin(q5))*sin(
 q3 + q4 + q5) +
 2.0*(1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*sin(q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +
 3.33333*(1.0*dotq4*cos(
 q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(0.6*cos(q3) -
 0.15*cos(q3 + q4 + q5) - 0.3*cos(2*q3 + q4 + q5) +
 0.6*sin(q3 + q5) + 2.0*sin(q3 + q4 + q5) -
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + ((-0.0113 +
 0.0489*cos(q3) - 0.01467*sin(q3) +
 0.00015*cos(q5)*sin(q4) +
 cos(q4)*(0.00075 +
 0.00015*sin(q5)))*(-0.015*dotq5*(1.0*cos(q3 + q5) +
 6.66667*cos(q3 + q4 + q5) + 1.0*cos(q3 + 2*q4 + q5) -
 2.0*sin(q3))*(3.33333 - 3.33333*cos(q3) +
 1.0*sin(q3)) -
 0.03*dotq4*(1.0*cos(q3 + q5) +
 3.33333*cos(q3 + q4 + q5) - 1.0*sin(q3))*(3.33333 -
 3.33333*cos(q3) + 1.0*sin(q3)) +
 0.1*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(1.0*cos(q3 + q4 + q5) -
 1.0*cos(2*q3 + q4 + q5) +
 0.3*sin(2*q3 + q4 + q5))))/pow((3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)))/(3.33333 +
 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) - (4.5*pow(10, -6)*cos(
 q3 + q4 +
 q5)*(1.0*dotq5*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4)*cos(q5) +
 3.33333*sin(q3) + 2.0*sin(q4) +
 2.0*sin(q4)*sin(q5))*(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) +

 2.0*(1.0*cos(q4)*cos(q5) - 1.0*sin(q4)*sin(q5))*(-0.25 -
 0.5*cos(q3) + 1.0*cos(q4)*cos(q5) - 1.66667*sin(q3) -
 1.0*sin(q4) - 1.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5) +
 1.0*(2.0*cos(q5)*sin(q4) + 2.0*cos(q4)*sin(q5))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)))*sin(
 q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +
 1.0*dotq4*(1.0*cos(
 q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4)*cos(q5) +
 3.33333*sin(q3) + 2.0*sin(q4) +
 2.0*sin(q4)*sin(q5))*(3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))) -
 2.0*(1.0*cos(q4)*cos(q5) +
 sin(q4)*(-1.0 - 1.0*sin(q5)))*(0.25 + 0.5*cos(q3) -
 1.0*cos(q4)*cos(q5) + 1.66667*sin(q3) + 1.0*sin(q4) +
 1.0*sin(q4)*sin(q5))*sin(q3 + q4 + q5) +
 2.0*(1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5)))*(3.33333 +
 1.0*cos(q5)*sin(q4) + cos(q4)*(1.0 + 1.0*sin(q5)))*sin(
 q3 + q4 + q5) -
 2.0*pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2)*sin(q3 + q4 + q5)) +
 3.33333*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(0.6*cos(q3) - 0.15*cos(q3 + q4 + q5) -
 0.3*cos(2*q3 + q4 + q5) + 0.6*sin(q3 + q5) +
 2.0*sin(q3 + q4 + q5) -
 1.0*sin(2*q3 + q4 + q5))))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2) + (0.00015*sin(
 q3 + q4 +
 q5)*(-0.015*dotq5*(1.0*cos(q3 + q5) +
 6.66667*cos(q3 + q4 + q5) + 1.0*cos(q3 + 2*q4 + q5) -
 2.0*sin(q3))*(3.33333 - 3.33333*cos(q3) + 1.0*sin(q3)) -
 0.03*dotq4*(1.0*cos(q3 + q5) + 3.33333*cos(q3 + q4 + q5) -
 1.0*sin(q3))*(3.33333 - 3.33333*cos(q3) + 1.0*sin(q3)) +
 0.1*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
 q3 + q4 + q5))*(1.0*cos(q3 + q4 + q5) -
 1.0*cos(2*q3 + q4 + q5) +
 0.3*sin(2*q3 + q4 + q5))))/pow((3.33333 + 1.0*cos(q5)*sin(q4) +
 cos(q4)*(1.0 + 1.0*sin(q5))), 2);
 gsl_matrix_set(C2, 1, 1, value);
}

void print_matrix(const gsl_matrix *m){

    int i, j;

    for (i = 0; i < m->size1; i++){
        for (j = 0; j < m->size2; j++)
            printf("%g ", gsl_matrix_get(m, i, j));
        printf("\n");
    }
}
/*
int main(){
	int i;
	state rob;
	rob.q1 = rob.q2 = rob.q3 = rob.q4 = rob.q5 = rob.q6 = 1;
	rob.dotq1 = rob.dotq2 = rob.dotq3 = rob.dotq4 = rob.dotq5 = rob.dotq6 = 1;

	gsl_matrix *Tsee, *M1, *C1, *S2, *M2, *C2;
	gsl_vector *G1, *G2;

	Tsee = gsl_matrix_alloc(4,4);
	M1 = gsl_matrix_alloc(2,2);
	M2 = gsl_matrix_alloc(2,2);
	C1 = gsl_matrix_alloc(2,2);
	S2 = gsl_matrix_alloc(4,2);
	C2 = gsl_matrix_alloc(2,2);
	G1 = gsl_vector_alloc(2);
	G2 = gsl_vector_alloc(2);

	update_S2(S2, rob);
	print_matrix(S2);

	for(i = 0; i < 1000000; i++){
		update_kyn(Tsee, rob);
		update_M1(M1, rob);
		update_C1(C1, rob);
		update_G1(G1, rob);
		update_S2(S2, rob);
		update_G2(G2, rob);
		update_M2(M2, rob);
		update_C2(C2, rob);
		rob.q1++;
		rob.q2++;
		rob.q3++;
		rob.q4++;
		rob.q5++;
		rob.q6++;
	}

	return 0;
}*/