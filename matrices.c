#include "matrices.h"
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/*state robot;
dot_state dot_robot;*/

void update_kyn(float Tsee[4][4], state robot){
	float value;
	float q1, q2, q3, q4, q5; 
	q1 = robot.q1;
	q2 = robot.q2;
	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;

	// Prima riga
	value = cos(q4 + q5)*sin(q3) + cos(q3)*sin(q4 + q5);
	Tsee[0][0] = value;

  	value = cos(q3)*cos(q4 + q5) - sin(q3)*sin(q4 + q5);
  	Tsee[0][1] =  value;

  	value = 0.0;
  	Tsee[0][2] = value;

  	value = -(15.0/2.0) + q1 + (15*cos(q3))/2.0 - 3.0*sin(q3) -
			sin(q3)*(3.0/2.0 - 6.0*cos(q4 + q5) + 6.0*sin(q4)) +
			cos(q3)*(15.0/2.0 + 6.0*cos(q4) + 6.0*sin(q4 + q5));
	Tsee[0][3] = value;

	// Seconda riga
	value = -cos(q3)*cos(q4 + q5) + sin(q3)*sin(q4 + q5);
	Tsee[1][0] = value;

	value = cos(q4 + q5)*sin(q3) + cos(q3)*sin(q4 + q5);
	Tsee[1][1] =  value;

	value = 0.0;
	Tsee[1][2] =  value;

	value = 1.5 + q2 + 3.0*cos(q3) + (15*sin(q3))/2 +
			cos(q3)*(3.0/2.0 - 6*cos(q4 + q5) + 6*sin(q4)) +
			sin(q3)*(15.0/2.0 + 6*cos(q4) + 6*sin(q4 + q5));
	Tsee[1][3] =  value;

	// Terza riga
	Tsee[2][0] =  0.0;
	Tsee[2][1] =  0.0;
	Tsee[2][2] =  1.0;
	Tsee[2][3] =  0.0;

	// Quarta riga
	Tsee[3][0] = 0.0;
	Tsee[3][1] = 0.0;
	Tsee[3][2] = 0.0;
	Tsee[3][3] = 1.0;
}

void update_S2(float S2[4][2], state robot){
	float value;
	float q3, q4, q5;
	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;

	//PRIMA RIGA
	value = -6*cos(q3)*(cos(q4 + q5) - sin(q4)) +
			6*sin(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5));
	S2[0][0] = value;

	value =   -6*cos(q3 + q4 + q5) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*sin(q3 + q4 + q5))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5));
	S2[0][1] = value;

	//SECONDA RIGA
	value = (7.5*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5));
	S2[1][0] = value;

	value = (7.5*(-1.0 + 1.0*cos(q3) - 0.4*sin(q3))*sin(q3 + q4 + q5))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5));
	S2[1][1] = value;

	//TERZA RIGA
	value = (-1.0*cos(q3 + q4) -1.0*sin(q3 + q4 + q5))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5));
	S2[2][0] = value;

	value = -((1.0*sin(q3 + q4 + q5))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)));
	S2[2][1] = value;

	//QUARTA RIGA
	value =(3.0*cos(q3 + q4) - 4.0*cos(q3 - q5) + 4.0*cos(q3 + q5) +
			10.0*cos(q3 + q4 + q5) - 8.0*sin(q3) - 10.0*sin(q3 + q4) +
			3.0*sin(q3 + q4 + q5))/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5));
	S2[3][0] = value;

	value =(4.0*cos(q3 + q5) + 10.0*cos(q3 + q4 + q5) - 4.0*sin(q3) +
			3.0*sin(q3 + q4 + q5))/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5));
	S2[3][1] = value;
}

void update_M1(float M1[2][2], state robot){
	float value;
	float q5;
	q5 = robot.q5;

	//PRIMA RIGA
	value =180.0*(4.57118 + sin(q5));
	M1[0][0] = value;

	value =231.406 + 90*sin(q5);
	M1[0][1] = value;

	//SECONDA RIGA
	value = 231.406 + 90*sin(q5);
	M1[1][0] = value;

	value =  231.406;
	M1[1][1] = value;
}

void update_G1(float G1[2], state robot){
	float value;
	float q3, q4, q5;
	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;
	
	//PRIMA RIGA
	value = 0.0 + 735.0*cos(q3 + q4) + 147.0*sin(q3 + q4 + q5);
	G1[0] = value;

	//SECONDA RIGA
	value = 0.0 + 147.0*sin(q3 + q4 + q5);
	G1[1] = value;
}

void update_C1(float C1[2][2], state robot, dot_state dot_robot){
	float value;
	float q3, q4, q5, dotq1, dotq2, dotq3, dotq4, dotq5;
	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;
	dotq1 = dot_robot.dq1;
	dotq2 = dot_robot.dq2;
	dotq3 = dot_robot.dq3;
	dotq4 = dot_robot.dq4;
	dotq5 = dot_robot.dq5;

	//PRIMA RIGA
	value = -112.5*dotq3*cos(q4) + (37.5*dotq1 - 56.25*dotq3)*cos(q3 + q4) +
			180.0*dotq5*cos(q5) - 56.25*dotq3*cos(q4 + q5) -
			7.5*dotq2*cos(q3 + q4 + q5) - 56.25*dotq3*cos(q3 + q4 + q5) +
			281.25*dotq3*sin(q4) + 37.5*dotq2*sin(q3 + q4) +
			281.25*dotq3*sin(q3 + q4) - 22.5*dotq3*sin(q4 + q5) +
			7.5*dotq1*sin(q3 + q4 + q5) - 11.25*dotq3*sin(q3 + q4 + q5);
	C1[0][0] = value;

	value = 90.0*dotq5*cos(q5) - 56.25*dotq3*cos(q4 + q5) -
			7.5*dotq2*cos(q3 + q4 + q5) - 56.25*dotq3*cos(q3 + q4 + q5) -
			22.5*dotq3*sin(q4 + q5) + 7.5*dotq1*sin(q3 + q4 + q5) -
			11.25*dotq3*sin(q3 + q4 + q5);
	C1[0][1] = value;

	//SECONDA RIGA
	value = -90.0*dotq4*cos(q5) +
			45*dotq5*cos(q5) - 7.5*dotq2*cos(q3 + q4 + q5) -
			45.0*dotq3*(1.0*cos(q3 - q5) + 1.0*cos(q3 + q5) + 1.25*cos(q4 + q5) +
			1.25*cos(q3 + q4 + q5) + 0.5*sin(q4 + q5) +
			0.25*sin(q3 + q4 + q5)) +
			7.5*dotq1*sin(q3 + q4 + q5);
	C1[1][0] = value;

	value =-45.0*dotq4*cos(q5) -
			45.0*dotq3*cos(q3 + q5) - 56.25*dotq3*cos(q4 + q5) -
			7.5*dotq2*cos(q3 + q4 + q5) - 56.25*dotq3*cos(q3 + q4 + q5) -
			22.5*dotq3*sin(q4 + q5) + 7.5*dotq1*sin(q3 + q4 + q5) -
			11.25*dotq3*sin(q3 + q4 + q5);
	C1[1][1] = value;
}

void update_G2(float G2[2], state robot){
	float value;
	float q3, q4, q5;
	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;

	//PRIMA RIGA
	value = 0.0 + 735.0*cos(q3 + q4) + (
			29767.5*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) + ((-8305.5 + 29767.5*cos(q3) + 735.0*cos(q4) -
			11907.0*sin(q3) + 147.0*sin(q4 + q5))*(-1.0*cos(q3 + q4) -
			1.0*sin(q3 + q4 + q5)))/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			147.0*sin(q3 + q4 + q5);
	G2[0] =  value;

	//SECONDA RIGA
	value = 0.0 + 147.0*sin(q3 + q4 + q5) + (
			29767.5*(-1.0 + 1.0*cos(q3) - 0.4*sin(q3))*sin(q3 + q4 + q5))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) - (
			1.0*(-8305.5 + 29767.5*cos(q3) + 735.0*cos(q4) - 11907.0*sin(q3) +
			147.0*sin(q4 + q5))*sin(q3 + q4 + q5))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5));
	G2[1] = value;
}

void update_M2(float M2[2][2], state robot){
 	float value;
 	float q3, q4, q5;
 	q3 = robot.q3;
 	q4 = robot.q4;
 	q5 = robot.q5;

	//M2(0,0)
	value = 1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(178771.0 -
			745.625*cos(q3) - 49842.0*cos(2*q3) + 126469.0*cos(q4) +
			186.406*cos(2*q4) - 1864.06*cos(q3 + q4) -
			41005.5*cos(2*(q3 + q4)) - 124605.0*cos(2*q3 + q4) -
			372.813*cos(q3 + 2*q4) + 8262.0*cos(2*(q3 - q5)) +
			7056.0*cos(q4 - q5) - 7056.0*cos(2*q3 + q4 - q5) -
			16524.0*cos(2*q5) + 8262.0*cos(2*(q3 + q5)) -
			21168.0*cos(q4 + q5) - 186.406*cos(2*(q4 + q5)) +
			41005.5*cos(2*(q3 + q4 + q5)) + 21168.0*cos(2*q3 + q4 + q5) +
			35280.0*cos(2*q3 + 2*q4 + q5) - 41310.0*cos(q4 + 2*q5) +
			41310.0*cos(2*q3 + q4 + 2*q5) + 372.813*cos(q3 + 2*(q4 + q5)) +
			21168.0*sin(q4) - 17640.0*sin(2*(q3 + q4)) -
			21168.0*sin(2*q3 + q4) + 372.813*sin(q3 - q5) +
			33183.0*sin(2*q3 - q5) - 41310.0*sin(q4 - q5) +
			41310.0*sin(2*q3 + q4 - q5) + 191278.0*sin(q5) -
			372.813*sin(q3 + q5) - 33183.0*sin(2*q3 + q5) +
			126469.0*sin(q4 + q5) - 1864.06*sin(q3 + q4 + q5) +
			17640.0*sin(2*(q3 + q4 + q5)) - 124605.0*sin(2*q3 + q4 + q5) +
			372.813*sin(2*q4 + q5) - 745.625*sin(q3 + 2*q4 + q5) -
			82010.9*sin(2*q3 + 2*q4 + q5) - 7056.0*sin(q4 + 2*q5) +
			7056.0*sin(2*q3 + q4 + 2*q5));
	M2[0][0] = value;

	//M2(0,1)
	value = 1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(90600.7 -
			372.813*cos(q3) - 25011.0*cos(2*q3) + 64078.3*cos(q4) +
			171.953*cos(2*q4) - 466.016*cos(q3 + q4) -
			78.75*cos(2*(q3 + q4)) - 41703.8*cos(2*q3 + q4) -
			93.2031*cos(q3 + 2*q4) + 3528.0*cos(q4 - q5) - 8262.0*cos(2*q5) +
			8262.0*cos(2*(q3 + q5)) - 10584.0*cos(q4 + q5) -
			104.453*cos(2*(q4 + q5)) + 41016.7*cos(2*(q3 + q4 + q5)) +
			14112.0*cos(2*q3 + q4 + q5) + 17640.0*cos(2*q3 + 2*q4 + q5) -
			20655.0*cos(q4 + 2*q5) + 41310.0*cos(2*q3 + q4 + 2*q5) +
			279.609*cos(q3 + 2*(q4 + q5)) + 10584.0*sin(q4) -
			7056.0*sin(2*q3 + q4) + 93.2031*sin(q3 - q5) +
			8340.75*sin(2*q3 - q5) - 20655.0*sin(q4 - q5) + 95728.9*sin(q5) -
			279.609*sin(q3 + q5) - 24932.3*sin(2*q3 + q5) +
			63740.8*sin(q4 + q5) - 1398.05*sin(q3 + q4 + q5) +
			17640.0*sin(2*(q3 + q4 + q5)) - 83351.3*sin(2*q3 + q4 + q5) +276.406*sin(2*q4 + q5) - 372.813*sin(q3 + 2*q4 + q5) -
			41095.5*sin(2*q3 + 2*q4 + q5) - 3528.0*sin(q4 + 2*q5) +
			7056.0*sin(2*q3 + q4 + 2*q5));
	M2[0][1] = value;

	//M2(1,0)
	value = 1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*(90600.7 - 372.812*cos(q3) -
			25011.0*cos(2*q3) + 64078.3*cos(q4) + 171.953*cos(2*q4) -
			466.016*cos(q3 + q4) - 78.75*cos(2*(q3 + q4)) -
			41703.8*cos(2*q3 + q4) - 93.2031*cos(q3 + 2*q4) +
			3528.0*cos(q4 - q5) - 8262.0*cos(2*q5) + 8262.0*cos(2*(q3 + q5)) -
			10584.0*cos(q4 + q5) - 104.453*cos(2*(q4 + q5)) +
			41016.7*cos(2*(q3 + q4 + q5)) + 14112.0*cos(2*q3 + q4 + q5) +
			17640.0*cos(2*q3 + 2*q4 + q5) - 20655.0*cos(q4 + 2*q5) +
			41310.0*cos(2*q3 + q4 + 2*q5) + 279.609*cos(q3 + 2*(q4 + q5)) +
			10584.0*sin(q4) - 7056.0*sin(2*q3 + q4) + 93.2031*sin(q3 - q5) +
			8340.75*sin(2*q3 - q5) - 20655.0*sin(q4 - q5) + 95728.9*sin(q5) -
			279.609*sin(q3 + q5) - 24932.3*sin(2*q3 + q5) +
			63740.8*sin(q4 + q5) - 1398.05*sin(q3 + q4 + q5) +
			17640.0*sin(2*(q3 + q4 + q5)) - 83351.3*sin(2*q3 + q4 + q5) +
			276.406*sin(2*q4 + q5) - 372.812*sin(q3 + 2*q4 + q5) -
			41095.5*sin(2*q3 + 2*q4 + q5) - 3528.0*sin(q4 + 2*q5) +
			7056.0*sin(2*q3 + q4 + 2*q5));
	M2[1][0] = value;

	//M2(1,1)
	value = 1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(82406.2 -
			186.406*cos(q3) - 8419.5*cos(2*q3) + 43817.0*cos(q4) +
			250.703*cos(2*q4) + 8262.0*cos(2*(q3 + q5)) -
			7056.0*cos(q4 + q5) - 93.2031*cos(2*(q4 + q5)) +
			41028.0*cos(2*(q3 + q4 + q5)) + 7056.0*cos(2*q3 + q4 + q5) +
			41310.0*cos(2*q3 + q4 + 2*q5) + 186.406*cos(q3 + 2*(q4 + q5)) +
			7056.0*sin(q4) + 16867.9*sin(q5) - 186.406*sin(q3 + q5) -
			16681.5*sin(2*q3 + q5) + 43029.5*sin(q4 + q5) -
			932.031*sin(q3 + q4 + q5) + 17640.0*sin(2*(q3 + q4 + q5)) -
			42097.5*sin(2*q3 + q4 + q5) + 343.906*sin(2*q4 + q5) -
			186.406*sin(q3 + 2*q4 + q5) - 157.5*sin(2*q3 + 2*q4 + q5) +
			7056.0*sin(2*q3 + q4 + 2*q5));
	M2[1][1] = value;
}

void update_C2(float C2[2][2], state robot, dot_state dot_robot){

	float value;
	float q3, q4, q5;
    float dotq1, dotq2, dotq3, dotq4, dotq5;

	q3 = robot.q3;
	q4 = robot.q4;
	q5 = robot.q5;
    dotq1 = dot_robot.dq1;
    dotq2 = dot_robot.dq2;
    dotq3 = dot_robot.dq3;
    dotq4 = dot_robot.dq4;
    dotq5 = dot_robot.dq5;

	//C2(0,0)
	value = 37.5*dotq1*cos(q3 + q4) + 180.0*dotq5*cos(q5) -
			7.5*dotq2*cos(q3 + q4 + q5) + 37.5*dotq2*sin(q3 + q4) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))*56.25*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(1.0*dotq3*cos(
			q4 + q5) + (-2.0*dotq3 - 1.0*dotq4 - 1.0*dotq5)*cos(
			q3 + q4 + q5) - 5.0*dotq3*sin(q4) + 10.0*dotq3*sin(q3 + q4) +
			5.0*dotq4*sin(q3 + q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))) +
			1/2*dotq3*(-225.0*cos(q4) - 112.5*cos(q3 + q4) -
			112.5*cos(q4 + q5) - 112.5*cos(q3 + q4 + q5) + 562.5*sin(q4) +
			562.5*sin(q3 + q4) - 45.0*sin(q4 + q5) -
			22.5*sin(q3 + q4 + q5)) +
			7.5*dotq1*sin(
			q3 + q4 + q5) + (-6*cos(q3)*(cos(q4 + q5) - sin(q4)) +
			6*sin(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 
			1.0*sin(q4 + q5)))*(37.5*dotq3*cos(
			q4) + (-75.0*dotq3 - 37.5*dotq4)*cos(q3 + q4) +
			7.5*dotq3*sin(q4 + q5) - 15.0*dotq3*sin(q3 + q4 + q5) -
			7.5*dotq4*sin(q3 + q4 + q5) - 7.5*dotq5*sin(q3 + q4 + q5)) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*37.5*(1.0*cos(q3 + q4) +
			1.0*sin(q3 + q4 + q5))*(-6.0*dotq3*cos(
			q3 - q4) + (1.0*dotq1 - 3.0*dotq3 + 3.0*dotq4)*cos(q4) +
			3.0*dotq3*cos(q3 + q4) + 1.5*dotq4*cos(q3 + q4) -
			2.4*dotq3*cos(q3 - q5) + 2.4*dotq5*cos(q3 - q5) -
			3.0*dotq3*cos(q3 - q4 - q5) + 2.4*dotq3*cos(q3 + q5) +
			2.4*dotq5*cos(q3 + q5) - 0.2*dotq2*cos(q4 + q5) -
			3.0*dotq3*cos(q4 + q5) + 1.5*dotq4*cos(q4 + q5) +
			3.0*dotq5*cos(q4 + q5) + 3.0*dotq3*cos(q3 + q4 + q5) +
			1.5*dotq4*cos(q3 + q4 + q5) + 1.5*dotq5*cos(q3 + q4 + q5) -
			1.5*dotq5*cos(q3)*cos(q3 + q4 + q5) - 12.0*dotq3*sin(q3) +
			0.6*dotq5*cos(q3 + q4 + q5)*sin(q3) - 15.0*dotq3*sin(q3 - q4) +
			1.0*dotq2*sin(q4) + 15.0*dotq3*sin(q4) - 7.5*dotq4*sin(q4) -
			15.0*dotq3*sin(q3 + q4) - 7.5*dotq4*sin(q3 + q4) +
			1.2*dotq3*sin(q3 - q4 - q5) + 0.2*dotq1*sin(q4 + q5) -
			0.6*dotq3*sin(q4 + q5) + 0.6*dotq4*sin(q4 + q5) +
			1.2*dotq5*sin(q4 + q5) + 0.6*dotq3*sin(q3 + q4 + q5) +
			0.3*dotq4*sin(q3 + q4 + q5) + 0.3*dotq5*sin(q3 + q4 + q5) -
			0.6*dotq5*cos(q3)*sin(q3 + q4 + q5) -
			1.5*dotq5*sin(q3)*sin(q3 + q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(372.813 + 450.0*cos(q3) +
			562.5*cos(q4) + 562.5*cos(q3 + q4) - 45.0*cos(q4 + q5) -
			22.5*cos(q3 + q4 + q5) + 225.0*sin(q4) + 112.5*sin(q3 + q4) -
			90.0*sin(q3 - q5) + 90.0*sin(q3 + q5) + 112.5*sin(q4 + q5) +
			112.5*sin(
			q3 + q4 + q5))*(dotq5*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			dotq4*(sin(
			q3)*(2.0 + 2.5*cos(q5)*sin(q4) + 2.0*sin(q5) +
			cos(q4)*(2.5 + 2.5*sin(q5))) +
			cos(q3)*(-2.5*cos(q4)*cos(q5) +
			sin(q4)*(2.5 + 2.5*sin(q5)))) +
			1.0*(1.0*cos(q3 + q4 + q5) -
			1.0*sin(q3 + q4))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 3)*(3.0*cos(q3 + q4) -
			4.0*cos(q3 - q5) + 4.0*cos(q3 + q5) + 10.0*cos(q3 + q4 + q5) -
			8.0*sin(q3) - 10.0*sin(q3 + q4) +
			3.0*sin(q3 + q4 +
			q5))*(149.625*(dotq5*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			dotq4*(sin(
			q3)*(2.0 + 2.5*cos(q5)*sin(q4) + 2.0*sin(q5) +
			cos(q4)*(2.5 + 2.5*sin(q5))) +
			cos(q3)*(-2.5*cos(q4)*cos(q5) +
			sin(q4)*(2.5 + 2.5*sin(q5)))) +
			1.0*(1.0*cos(q3 + q4 + q5) -
			1.0*sin(q3 + q4))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			149.625*(dotq5*(cos(
			q3)*(-10.0 - 25.0*cos(q5)*sin(q4) +
			cos(q4)*(7.5*cos(q5) - 25.0*sin(q5)) - 10.0*sin(q5) -
			7.5*sin(q4)*sin(q5)) +
			sin(q3)*(-3.0 +
			cos(q5)*(-20.0 - 25.0*cos(q4) - 7.5*sin(q4)) - 
			3.0*sin(q5) - 7.5*cos(q4)*sin(q5) +
			sin(q4)*(-8.0 + 17.0*sin(q5)))) +
			8.0*(1.0*cos(q3) + 1.25*cos(q3 + q4) -
			0.375*cos(q3 + q4 + q5) + 0.375*sin(q3 + q4) -
			0.5*sin(q3 - q5) + 0.5*sin(q3 + q5) +
			1.25*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5)) +
			dotq4*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(-10.0*cos(q3 + q4) +
			3.0*cos(q3 + q4 + q5) - 3.0*sin(q3 + q4) -
			10.0*sin(q3 + q4 + q5)) -
			3.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(1.0*cos(q3 + q4) -
			1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
			3.33333*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
			3.33333*sin(q3 + q4) + 1.0*sin(q3 + q4 + q5))))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*(75.0*cos(q3 + q4) +
			15.0*sin(q3 + q4 + q5))*(-7.5*dotq5*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			18.75*dotq4*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(cos(
			q3)*(1.0*cos(q4)*cos(q5) + sin(q4)*(-1.0 - 1.0*sin(q5))) +
			sin(q3)*(-0.8 - 1.0*cos(q5)*sin(q4) +
			cos(q4)*(-1.0 - 1.0*sin(q5)) - 0.8*sin(q5))) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(1.0*cos(2*q3 + q4) +
			2.5*cos(q3 + q4 + q5) - 2.5*cos(2*q3 + q4 + q5) -
			2.5*sin(q3 + q4) + 2.5*sin(2*q3 + q4) +
			1.0*sin(2*q3 + q4 + q5))) + (15.0*cos(q3 + q4 + q5) -
			75.0*sin(q3 + q4))*(dotq5*(6*cos(q4 + q5)*sin(q3) +
			6*cos(q3)*sin(q4 + q5) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(cos(q3)*cos(q4 + q5) -
			1.0*sin(q3)*sin(q4 + q5)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) + (
			3.0*cos(q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*- (
			6.0*sin(q4 + q5)*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			dotq4*(6*sin(q3)*(cos(q4 + q5) - sin(q4)) +
			6*cos(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(1.0*cos(q3 + q4 + q5) - 1.0*sin(q3 + q4)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*3.0*(1.0*cos(q4 + q5) -
			1.0*sin(q4))*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))) - (
			3.0*(2.0*cos(q4) +
			2.0*sin(q4 + q5))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) + 
			1.0*sin(q4 + q5), 2)*7.5*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(-1.6*cos(q3) - 2.0*cos(q3 + q4) +
			1.0*cos(2*q3 + q4) + 0.2*cos(q3 + q4 + q5) +
			0.4*cos(2*q3 + q4 + q5) - 0.2*sin(q3 + q4) -
			0.4*sin(2*q3 + q4) + 0.8*sin(q3 - q5) - 0.8*sin(q3 + q5) -
			2.0*sin(q3 + q4 + q5) + 1.0*sin(2*q3 + q4 + q5))) +
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*7.5*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5)))*(15.0*dotq5*cos(q3 + q4 + q5) +
			dotq3*(15.0*cos(q3 + q4 + q5) - 75.0*sin(q3 + q4)) +
			dotq4*(15.0*cos(q3 + q4 + q5) - 75.0*sin(q3 + q4)) + (
			1215.0*(1.0*dotq3*cos(
			q3) + (-0.0123457*dotq4 - 0.0123457*dotq5)*cos(q4 + q5) +
			2.5*dotq3*sin(q3) +
			0.0617284*dotq4*sin(q4))*(1.0*cos(q3 + q4) +
			1.0*sin(q3 + q4 + q5)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(-847.5 +
			3037.5*cos(q3) + 75.0*cos(q4) - 1215.0*sin(q3) +
			15.0*sin(q4 + q5))*(dotq5*(0.5*cos(q3 - q5) -
			0.5*cos(q3 + q5) - 2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			dotq4*(sin(
			q3)*(2.0 + 2.5*cos(q5)*sin(q4) + 2.0*sin(q5) +
			cos(q4)*(2.5 + 2.5*sin(q5))) +
			cos(q3)*(-2.5*cos(q4)*cos(q5) +
			sin(q4)*(2.5 + 2.5*sin(q5)))) +
			1.0*(1.0*cos(q3 + q4 + q5) -
			1.0*sin(q3 + q4))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 +
			q5), 2)*405.0*(-7.5*dotq5*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			18.75*dotq4*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(cos(
			q3)*(1.0*cos(q4)*cos(q5) + sin(q4)*(-1.0 - 1.0*sin(q5))) +
			sin(q3)*(-0.8 - 1.0*cos(q5)*sin(q4) +
			cos(q4)*(-1.0 - 1.0*sin(q5)) - 0.8*sin(q5))) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(1.0*cos(2*q3 + q4) +
			2.5*cos(q3 + q4 + q5) - 2.5*cos(2*q3 + q4 + q5) -
			2.5*sin(q3 + q4) + 2.5*sin(2*q3 + q4) +
			1.0*sin(2*q3 + q4 + q5)))) + (-6*cos(
			q3)*(cos(q4 + q5) - sin(q4)) +
			6*sin(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)))*(dotq3*(-75.0*cos(q3 + q4) -
			15.0*sin(q3 + q4 + q5)) +
			dotq4*(-75.0*cos(q3 + q4) - 15.0*sin(q3 + q4 + q5)) -
			15.0*dotq5*sin(q3 + q4 + q5) + (
			75.0*(40.5*dotq3*cos(q3) + 1.0*dotq4*cos(q4) -
			16.2*dotq3*sin(q3) + 0.2*dotq4*sin(q4 + q5) +
			0.2*dotq5*sin(q4 + q5))*(1.0*cos(q3 + q4) + 
			1.0*sin(q3 + q4 + q5)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(369.0 -
			1215.0*cos(q3) + 15.0*cos(q4 + q5) - 3037.5*sin(q3) -
			75.0*sin(q4))*(dotq5*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			dotq4*(sin(
			q3)*(2.0 + 2.5*cos(q5)*sin(q4) + 2.0*sin(q5) +
			cos(q4)*(2.5 + 2.5*sin(q5))) +
			cos(q3)*(-2.5*cos(q4)*cos(q5) +
			sin(q4)*(2.5 + 2.5*sin(q5)))) +
			1.0*(1.0*cos(q3 + q4 + q5) -
			1.0*sin(q3 + q4))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			405.0*(dotq5*(6*cos(q4 + q5)*sin(q3) + 6*cos(q3)*sin(q4 + q5) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(cos(q3)*cos(q4 + q5) -
			1.0*sin(q3)*sin(q4 + q5)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) + (
			3.0*cos(q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*- (
			6.0*sin(q4 + q5)*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			dotq4*(6*sin(q3)*(cos(q4 + q5) - sin(q4)) +
			6*cos(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(1.0*cos(q3 + q4 + q5) - 1.0*sin(q3 + q4)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*3.0*(1.0*cos(q4 + q5) -
			1.0*sin(q4))*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))) - (
			3.0*(2.0*cos(q4) +
			2.0*sin(q4 + q5))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*7.5*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(-1.6*cos(q3) - 2.0*cos(q3 + q4) +
			1.0*cos(2*q3 + q4) + 0.2*cos(q3 + q4 + q5) +
			0.4*cos(2*q3 + q4 + q5) - 0.2*sin(q3 + q4) -
			0.4*sin(2*q3 + q4) + 0.8*sin(q3 - q5) -
			0.8*sin(q3 + q5) - 2.0*sin(q3 + q4 + q5) +
			1.0*sin(2*q3 + q4 + q5)))) +
			1/(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(-1.0*cos(q3 + q4) -
			1.0*sin(q3 + q4 + q5))*(37.5*dotq1*cos(q3 + q4) -
			7.5*dotq2*cos(q3 + q4 + q5) + 37.5*dotq2*sin(q3 + q4) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*4556.25*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(1.0*dotq3*cos(
			q3) + (-0.0246914*dotq4 - 0.0246914*dotq5)*cos(q4 + q5) +
			0.0123457*dotq4*cos(q3 + q4 + q5) + 
			0.0123457*dotq5*cos(q3 + q4 + q5) + 2.5*dotq3*sin(q3) +
			0.123457*dotq4*sin(q4) -
			0.0617284*dotq4*sin(q3 + q4))*(sin(
			q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))) +
			56.25*dotq3*(1.0*cos(q3 + q4) - 0.8*cos(q3 - q5) +
			0.8*cos(q3 + q5) + 1.0*cos(q3 + q4 + q5) - 4.0*sin(q3) -
			5.0*sin(q3 + q4) + 0.2*sin(q3 + q4 + q5)) +
			7.5*dotq1*sin(q3 + q4 + q5) +
			dotq5*(90.0*cos(q3 - q5) + 90.0*cos(q3 + q5) +
			112.5*cos(q4 + q5) + 112.5*cos(q3 + q4 + q5) +
			45.0*sin(q4 + q5) + 22.5*sin(q3 + q4 + q5)) +
			dotq4*(225.0*cos(q4) + 112.5*cos(q3 + q4) + 112.5*cos(q4 + q5) +
			112.5*cos(q3 + q4 + q5) - 562.5*sin(q4) -
			562.5*sin(q3 + q4) + 45.0*sin(q4 + q5) +
			22.5*sin(q3 + q4 + q5)) +
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*56.25*(1.0*cos(q3 + q4) +
			1.0*sin(q3 + q4 + q5))*((-27.0*dotq1 - 10.8*dotq2 +
			4.0*dotq3)*cos(q3) + (4.0*dotq3 - 8.0*dotq4)*cos(q3 - q4) -
			4.0*dotq4*cos(q4) + 1.0*dotq4*cos(q3 + q4) -
			0.8*dotq4*cos(q3 - q5) + 2.0*dotq3*cos(q3 - q4 - q5) -
			4.0*dotq4*cos(q3 - q4 - q5) - 4.0*dotq5*cos(q3 - q4 - q5) -
			3.2*dotq5*cos(q5) + 0.8*dotq4*cos(q3 + q5) +
			0.8*dotq5*cos(q3 + q5) - 4.0*dotq4*cos(q4 + q5) -
			4.0*dotq5*cos(q4 + q5) + 1.0*dotq4*cos(q3 + q4 + q5) +
			1.0*dotq5*cos(q3 + q4 + q5) + 10.8*dotq1*sin(q3) -
			27.0*dotq2*sin(q3) - 132.68*dotq3*sin(q3) -
			4.0*dotq4*sin(q3) - 0.4*dotq5*sin(q3) +
			10.0*dotq3*sin(q3 - q4) - 20.0*dotq4*sin(q3 - q4) +
			20.0*dotq4*sin(q4) - 5.0*dotq4*sin(q3 + q4) -
			0.8*dotq3*sin(q3 - q4 - q5) + 1.6*dotq4*sin(q3 - q4 - q5) +
			1.6*dotq5*sin(q3 - q4 - q5) - 0.8*dotq4*sin(q4 + q5) -
			0.8*dotq5*sin(q4 + q5) + 0.2*dotq4*sin(q3 + q4 + q5) +
			0.2*dotq5*sin(q3 + q4 + q5)) + (-6*cos(
			q3)*(cos(q4 + q5) - sin(q4)) +
			6*sin(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)))*(-1518.75*dotq3*cos(q3) -
			75.0*dotq4*cos(q4) + 37.5*dotq4*cos(q3 + q4) +
			607.5*dotq3*sin(q3) - 15.0*dotq4*sin(q4 + q5) -
			15.0*dotq5*sin(q4 + q5) + 7.5*dotq4*sin(q3 + q4 + q5) +
			7.5*dotq5*sin(q3 + q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(38416.9 -
			14926.5*cos(q3) + 1125.0*cos(q3 - q4) + 1125.0*cos(q4) -
			90.0*cos(q3 - q4 - q5) - 45.0*cos(q4 + q5) - 450.0*sin(q3) -
			450.0*sin(q3 - q4) + 225.0*sin(q4) - 225.0*sin(q3 - q4 - q5) +
			180.0*sin(q5) +
			225.0*sin(
			q4 + q5))*(dotq5*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			dotq4*(sin(
			q3)*(2.0 + 2.5*cos(q5)*sin(q4) + 2.0*sin(q5) +
			cos(q4)*(2.5 + 2.5*sin(q5))) +
			cos(q3)*(-2.5*cos(q4)*cos(q5) +
			sin(q4)*(2.5 + 2.5*sin(q5)))) +
			1.0*(1.0*cos(q3 + q4 + q5) - 
			1.0*sin(q3 + q4))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 +
			q5), 2)*149.625*(dotq5*(cos(
			q3)*(-10.0 - 25.0*cos(q5)*sin(q4) +
			cos(q4)*(7.5*cos(q5) - 25.0*sin(q5)) - 10.0*sin(q5) -
			7.5*sin(q4)*sin(q5)) +
			sin(q3)*(-3.0 +
			cos(q5)*(-20.0 - 25.0*cos(q4) - 7.5*sin(q4)) -
			3.0*sin(q5) - 7.5*cos(q4)*sin(q5) +
			sin(q4)*(-8.0 + 17.0*sin(q5)))) +
			8.0*(1.0*cos(q3) + 1.25*cos(q3 + q4) -
			0.375*cos(q3 + q4 + q5) + 0.375*sin(q3 + q4) -
			0.5*sin(q3 - q5) + 0.5*sin(q3 + q5) +
			1.25*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5)) +
			dotq4*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(-10.0*cos(q3 + q4) +
			3.0*cos(q3 + q4 + q5) - 3.0*sin(q3 + q4) -
			10.0*sin(q3 + q4 + q5)) -
			3.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(1.0*cos(q3 + q4) -
			1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
			3.33333*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
			3.33333*sin(q3 + q4) + 1.0*sin(q3 + q4 + q5)))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*(-847.5 + 3037.5*cos(q3) + 75.0*cos(q4) -
			1215.0*sin(q3) +
			15.0*sin(q4 + q5))*(-7.5*dotq5*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			18.75*dotq4*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(cos(
			q3)*(1.0*cos(q4)*cos(q5) + sin(q4)*(-1.0 - 1.0*sin(q5))) +
			sin(q3)*(-0.8 - 1.0*cos(q5)*sin(q4) +
			cos(q4)*(-1.0 - 1.0*sin(q5)) - 0.8*sin(q5))) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(1.0*cos(2*q3 + q4) +
			2.5*cos(q3 + q4 + q5) - 2.5*cos(2*q3 + q4 + q5) -
			2.5*sin(q3 + q4) + 2.5*sin(2*q3 + q4) +
			1.0*sin(2*q3 + q4 + q5))) + (369.0 - 1215.0*cos(q3) +
			15.0*cos(q4 + q5) - 3037.5*sin(q3) -
			75.0*sin(
			q4))*(dotq5*(6*cos(q4 + q5)*sin(q3) +
			6*cos(q3)*sin(q4 + q5) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(cos(q3)*cos(q4 + q5) -
			1.0*sin(q3)*sin(q4 + q5)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) + (
			3.0*cos(q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*- (
			6.0*sin(q4 + q5)*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			dotq4*(6*sin(q3)*(cos(q4 + q5) - sin(q4)) +
			6*cos(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) + 
			2.0*sin(q4))*(1.0*cos(q3 + q4 + q5) - 1.0*sin(q3 + q4)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*3.0*(1.0*cos(q4 + q5) -
			1.0*sin(q4))*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))) - (
			3.0*(2.0*cos(q4) +
			2.0*sin(q4 + q5))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*7.5*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(-1.6*cos(q3) - 2.0*cos(q3 + q4) +
			1.0*cos(2*q3 + q4) + 0.2*cos(q3 + q4 + q5) +
			0.4*cos(2*q3 + q4 + q5) - 0.2*sin(q3 + q4) -
			0.4*sin(2*q3 + q4) + 0.8*sin(q3 - q5) -
			0.8*sin(q3 + q5) - 2.0*sin(q3 + q4 + q5) +
			1.0*sin(2*q3 + q4 + q5))));
	C2[0][0] = value;

	//C2(0,1)
	value = 90*dotq5*cos(q5) - 7.5*dotq2*cos(q3 + q4 + q5) +
			dotq3*(-56.25*cos(q4 + q5) - 56.25*cos(q3 + q4 + q5) -
			22.5*sin(q4 + q5) - 11.25*sin(q3 + q4 + q5)) +
			7.5*dotq1*sin(q3 + q4 + q5) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))*56.25*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(1.0*dotq3*cos(
			q4 + q5) + (-2.0*dotq3 - 1.0*dotq4 - 1.0*dotq5)*cos(
			q3 + q4 + q5) - 5.0*dotq3*sin(q4) + 10.0*dotq3*sin(q3 + q4) +
			5.0*dotq4*sin(q3 + q4))*sin(q3 + q4 + q5) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*22.5*sin(
			q3 + q4 +
			q5)*(-10.0*dotq3*cos(
			q3 - q4) + (1.66667*dotq1 - 5.0*dotq3 + 5.0*dotq4)*cos(q4) +
			5.0*dotq3*cos(q3 + q4) + 2.5*dotq4*cos(q3 + q4) -
			4.0*dotq3*cos(q3 - q5) + 4.0*dotq5*cos(q3 - q5) -
			5.0*dotq3*cos(q3 - q4 - q5) + 4.0*dotq3*cos(q3 + q5) +
			4.0*dotq5*cos(q3 + q5) - 0.333333*dotq2*cos(q4 + q5) -
			5.0*dotq3*cos(q4 + q5) + 2.5*dotq4*cos(q4 + q5) +
			2.5*dotq5*cos(q4 + q5) + 5.0*dotq3*cos(q3 + q4 + q5) +
			2.5*dotq4*cos(q3 + q4 + q5) + 2.5*dotq5*cos(q3 + q4 + q5) -
			20.0*dotq3*sin(q3) - 25.0*dotq3*sin(q3 - q4) +
			1.66667*dotq2*sin(q4) + 25.0*dotq3*sin(q4) -
			12.5*dotq4*sin(q4) - 25.0*dotq3*sin(q3 + q4) -
			12.5*dotq4*sin(q3 + q4) + 2.0*dotq3*sin(q3 - q4 - q5) +
			0.333333*dotq1*sin(q4 + q5) - 1.0*dotq3*sin(q4 + q5) +
			1.0*dotq4*sin(q4 + q5) + 1.0*dotq5*sin(q4 + q5) +
			1.0*dotq3*sin(q3 + q4 + q5) + 0.5*dotq4*sin(q3 + q4 + q5) +
			0.5*dotq5*sin(q3 + q4 + q5)) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*225.0*(1.0*dotq3*cos(q4) + (-2.0*dotq3 - 1.0*dotq4)*cos(q3 + q4) +
			0.2*dotq3*sin(q4 + q5) - 0.4*dotq3*sin(q3 + q4 + q5) -
			0.2*dotq4*sin(q3 + q4 + q5) -
			0.2*dotq5*sin(q3 + q4 + q5))*(cos(
			q3 + q4 + q5)*(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) + (0.25 +
			0.5*cos(q3) - 1.0*cos(q4 + q5) + 1.25*sin(q3) + 
			1.0*sin(q4))*sin(q3 + q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(372.813 + 450.0*cos(q3) +
			562.5*cos(q4) + 562.5*cos(q3 + q4) - 45.0*cos(q4 + q5) -
			22.5*cos(q3 + q4 + q5) + 225.0*sin(q4) + 112.5*sin(q3 + q4) -
			90.0*sin(q3 - q5) + 90.0*sin(q3 + q5) + 112.5*sin(q4 + q5) +
			112.5*sin(q3 + q4 + q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
			q3 + q5) + (-2.5*dotq4 - 2.5*dotq5 +
			1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
			0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
			1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
			0.5*dotq5*sin(2*(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 3)*(3.0*cos(q3 + q4) -
			4.0*cos(q3 - q5) + 4.0*cos(q3 + q5) + 10.0*cos(q3 + q4 + q5) -
			8.0*sin(q3) - 10.0*sin(q3 + q4) +
			3.0*sin(q3 + q4 +
			q5))*(149.625*(dotq5*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(3.0*cos(q3 + q4 + q5) -
			4.0*sin(q3 + q5) - 10.0*sin(q3 + q4 + q5)) -
			4.0*cos(q4 + q5)*(1.0*cos(q3 + q5) +
			2.5*cos(q3 + q4 + q5) - 1.0*sin(q3) +
			0.75*sin(q3 + q4 + q5))) +
			dotq4*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(3.0*cos(q3 + q4 + q5) -
			10.0*sin(q3 + q4 + q5)) -
			4.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(1.0*cos(q3 + q5) +
			2.5*cos(q3 + q4 + q5) - 1.0*sin(q3) +
			0.75*sin(q3 + q4 + q5))) +
			4.0*(1.0*cos(q3) - 0.75*cos(q3 + q4 + q5) + 1.0*sin(q3 + q5) +
			2.5*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			149.625*((-1.0*dotq4 - 0.5*dotq5)*cos(
			q3 + q5) + (-2.5*dotq4 - 2.5*dotq5 +
			1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
			0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
			1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
			0.5*dotq5*sin(2*(q3 + q4 + q5)))) + (-6*cos(
			q3)*(cos(q4 + q5) - sin(q4)) +
			6*sin(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)))*(-15.0*dotq3*sin(q3 + q4 + q5) -
			15.0*dotq4*sin(q3 + q4 + q5) - 15.0*dotq5*sin(q3 + q4 + q5) + (
			75.0*(40.5*dotq3*cos(q3) + 1.0*dotq4*cos(q4) -
			16.2*dotq3*sin(q3) + 0.2*dotq4*sin(q4 + q5) +
			0.2*dotq5*sin(q4 + q5))*sin(q3 + q4 + q5))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(369.0 -
			1215.0*cos(q3) + 15.0*cos(q4 + q5) - 3037.5*sin(q3) -
			75.0*sin(q4))*((-1.0*dotq4 - 0.5*dotq5)*cos(
			q3 + q5) + (-2.5*dotq4 - 2.5*dotq5 +
			1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
			0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
			1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
			0.5*dotq5*sin(2*(q3 + q4 + q5))) -
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 +
			q5), 2)*1215.0*(1.0*dotq5*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) +
			2.0*cos(q4 + q5)*(-0.25 - 0.5*cos(q3) + 1.0*cos(q4 + q5) -
			1.25*sin(q3) - 1.0*sin(q4))*sin(q3 + q4 + q5) +
			2.0*sin(q4 + q5)*(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))*sin(
			q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			1.0*dotq4*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) +
			2.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(-0.25 - 0.5*cos(q3) +
			1.0*cos(q4 + q5) - 1.25*sin(q3) - 1.0*sin(q4))*sin(
			q3 + q4 + q5) +
			2.0*(1.0*cos(q4) + 1.0*sin(q4 + q5))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*sin(q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			2.5*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(0.8*cos(q3) - 0.2*cos(q3 + q4 + q5) -
			0.4*cos(2*q3 + q4 + q5) + 0.8*sin(q3 + q5) +
			2.0*sin(q3 + q4 + q5) - 1.0*sin(2*q3 + q4 + q5)))) -
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*3.0*(15.0*cos(q3 + q4 + q5) -
			75.0*sin(q3 +
			q4))*(1.0*dotq5*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) +
			2.0*cos(q4 + q5)*(-0.25 - 0.5*cos(q3) + 1.0*cos(q4 + q5) -
			1.25*sin(q3) - 1.0*sin(q4))*sin(q3 + q4 + q5) +
			2.0*sin(q4 + q5)*(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))*sin(
			q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			1.0*dotq4*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) +
			2.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(-0.25 - 0.5*cos(q3) +
			1.0*cos(q4 + q5) - 1.25*sin(q3) - 1.0*sin(q4))*sin(
			q3 + q4 + q5) +
			2.0*(1.0*cos(q4) + 1.0*sin(q4 + q5))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*sin(q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			2.5*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(0.8*cos(q3) - 0.2*cos(q3 + q4 + q5) -
			0.4*cos(2*q3 + q4 + q5) + 0.8*sin(q3 + q5) +
			2.0*sin(q3 + q4 + q5) - 1.0*sin(2*q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(75.0*cos(q3 + q4) +
			15.0*sin(q3 + q4 + q5))*(-1.5*dotq5*(1.0*cos(q3 + q5) +
			5.0*cos(q3 + q4 + q5) + 1.0*cos(q3 + 2*q4 + q5) -
			2.0*sin(q3))*(2.5 - 2.5*cos(q3) + 1.0*sin(q3)) -
			3.0*dotq4*(1.0*cos(q3 + q5) + 2.5*cos(q3 + q4 + q5) -
			1.0*sin(q3))*(2.5 - 2.5*cos(q3) + 1.0*sin(q3)) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(2.5*cos(q3 + q4 + q5) -
			2.5*cos(2*q3 + q4 + q5) + 1.0*sin(2*q3 + q4 + q5))) + 
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*7.5*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5)))*(15.0*dotq3*cos(
			q3 + q4 + q5) + 15.0*dotq4*cos(q3 + q4 + q5) +
			15.0*dotq5*cos(q3 + q4 + q5) + (
			1215.0*(1.0*dotq3*cos(
			q3) + (-0.0123457*dotq4 - 0.0123457*dotq5)*cos(q4 + q5) +
			2.5*dotq3*sin(q3) + 0.0617284*dotq4*sin(q4))*sin(
			q3 + q4 + q5))/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*(-847.5 + 3037.5*cos(q3) + 75.0*cos(q4) -
			1215.0*sin(q3) +
			15.0*sin(q4 + q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
			q3 + q5) + (-2.5*dotq4 - 2.5*dotq5 +
			1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
			0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
			1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
			0.5*dotq5*sin(2*(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*405.0*(-1.5*dotq5*(1.0*cos(q3 + q5) +
			5.0*cos(q3 + q4 + q5) + 1.0*cos(q3 + 2*q4 + q5) -
			2.0*sin(q3))*(2.5 - 2.5*cos(q3) + 1.0*sin(q3)) -
			3.0*dotq4*(1.0*cos(q3 + q5) + 2.5*cos(q3 + q4 + q5) -
			1.0*sin(q3))*(2.5 - 2.5*cos(q3) + 1.0*sin(q3)) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(2.5*cos(q3 + q4 + q5) -
			2.5*cos(2*q3 + q4 + q5) + 1.0*sin(2*q3 + q4 + q5)))) +
			1/(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(-1.0*cos(q3 + q4) -
			1.0*sin(q3 + q4 + q5))*(-7.5*dotq2*cos(q3 + q4 + q5) +
			7.5*dotq1*sin(q3 + q4 + q5) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*4556.25*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(1.0*dotq3*cos(
			q3) + (-0.0246914*dotq4 - 0.0246914*dotq5)*cos(q4 + q5) +
			0.0123457*dotq4*cos(q3 + q4 + q5) +
			0.0123457*dotq5*cos(q3 + q4 + q5) + 2.5*dotq3*sin(q3) +
			0.123457*dotq4*sin(q4) - 0.0617284*dotq4*sin(q3 + q4))*sin(
			q3 + q4 + q5) +
			dotq3*(45.0*cos(q3 + q5) + 56.25*cos(q3 + q4 + q5) -
			22.5*sin(q3) + 11.25*sin(q3 + q4 + q5)) +
			dotq4*(112.5*cos(q4 + q5) + 112.5*cos(q3 + q4 + q5) +
			45.0*sin(q4 + q5) + 22.5*sin(q3 + q4 + q5)) +
			dotq5*(90.0*cos(q3 + q5) + 112.5*cos(q4 + q5) +
			112.5*cos(q3 + q4 + q5) + 45.0*sin(q4 + q5) +
			22.5*sin(q3 + q4 + q5)) +
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*11.25*sin(
			q3 + q4 +
			q5)*((-135.0*dotq1 - 54.0*dotq2 + 20.0*dotq3)*cos(
			q3) + (20.0*dotq3 - 40.0*dotq4)*cos(q3 - q4) -
			20.0*dotq4*cos(q4) + 5.0*dotq4*cos(q3 + q4) -
			4.0*dotq4*cos(q3 - q5) + 10.0*dotq3*cos(q3 - q4 - q5) -
			20.0*dotq4*cos(q3 - q4 - q5) - 20.0*dotq5*cos(q3 - q4 - q5) -
			16.0*dotq5*cos(q5) + 4.0*dotq4*cos(q3 + q5) +
			4.0*dotq5*cos(q3 + q5) - 20.0*dotq4*cos(q4 + q5) -
			20.0*dotq5*cos(q4 + q5) + 5.0*dotq4*cos(q3 + q4 + q5) +
			5.0*dotq5*cos(q3 + q4 + q5) + 54.0*dotq1*sin(q3) -
			135.0*dotq2*sin(q3) - 663.4*dotq3*sin(q3) - 
			20.0*dotq4*sin(q3) - 2.0*dotq5*sin(q3) +
			50.0*dotq3*sin(q3 - q4) - 100.0*dotq4*sin(q3 - q4) +
			100.0*dotq4*sin(q4) - 25.0*dotq4*sin(q3 + q4) -
			4.0*dotq3*sin(q3 - q4 - q5) + 8.0*dotq4*sin(q3 - q4 - q5) +
			8.0*dotq5*sin(q3 - q4 - q5) - 4.0*dotq4*sin(q4 + q5) -
			4.0*dotq5*sin(q4 + q5) + 1.0*dotq4*sin(q3 + q4 + q5) +
			1.0*dotq5*sin(q3 + q4 + q5)) +
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*4556.25*(1.0*dotq3*cos(q3) + 0.0493827*dotq4*cos(q4) -
			0.0246914*dotq4*cos(q3 + q4) - 0.4*dotq3*sin(q3) +
			0.00987654*dotq4*sin(q4 + q5) +
			0.00987654*dotq5*sin(q4 + q5) -
			0.00493827*dotq4*sin(q3 + q4 + q5) -
			0.00493827*dotq5*sin(q3 + q4 + q5))*(cos(
			q3 + q4 + q5)*(5.0 + 2.0*cos(q4) + 2.0*sin(q4 + q5)) + (0.5 +
			1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*sin(q3 + q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 +
			q5), 2)*149.625*(dotq5*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(3.0*cos(q3 + q4 + q5) -
			4.0*sin(q3 + q5) - 10.0*sin(q3 + q4 + q5)) -
			4.0*cos(q4 + q5)*(1.0*cos(q3 + q5) +
			2.5*cos(q3 + q4 + q5) - 1.0*sin(q3) +
			0.75*sin(q3 + q4 + q5))) +
			dotq4*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(3.0*cos(q3 + q4 + q5) -
			10.0*sin(q3 + q4 + q5)) -
			4.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(1.0*cos(q3 + q5) +
			2.5*cos(q3 + q4 + q5) - 1.0*sin(q3) +
			0.75*sin(q3 + q4 + q5))) +
			4.0*(1.0*cos(q3) - 0.75*cos(q3 + q4 + q5) + 1.0*sin(q3 + q5) +
			2.5*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(38416.9 -
			14926.5*cos(q3) + 1125.0*cos(q3 - q4) + 1125.0*cos(q4) -
			90.0*cos(q3 - q4 - q5) - 45.0*cos(q4 + q5) - 450.0*sin(q3) -
			450.0*sin(q3 - q4) + 225.0*sin(q4) - 225.0*sin(q3 - q4 - q5) +
			180.0*sin(q5) +
			225.0*sin(q4 + q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
			q3 + q5) + (-2.5*dotq4 - 2.5*dotq5 +
			1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
			0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
			1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
			0.5*dotq5*sin(2*(q3 + q4 + q5))) -
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*3.0*(369.0 -
			1215.0*cos(q3) + 15.0*cos(q4 + q5) - 3037.5*sin(q3) -
			75.0*sin(q4))*(1.0*dotq5*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) +
			2.0*cos(q4 + q5)*(-0.25 - 0.5*cos(q3) + 1.0*cos(q4 + q5) -
			1.25*sin(q3) - 1.0*sin(q4))*sin(q3 + q4 + q5) +
			2.0*sin(q4 + q5)*(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))*sin(
			q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			1.0*dotq4*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) +
			2.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(-0.25 - 0.5*cos(q3) +
			1.0*cos(q4 + q5) - 1.25*sin(q3) - 1.0*sin(q4))*sin(
			q3 + q4 + q5) +
			2.0*(1.0*cos(q4) + 1.0*sin(q4 + q5))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*sin(q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			2.5*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(0.8*cos(q3) - 0.2*cos(q3 + q4 + q5) -
			0.4*cos(2*q3 + q4 + q5) + 0.8*sin(q3 + q5) +
			2.0*sin(q3 + q4 + q5) - 1.0*sin(2*q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(-847.5 +
			3037.5*cos(q3) + 75.0*cos(q4) - 1215.0*sin(q3) +
			15.0*sin(q4 + q5))*(-1.5*dotq5*(1.0*cos(q3 + q5) +
			5.0*cos(q3 + q4 + q5) + 1.0*cos(q3 + 2*q4 + q5) -
			2.0*sin(q3))*(2.5 - 2.5*cos(q3) + 1.0*sin(q3)) -
			3.0*dotq4*(1.0*cos(q3 + q5) + 2.5*cos(q3 + q4 + q5) -
			1.0*sin(q3))*(2.5 - 2.5*cos(q3) + 1.0*sin(q3)) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(2.5*cos(q3 + q4 + q5) -
			2.5*cos(2*q3 + q4 + q5) +
			1.0*sin(2*q3 + q4 + q5))));
 	C2[0][1] = value;

	//C2(1,0)
	value = -90.0*dotq4*cos(q5) +
			45*dotq5*cos(q5) - 7.5*dotq2*cos(q3 + q4 + q5) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))*56.25*(1.0*dotq3*cos(
			q4 + q5) + (-2.0*dotq3 - 1.0*dotq4 - 1.0*dotq5)*cos(
			q3 + q4 + q5))*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))) -
			45.0*dotq3*(1.0*cos(q3 - q5) + 1.0*cos(q3 + q5) + 1.25*cos(q4 + q5) +
			1.25*cos(q3 + q4 + q5) + 0.5*sin(q4 + q5) +
			0.25*sin(q3 + q4 + q5)) +
			7.5*dotq1*sin(
			q3 + q4 + q5) + (-6*cos(q3)*(cos(q4 + q5) - sin(q4)) +
			6*sin(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)))*(7.5*dotq3*sin(
			q4 + q5) + (-15.0*dotq3 - 7.5*dotq4 - 7.5*dotq5)*sin(
			q3 + q4 + q5)) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*45.0*(1.0*cos(q3 + q4) +
			1.0*sin(q3 + q4 + q5))*(-1.0*dotq4*cos(q3 - q5) -
			2.5*dotq3*cos(q3 - q4 - q5) - 2.0*dotq3*cos(q5) +
			2.0*dotq3*cos(q3 + q5) - 1.0*dotq4*cos(q3 + q5) +
			1.0*dotq5*cos(q3 + q5) - 0.166667*dotq2*cos(q4 + q5) -
			2.5*dotq3*cos(q4 + q5) + 1.25*dotq4*cos(q4 + q5) +
			1.25*dotq5*cos(q4 + q5) + 2.5*dotq3*cos(q3 + q4 + q5) +
			1.25*dotq4*cos(q3 + q4 + q5) + 1.25*dotq5*cos(q3 + q4 + q5) -
			1.0*dotq3*sin(q3) + 1.0*dotq3*sin(q3 - q4 - q5) +
			0.166667*dotq1*sin(q4 + q5) - 0.5*dotq3*sin(q4 + q5) +
			0.5*dotq4*sin(q4 + q5) + 0.5*dotq5*sin(q4 + q5) +
			0.5*dotq3*sin(q3 + q4 + q5) + 0.25*dotq4*sin(q3 + q4 + q5) +
			0.25*dotq5*sin(q3 + q4 + q5)) + 
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(186.406 +
			cos(q3 + q4 + q5)*(-22.5 - 45.0*cos(q3) + 45.0*cos(q4 + q5) -
			112.5*sin(q3) - 90.0*sin(q4)) + (112.5 + 112.5*cos(q3) +
			90.0*cos(q4) - 45.0*sin(q3) + 45.0*sin(q4 + q5))*sin(
			q3 + q4 + q5))*(dotq5*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			dotq4*(sin(
			q3)*(2.0 + 2.5*cos(q5)*sin(q4) + 2.0*sin(q5) +
			cos(q4)*(2.5 + 2.5*sin(q5))) +
			cos(q3)*(-2.5*cos(q4)*cos(q5) +
			sin(q4)*(2.5 + 2.5*sin(q5)))) +
			1.0*(1.0*cos(q3 + q4 + q5) -
			1.0*sin(q3 + q4))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 3)*(4.0*cos(q3 + q5) +
			10.0*cos(q3 + q4 + q5) - 4.0*sin(q3) +
			3.0*sin(q3 + q4 +
			q5))*(149.625*(dotq5*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			dotq4*(sin(
			q3)*(2.0 + 2.5*cos(q5)*sin(q4) + 2.0*sin(q5) +
			cos(q4)*(2.5 + 2.5*sin(q5))) +
			cos(q3)*(-2.5*cos(q4)*cos(q5) +
			sin(q4)*(2.5 + 2.5*sin(q5)))) +
			1.0*(1.0*cos(q3 + q4 + q5) -
			1.0*sin(q3 + q4))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			149.625*(dotq5*(cos(
			q3)*(-10.0 - 25.0*cos(q5)*sin(q4) +
			cos(q4)*(7.5*cos(q5) - 25.0*sin(q5)) - 10.0*sin(q5) -
			7.5*sin(q4)*sin(q5)) +
			sin(q3)*(-3.0 +
			cos(q5)*(-20.0 - 25.0*cos(q4) - 7.5*sin(q4)) -
			3.0*sin(q5) - 7.5*cos(q4)*sin(q5) +
			sin(q4)*(-8.0 + 17.0*sin(q5)))) +
			8.0*(1.0*cos(q3) + 1.25*cos(q3 + q4) -
			0.375*cos(q3 + q4 + q5) + 0.375*sin(q3 + q4) -
			0.5*sin(q3 - q5) + 0.5*sin(q3 + q5) +
			1.25*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5)) +
			dotq4*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(-10.0*cos(q3 + q4) +
			3.0*cos(q3 + q4 + q5) - 3.0*sin(q3 + q4) -
			10.0*sin(q3 + q4 + q5)) -
			3.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(1.0*cos(q3 + q4) -
			1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
			3.33333*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
			3.33333*sin(q3 + q4) + 1.0*sin(q3 + q4 + q5))))) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*15.0*sin(
			q3 + q4 +
			q5)*(-7.5*dotq5*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			18.75*dotq4*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(cos(
			q3)*(1.0*cos(q4)*cos(q5) + sin(q4)*(-1.0 - 1.0*sin(q5))) +
			sin(q3)*(-0.8 - 1.0*cos(q5)*sin(q4) +
			cos(q4)*(-1.0 - 1.0*sin(q5)) - 0.8*sin(q5))) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(1.0*cos(2*q3 + q4) + 
			2.5*cos(q3 + q4 + q5) - 2.5*cos(2*q3 + q4 + q5) -
			2.5*sin(q3 + q4) + 2.5*sin(2*q3 + q4) +
			1.0*sin(2*q3 + q4 + q5))) +
			15.0*cos(q3 + q4 +
			q5)*(dotq5*(6*cos(q4 + q5)*sin(q3) + 6*cos(q3)*sin(q4 + q5) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(cos(q3)*cos(q4 + q5) -
			1.0*sin(q3)*sin(q4 + q5)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) + (
			3.0*cos(
			q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*- (
			6.0*sin(q4 + q5)*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			dotq4*(6*sin(q3)*(cos(q4 + q5) - sin(q4)) +
			6*cos(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(1.0*cos(q3 + q4 + q5) - 1.0*sin(q3 + q4)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*3.0*(1.0*cos(q4 + q5) -
			1.0*sin(q4))*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))) - (
			3.0*(2.0*cos(q4) +
			2.0*sin(q4 + q5))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*7.5*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(-1.6*cos(q3) - 2.0*cos(q3 + q4) +
			1.0*cos(2*q3 + q4) + 0.2*cos(q3 + q4 + q5) +
			0.4*cos(2*q3 + q4 + q5) - 0.2*sin(q3 + q4) -
			0.4*sin(2*q3 + q4) + 0.8*sin(q3 - q5) - 0.8*sin(q3 + q5) -
			2.0*sin(q3 + q4 + q5) + 1.0*sin(2*q3 + q4 + q5))) +
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*7.5*(-1.0 + 1.0*cos(q3) - 0.4*sin(q3))*sin(
			q3 + q4 + q5)*(15.0*dotq5*cos(q3 + q4 + q5) +
			dotq3*(15.0*cos(q3 + q4 + q5) - 75.0*sin(q3 + q4)) +
			dotq4*(15.0*cos(q3 + q4 + q5) - 75.0*sin(q3 + q4)) + (
			1215.0*(1.0*dotq3*cos(
			q3) + (-0.0123457*dotq4 - 0.0123457*dotq5)*cos(q4 + q5) +
			2.5*dotq3*sin(q3) +
			0.0617284*dotq4*sin(q4))*(1.0*cos(q3 + q4) +
			1.0*sin(q3 + q4 + q5)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(-847.5 +
			3037.5*cos(q3) + 75.0*cos(q4) - 1215.0*sin(q3) +
			15.0*sin(q4 + q5))*(dotq5*(0.5*cos(q3 - q5) -
			0.5*cos(q3 + q5) - 2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			dotq4*(sin(
			q3)*(2.0 + 2.5*cos(q5)*sin(q4) + 2.0*sin(q5) +
			cos(q4)*(2.5 + 2.5*sin(q5))) +
			cos(q3)*(-2.5*cos(q4)*cos(q5) + 
			sin(q4)*(2.5 + 2.5*sin(q5)))) +
			1.0*(1.0*cos(q3 + q4 + q5) -
			1.0*sin(q3 + q4))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 +
			q5), 2)*405.0*(-7.5*dotq5*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			18.75*dotq4*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(cos(
			q3)*(1.0*cos(q4)*cos(q5) + sin(q4)*(-1.0 - 1.0*sin(q5))) +
			sin(q3)*(-0.8 - 1.0*cos(q5)*sin(q4) +
			cos(q4)*(-1.0 - 1.0*sin(q5)) - 0.8*sin(q5))) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(1.0*cos(2*q3 + q4) +
			2.5*cos(q3 + q4 + q5) - 2.5*cos(2*q3 + q4 + q5) -
			2.5*sin(q3 + q4) + 2.5*sin(2*q3 + q4) +
			1.0*sin(2*q3 + q4 + q5)))) + (-6*cos(q3 + q4 + q5) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*sin(q3 + q4 + q5))/(
			2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)))*(dotq3*(-75.0*cos(q3 + q4) -
			15.0*sin(q3 + q4 + q5)) +
			dotq4*(-75.0*cos(q3 + q4) - 15.0*sin(q3 + q4 + q5)) -
			15.0*dotq5*sin(q3 + q4 + q5) + (
			75.0*(40.5*dotq3*cos(q3) + 1.0*dotq4*cos(q4) -
			16.2*dotq3*sin(q3) + 0.2*dotq4*sin(q4 + q5) +
			0.2*dotq5*sin(q4 + q5))*(1.0*cos(q3 + q4) +
			1.0*sin(q3 + q4 + q5)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(369.0 -
			1215.0*cos(q3) + 15.0*cos(q4 + q5) - 3037.5*sin(q3) -
			75.0*sin(q4))*(dotq5*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			dotq4*(sin(
			q3)*(2.0 + 2.5*cos(q5)*sin(q4) + 2.0*sin(q5) +
			cos(q4)*(2.5 + 2.5*sin(q5))) +
			cos(q3)*(-2.5*cos(q4)*cos(q5) +
			sin(q4)*(2.5 + 2.5*sin(q5)))) +
			1.0*(1.0*cos(q3 + q4 + q5) -
			1.0*sin(q3 + q4))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			405.0*(dotq5*(6*cos(q4 + q5)*sin(q3) + 6*cos(q3)*sin(q4 + q5) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(cos(q3)*cos(q4 + q5) -
			1.0*sin(q3)*sin(q4 + q5)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) + (
			3.0*cos(q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*- (
			6.0*sin(q4 + q5)*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			dotq4*(6*sin(q3)*(cos(q4 + q5) - sin(q4)) +
			6*cos(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(1.0*cos(q3 + q4 + q5) - 1.0*sin(q3 + q4)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*3.0*(1.0*cos(q4 + q5) -
			1.0*sin(q4))*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))) - (
			3.0*(2.0*cos(q4) +
			2.0*sin(q4 + q5))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*7.5*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(-1.6*cos(q3) - 2.0*cos(q3 + q4) +
			1.0*cos(2*q3 + q4) + 0.2*cos(q3 + q4 + q5) +
			0.4*cos(2*q3 + q4 + q5) - 0.2*sin(q3 + q4) -
			0.4*sin(2*q3 + q4) + 0.8*sin(q3 - q5) -
			0.8*sin(q3 + q5) - 2.0*sin(q3 + q4 + q5) +
			1.0*sin(2*q3 + q4 + q5)))) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*1.0*sin(q3 + q4 + q5)*(37.5*dotq1*cos(q3 + q4) -
			7.5*dotq2*cos(q3 + q4 + q5) + 37.5*dotq2*sin(q3 + q4) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*4556.25*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(1.0*dotq3*cos(
			q3) + (-0.0246914*dotq4 - 0.0246914*dotq5)*cos(q4 + q5) +
			0.0123457*dotq4*cos(q3 + q4 + q5) +
			0.0123457*dotq5*cos(q3 + q4 + q5) + 2.5*dotq3*sin(q3) +
			0.123457*dotq4*sin(q4) -
			0.0617284*dotq4*sin(q3 + q4))*(sin(
			q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))) +
			56.25*dotq3*(1.0*cos(q3 + q4) - 0.8*cos(q3 - q5) +
			0.8*cos(q3 + q5) + 1.0*cos(q3 + q4 + q5) - 4.0*sin(q3) -
			5.0*sin(q3 + q4) + 0.2*sin(q3 + q4 + q5)) +
			7.5*dotq1*sin(q3 + q4 + q5) +
			dotq5*(90.0*cos(q3 - q5) + 90.0*cos(q3 + q5) +
			112.5*cos(q4 + q5) + 112.5*cos(q3 + q4 + q5) +
			45.0*sin(q4 + q5) + 22.5*sin(q3 + q4 + q5)) +
			dotq4*(225.0*cos(q4) + 112.5*cos(q3 + q4) + 112.5*cos(q4 + q5) +
			112.5*cos(q3 + q4 + q5) - 562.5*sin(q4) -
			562.5*sin(q3 + q4) + 45.0*sin(q4 + q5) +
			22.5*sin(q3 + q4 + q5)) +
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*56.25*(1.0*cos(q3 + q4) +
			1.0*sin(q3 + q4 + q5))*((-27.0*dotq1 - 10.8*dotq2 +
			4.0*dotq3)*cos(q3) + (4.0*dotq3 - 8.0*dotq4)*cos(q3 - q4) -
			4.0*dotq4*cos(q4) + 1.0*dotq4*cos(q3 + q4) -
			0.8*dotq4*cos(q3 - q5) + 2.0*dotq3*cos(q3 - q4 - q5) -
			4.0*dotq4*cos(q3 - q4 - q5) - 4.0*dotq5*cos(q3 - q4 - q5) -
			3.2*dotq5*cos(q5) + 0.8*dotq4*cos(q3 + q5) +
			0.8*dotq5*cos(q3 + q5) - 4.0*dotq4*cos(q4 + q5) -
			4.0*dotq5*cos(q4 + q5) + 1.0*dotq4*cos(q3 + q4 + q5) +
			1.0*dotq5*cos(q3 + q4 + q5) + 10.8*dotq1*sin(q3) -
			27.0*dotq2*sin(q3) - 132.68*dotq3*sin(q3) -
			4.0*dotq4*sin(q3) - 0.4*dotq5*sin(q3) +
			10.0*dotq3*sin(q3 - q4) - 20.0*dotq4*sin(q3 - q4) +
			20.0*dotq4*sin(q4) - 5.0*dotq4*sin(q3 + q4) -
			0.8*dotq3*sin(q3 - q4 - q5) + 1.6*dotq4*sin(q3 - q4 - q5) + 
			1.6*dotq5*sin(q3 - q4 - q5) - 0.8*dotq4*sin(q4 + q5) -
			0.8*dotq5*sin(q4 + q5) + 0.2*dotq4*sin(q3 + q4 + q5) +
			0.2*dotq5*sin(q3 + q4 + q5)) + (-6*cos(
			q3)*(cos(q4 + q5) - sin(q4)) +
			6*sin(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)))*(-1518.75*dotq3*cos(q3) -
			75.0*dotq4*cos(q4) + 37.5*dotq4*cos(q3 + q4) +
			607.5*dotq3*sin(q3) - 15.0*dotq4*sin(q4 + q5) -
			15.0*dotq5*sin(q4 + q5) + 7.5*dotq4*sin(q3 + q4 + q5) +
			7.5*dotq5*sin(q3 + q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(38416.9 -
			14926.5*cos(q3) + 1125.0*cos(q3 - q4) + 1125.0*cos(q4) -
			90.0*cos(q3 - q4 - q5) - 45.0*cos(q4 + q5) - 450.0*sin(q3) -

			450.0*sin(q3 - q4) + 225.0*sin(q4) - 225.0*sin(q3 - q4 - q5) +
			180.0*sin(q5) +
			225.0*sin(
			q4 + q5))*(dotq5*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			dotq4*(sin(
			q3)*(2.0 + 2.5*cos(q5)*sin(q4) + 2.0*sin(q5) +
			cos(q4)*(2.5 + 2.5*sin(q5))) +
			cos(q3)*(-2.5*cos(q4)*cos(q5) +
			sin(q4)*(2.5 + 2.5*sin(q5)))) +
			1.0*(1.0*cos(q3 + q4 + q5) -
			1.0*sin(q3 + q4))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 +
			q5), 2)*149.625*(dotq5*(cos(
			q3)*(-10.0 - 25.0*cos(q5)*sin(q4) +
			cos(q4)*(7.5*cos(q5) - 25.0*sin(q5)) - 10.0*sin(q5) -
			7.5*sin(q4)*sin(q5)) +
			sin(q3)*(-3.0 +
			cos(q5)*(-20.0 - 25.0*cos(q4) - 7.5*sin(q4)) -
			3.0*sin(q5) - 7.5*cos(q4)*sin(q5) +
			sin(q4)*(-8.0 + 17.0*sin(q5)))) +
			8.0*(1.0*cos(q3) + 1.25*cos(q3 + q4) -
			0.375*cos(q3 + q4 + q5) + 0.375*sin(q3 + q4) -
			0.5*sin(q3 - q5) + 0.5*sin(q3 + q5) +
			1.25*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5)) +
			dotq4*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(-10.0*cos(q3 + q4) +
			3.0*cos(q3 + q4 + q5) - 3.0*sin(q3 + q4) -
			10.0*sin(q3 + q4 + q5)) -
			3.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(1.0*cos(q3 + q4) -
			1.33333*cos(q3 - q5) + 1.33333*cos(q3 + q5) +
			3.33333*cos(q3 + q4 + q5) - 2.66667*sin(q3) -
			3.33333*sin(q3 + q4) + 1.0*sin(q3 + q4 + q5)))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*(-847.5 + 3037.5*cos(q3) + 75.0*cos(q4) -
			1215.0*sin(q3) +
			15.0*sin(q4 + q5))*(-7.5*dotq5*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(0.5*cos(q3 - q5) - 0.5*cos(q3 + q5) -
			2.5*cos(q3 + q4 + q5) + 1.0*sin(q3)) +
			18.75*dotq4*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(cos(
			q3)*(1.0*cos(q4)*cos(q5) + sin(q4)*(-1.0 - 1.0*sin(q5))) +
			sin(q3)*(-0.8 - 1.0*cos(q5)*sin(q4) +
			cos(q4)*(-1.0 - 1.0*sin(q5)) - 0.8*sin(q5))) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(1.0*cos(2*q3 + q4) +
			2.5*cos(q3 + q4 + q5) - 2.5*cos(2*q3 + q4 + q5) -
			2.5*sin(q3 + q4) + 2.5*sin(2*q3 + q4) +
			1.0*sin(2*q3 + q4 + q5))) + (369.0 - 1215.0*cos(q3) +
			15.0*cos(q4 + q5) - 3037.5*sin(q3) -
			75.0*sin(q4))*(dotq5*(6*cos(q4 + q5)*sin(q3) +
			6*cos(q3)*sin(q4 + q5) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(cos(q3)*cos(q4 + q5) -
			1.0*sin(q3)*sin(q4 + q5)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) + (
			3.0*cos(q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*- (
			6.0*sin(q4 + q5)*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			dotq4*(6*sin(q3)*(cos(q4 + q5) - sin(q4)) +
			6*cos(q3)*(cos(q4) + sin(q4 + q5)) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*(1.0*cos(q3 + q4 + q5) - 1.0*sin(q3 + q4)))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*3.0*(1.0*cos(q4 + q5) -
			1.0*sin(q4))*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) +
			2.0*sin(q4))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))) - (
			3.0*(2.0*cos(q4) +
			2.0*sin(q4 + q5))*(sin(q3)*(cos(q4 + q5) - 1.0*sin(q4)) +
			cos(q3)*(cos(q4) + sin(q4 + q5))))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*7.5*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(-1.6*cos(q3) - 2.0*cos(q3 + q4) +
			1.0*cos(2*q3 + q4) + 0.2*cos(q3 + q4 + q5) +
			0.4*cos(2*q3 + q4 + q5) - 0.2*sin(q3 + q4) -
			0.4*sin(2*q3 + q4) + 0.8*sin(q3 - q5) -
			0.8*sin(q3 + q5) - 2.0*sin(q3 + q4 + q5) +
			1.0*sin(2*q3 + q4 + q5))));
	C2[1][0] = value;

	//C2(1,1)

	value = -45*dotq4*cos(q5) -
			7.5*dotq2*cos(q3 + q4 + q5) +
			dotq3*(-45.0*cos(q3 + q5) - 56.25*cos(q4 + q5) -
			56.25*cos(q3 + q4 + q5) - 22.5*sin(q4 + q5) -
			11.25*sin(q3 + q4 + q5)) + 7.5*dotq1*sin(q3 + q4 + q5) - (
			56.25*(1.0*dotq3*cos(
			q4 + q5) + (-2.0*dotq3 - 1.0*dotq4 - 1.0*dotq5)*cos(
			q3 + q4 + q5))*(-1.0 + 1.0*cos(q3) - 0.4*sin(q3))*sin(
			q3 + q4 + q5))/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))*45.0*sin(
			q3 + q4 + q5)*(-1.0*dotq4*cos(q3 - q5) -
			2.5*dotq3*cos(q3 - q4 - q5) - 2.0*dotq3*cos(q5) +
			2.0*dotq3*cos(q3 + q5) - 1.0*dotq4*cos(q3 + q5) +
			1.0*dotq5*cos(q3 + q5) - 0.166667*dotq2*cos(q4 + q5) -
			2.5*dotq3*cos(q4 + q5) + 1.25*dotq4*cos(q4 + q5) +
			1.25*dotq5*cos(q4 + q5) + 2.5*dotq3*cos(q3 + q4 + q5) +
			1.25*dotq4*cos(q3 + q4 + q5) + 1.25*dotq5*cos(q3 + q4 + q5) -
			1.0*dotq3*sin(q3) + 1.0*dotq3*sin(q3 - q4 - q5) +
			0.166667*dotq1*sin(q4 + q5) - 0.5*dotq3*sin(q4 + q5) +
			0.5*dotq4*sin(q4 + q5) + 0.5*dotq5*sin(q4 + q5) +
			0.5*dotq3*sin(q3 + q4 + q5) + 0.25*dotq4*sin(q3 + q4 + q5) +
			0.25*dotq5*sin(q3 + q4 + q5)) + (7.5*dotq3*sin(
			q4 + q5) + (-15.0*dotq3 - 7.5*dotq4 - 7.5*dotq5)*sin(
			q3 + q4 + q5))*(-6*cos(q3 + q4 + q5) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*sin(q3 + q4 + q5))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))) -
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*45.0*(-4.14236 - 1.0*cos(q3) +
			1.0*cos(q4 + q5) + 0.5*cos(q3 + q4 + q5) - 2.0*sin(q3 + q5) -
			2.5*sin(q4 + q5) -
			2.5*sin(q3 + q4 + q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
			q3 + q5) + (-2.5*dotq4 - 2.5*dotq5 +
			1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
			0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
			1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
			0.5*dotq5*sin(2*(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 3)*(4.0*cos(q3 + q5) +
			10.0*cos(q3 + q4 + q5) - 4.0*sin(q3) +
			3.0*sin(q3 + q4 +
			q5))*(149.625*(dotq5*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(3.0*cos(q3 + q4 + q5) -
			4.0*sin(q3 + q5) - 10.0*sin(q3 + q4 + q5)) -
			4.0*cos(q4 + q5)*(1.0*cos(q3 + q5) +
			2.5*cos(q3 + q4 + q5) - 1.0*sin(q3) +
			0.75*sin(q3 + q4 + q5))) +
			dotq4*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(3.0*cos(q3 + q4 + q5) -
			10.0*sin(q3 + q4 + q5)) -
			4.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(1.0*cos(q3 + q5) +
			2.5*cos(q3 + q4 + q5) - 1.0*sin(q3) +
			0.75*sin(q3 + q4 + q5))) +
			4.0*(1.0*cos(q3) - 0.75*cos(q3 + q4 + q5) + 1.0*sin(q3 + q5) +
			2.5*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			149.625*((-1.0*dotq4 - 0.5*dotq5)*cos(
			q3 + q5) + (-2.5*dotq4 - 2.5*dotq5 +
			1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
			0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
			1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
			0.5*dotq5*sin(2*(q3 + q4 + q5)))) + (-6*cos(q3 + q4 + q5) - (
			3.0*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*sin(q3 + q4 + q5))/(
			2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)))*(-15.0*dotq3*sin(q3 + q4 + q5) -
			15.0*dotq4*sin(q3 + q4 + q5) - 15.0*dotq5*sin(q3 + q4 + q5) + (
			75.0*(40.5*dotq3*cos(q3) + 1.0*dotq4*cos(q4) -
			16.2*dotq3*sin(q3) + 0.2*dotq4*sin(q4 + q5) +
			0.2*dotq5*sin(q4 + q5))*sin(q3 + q4 + q5))/(
			2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(369.0 -
			1215.0*cos(q3) + 15.0*cos(q4 + q5) - 3037.5*sin(q3) -
			75.0*sin(q4))*((-1.0*dotq4 - 0.5*dotq5)*cos(
			q3 + q5) + (-2.5*dotq4 - 2.5*dotq5 +
			1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
			0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
			1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
			0.5*dotq5*sin(2*(q3 + q4 + q5))) -
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 +
			q5), 2)*1215.0*(1.0*dotq5*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) +
			2.0*cos(q4 + q5)*(-0.25 - 0.5*cos(q3) + 1.0*cos(q4 + q5) -
			1.25*sin(q3) - 1.0*sin(q4))*sin(q3 + q4 + q5) +
			2.0*sin(q4 + q5)*(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))*sin(
			q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			1.0*dotq4*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) + 2.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(-0.25 - 0.5*cos(q3) +
			1.0*cos(q4 + q5) - 1.25*sin(q3) - 1.0*sin(q4))*sin(
			q3 + q4 + q5) +
			2.0*(1.0*cos(q4) + 1.0*sin(q4 + q5))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*sin(q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			2.5*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(0.8*cos(q3) - 0.2*cos(q3 + q4 + q5) -
			0.4*cos(2*q3 + q4 + q5) + 0.8*sin(q3 + q5) +
			2.0*sin(q3 + q4 + q5) - 1.0*sin(2*q3 + q4 + q5)))) -
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*45.0*cos(
			q3 + q4 +
			q5)*(1.0*dotq5*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) +
			2.0*cos(q4 + q5)*(-0.25 - 0.5*cos(q3) + 1.0*cos(q4 + q5) -
			1.25*sin(q3) - 1.0*sin(q4))*sin(q3 + q4 + q5) +
			2.0*sin(q4 + q5)*(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))*sin(
			q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			1.0*dotq4*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) +
			2.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(-0.25 - 0.5*cos(q3) +
			1.0*cos(q4 + q5) - 1.25*sin(q3) - 1.0*sin(q4))*sin(
			q3 + q4 + q5) +
			2.0*(1.0*cos(q4) + 1.0*sin(q4 + q5))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*sin(q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			2.5*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(0.8*cos(q3) - 0.2*cos(q3 + q4 + q5) -
			0.4*cos(2*q3 + q4 + q5) + 0.8*sin(q3 + q5) +
			2.0*sin(q3 + q4 + q5) - 1.0*sin(2*q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*15.0*sin(
			q3 + q4 +
			q5)*(-1.5*dotq5*(1.0*cos(q3 + q5) + 5.0*cos(q3 + q4 + q5) +
			1.0*cos(q3 + 2*q4 + q5) - 2.0*sin(q3))*(2.5 - 2.5*cos(q3) +
			1.0*sin(q3)) -
			3.0*dotq4*(1.0*cos(q3 + q5) + 2.5*cos(q3 + q4 + q5) -
			1.0*sin(q3))*(2.5 - 2.5*cos(q3) + 1.0*sin(q3)) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(2.5*cos(q3 + q4 + q5) -
			2.5*cos(2*q3 + q4 + q5) + 1.0*sin(2*q3 + q4 + q5))) +
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*7.5*(-1.0 + 1.0*cos(q3) - 0.4*sin(q3))*sin(
			q3 + q4 + q5)*(15.0*dotq3*cos(q3 + q4 + q5) +
			15.0*dotq4*cos(q3 + q4 + q5) + 15.0*dotq5*cos(q3 + q4 + q5) + (
			1215.0*(1.0*dotq3*cos(
			q3) + (-0.0123457*dotq4 - 0.0123457*dotq5)*cos(q4 + q5) +
			2.5*dotq3*sin(q3) + 0.0617284*dotq4*sin(q4))*sin(
			q3 + q4 + q5))/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*(-847.5 + 3037.5*cos(q3) + 75.0*cos(q4) -
			1215.0*sin(q3) +
			15.0*sin(q4 + q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
			q3 + q5) + (-2.5*dotq4 - 2.5*dotq5 +
			1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
			0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
			1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
			0.5*dotq5*sin(2*(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5), 2)*405.0*(-1.5*dotq5*(1.0*cos(q3 + q5) +
			5.0*cos(q3 + q4 + q5) + 1.0*cos(q3 + 2*q4 + q5) -
			2.0*sin(q3))*(2.5 - 2.5*cos(q3) + 1.0*sin(q3)) -
			3.0*dotq4*(1.0*cos(q3 + q5) + 2.5*cos(q3 + q4 + q5) -
			1.0*sin(q3))*(2.5 - 2.5*cos(q3) + 1.0*sin(q3)) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(2.5*cos(q3 + q4 + q5) -
			2.5*cos(2*q3 + q4 + q5) + 1.0*sin(2*q3 + q4 + q5)))) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*1.0*sin(q3 + q4 + q5)*(-7.5*dotq2*cos(q3 + q4 + q5) +
			7.5*dotq1*sin(q3 + q4 + q5) -
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*4556.25*(-1.0 + 1.0*cos(q3) -
			0.4*sin(q3))*(1.0*dotq3*cos(
			q3) + (-0.0246914*dotq4 - 0.0246914*dotq5)*cos(q4 + q5) +
			0.0123457*dotq4*cos(q3 + q4 + q5) +
			0.0123457*dotq5*cos(q3 + q4 + q5) + 2.5*dotq3*sin(q3) +
			0.123457*dotq4*sin(q4) - 0.0617284*dotq4*sin(q3 + q4))*sin(
			q3 + q4 + q5) +
			dotq3*(45.0*cos(q3 + q5) + 56.25*cos(q3 + q4 + q5) -
			22.5*sin(q3) + 11.25*sin(q3 + q4 + q5)) +
			dotq4*(112.5*cos(q4 + q5) + 112.5*cos(q3 + q4 + q5) +
			45.0*sin(q4 + q5) + 22.5*sin(q3 + q4 + q5)) +
			dotq5*(90.0*cos(q3 + q5) + 112.5*cos(q4 + q5) +
			112.5*cos(q3 + q4 + q5) + 45.0*sin(q4 + q5) +
			22.5*sin(q3 + q4 + q5)) +
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*11.25*sin(
			q3 + q4 + 
			q5)*((-135.0*dotq1 - 54.0*dotq2 + 20.0*dotq3)*cos(q3) + (20.0*dotq3 - 40.0*dotq4)*cos(q3 - q4) -
			20.0*dotq4*cos(q4) + 5.0*dotq4*cos(q3 + q4) -
			4.0*dotq4*cos(q3 - q5) + 10.0*dotq3*cos(q3 - q4 - q5) -
			20.0*dotq4*cos(q3 - q4 - q5) - 20.0*dotq5*cos(q3 - q4 - q5) -
			16.0*dotq5*cos(q5) + 4.0*dotq4*cos(q3 + q5) +
			4.0*dotq5*cos(q3 + q5) - 20.0*dotq4*cos(q4 + q5) -
			20.0*dotq5*cos(q4 + q5) + 5.0*dotq4*cos(q3 + q4 + q5) +
			5.0*dotq5*cos(q3 + q4 + q5) + 54.0*dotq1*sin(q3) -
			135.0*dotq2*sin(q3) - 663.4*dotq3*sin(q3) -
			20.0*dotq4*sin(q3) - 2.0*dotq5*sin(q3) +
			50.0*dotq3*sin(q3 - q4) - 100.0*dotq4*sin(q3 - q4) +
			100.0*dotq4*sin(q4) - 25.0*dotq4*sin(q3 + q4) -
			4.0*dotq3*sin(q3 - q4 - q5) + 8.0*dotq4*sin(q3 - q4 - q5) +
			8.0*dotq5*sin(q3 - q4 - q5) - 4.0*dotq4*sin(q4 + q5) -
			4.0*dotq5*sin(q4 + q5) + 1.0*dotq4*sin(q3 + q4 + q5) +
			1.0*dotq5*sin(q3 + q4 + q5)) +
			1/(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))
			*4556.25*(1.0*dotq3*cos(q3) + 0.0493827*dotq4*cos(q4) -
			0.0246914*dotq4*cos(q3 + q4) - 0.4*dotq3*sin(q3) +
			0.00987654*dotq4*sin(q4 + q5) +
			0.00987654*dotq5*sin(q4 + q5) -
			0.00493827*dotq4*sin(q3 + q4 + q5) -
			0.00493827*dotq5*sin(q3 + q4 + q5))*(cos(

			q3 + q4 + q5)*(5.0 + 2.0*cos(q4) + 2.0*sin(q4 + q5)) + (0.5 +
			1.0*cos(q3) - 2.0*cos(q4 + q5) + 2.5*sin(q3) +
			2.0*sin(q4))*sin(q3 + q4 + q5)) +
			1/pow(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 +
			q5), 2)*149.625*(dotq5*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(3.0*cos(q3 + q4 + q5) -
			4.0*sin(q3 + q5) - 10.0*sin(q3 + q4 + q5)) -
			4.0*cos(q4 + q5)*(1.0*cos(q3 + q5) +
			2.5*cos(q3 + q4 + q5) - 1.0*sin(q3) +
			0.75*sin(q3 + q4 + q5))) +
			dotq4*((2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*(3.0*cos(q3 + q4 + q5) -
			10.0*sin(q3 + q4 + q5)) -
			4.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(1.0*cos(q3 + q5) +
			2.5*cos(q3 + q4 + q5) - 1.0*sin(q3) +
			0.75*sin(q3 + q4 + q5))) +
			4.0*(1.0*cos(q3) - 0.75*cos(q3 + q4 + q5) + 1.0*sin(q3 + q5) +
			2.5*sin(q3 + q4 + q5))*(1.0*dotq4*cos(
			q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(38416.9 -
			14926.5*cos(q3) + 1125.0*cos(q3 - q4) + 1125.0*cos(q4) -
			90.0*cos(q3 - q4 - q5) - 45.0*cos(q4 + q5) - 450.0*sin(q3) -
			450.0*sin(q3 - q4) + 225.0*sin(q4) - 225.0*sin(q3 - q4 - q5) +
			180.0*sin(q5) +
			225.0*sin(q4 + q5))*((-1.0*dotq4 - 0.5*dotq5)*cos(
			q3 + q5) + (-2.5*dotq4 - 2.5*dotq5 +
			1.0*dotq4*cos(q3 + q4))*cos(q3 + q4 + q5) -
			0.5*dotq5*cos(q3 + 2*q4 + q5) + 1.0*dotq4*sin(q3) +
			1.0*dotq5*sin(q3) + 0.5*dotq4*sin(2*(q3 + q4 + q5)) +
			0.5*dotq5*sin(2*(q3 + q4 + q5))) -
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*3.0*(369.0 -
			1215.0*cos(q3) + 15.0*cos(q4 + q5) - 3037.5*sin(q3) -
			75.0*sin(q4))*(1.0*dotq5*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) + 
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) +
			2.0*cos(q4 + q5)*(-0.25 - 0.5*cos(q3) + 1.0*cos(q4 + q5) -
			1.25*sin(q3) - 1.0*sin(q4))*sin(q3 + q4 + q5) +
			2.0*sin(q4 + q5)*(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5))*sin(
			q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			1.0*dotq4*(1.0*cos(
			q3 + q4 + q5)*(0.5 + 1.0*cos(q3) - 2.0*cos(q4 + q5) +
			2.5*sin(q3) + 2.0*sin(q4))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5)) +
			2.0*(1.0*cos(q4 + q5) - 1.0*sin(q4))*(-0.25 - 0.5*cos(q3) +
			1.0*cos(q4 + q5) - 1.25*sin(q3) - 1.0*sin(q4))*sin(
			q3 + q4 + q5) +
			2.0*(1.0*cos(q4) + 1.0*sin(q4 + q5))*(2.5 + 1.0*cos(q4) +
			1.0*sin(q4 + q5))*sin(q3 + q4 + q5) -
			2.0*pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*sin(
			q3 + q4 + q5)) +
			2.5*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(0.8*cos(q3) - 0.2*cos(q3 + q4 + q5) -
			0.4*cos(2*q3 + q4 + q5) + 0.8*sin(q3 + q5) +
			2.0*sin(q3 + q4 + q5) - 1.0*sin(2*q3 + q4 + q5))) +
			1/pow(2.5 + 1.0*cos(q4) + 1.0*sin(q4 + q5), 2)*(-847.5 +
			3037.5*cos(q3) + 75.0*cos(q4) - 1215.0*sin(q3) +
			15.0*sin(q4 + q5))*(-1.5*dotq5*(1.0*cos(q3 + q5) +
			5.0*cos(q3 + q4 + q5) + 1.0*cos(q3 + 2*q4 + q5) -
			2.0*sin(q3))*(2.5 - 2.5*cos(q3) + 1.0*sin(q3)) -
			3.0*dotq4*(1.0*cos(q3 + q5) + 2.5*cos(q3 + q4 + q5) -
			1.0*sin(q3))*(2.5 - 2.5*cos(q3) + 1.0*sin(q3)) +
			3.0*(1.0*dotq4*cos(q3 + q4) + (1.0*dotq4 + 1.0*dotq5)*sin(
			q3 + q4 + q5))*(2.5*cos(q3 + q4 + q5) -
			2.5*cos(2*q3 + q4 + q5) + 1.0*sin(2*q3 + q4 + q5))));
	C2[1][1] = value;
}

void vector_sum(float *a, float *b, float *c, int dim){
  int i;

  for(i=0; i<dim; i++)
    c[i] = a[i] + b[i];
}

void vector_sub(float *a, float *b, float *c, int dim){
  int i;

  for(i=0; i<dim; i++)
    c[i] = a[i] - b[i];
}

void vector_scal(float *a, float b, float *c, int dim){
  int i;

  for(i=0; i<dim; i++)
    c[i] = a[i]*b;
}

void vector_copy(float* src, float* dest, int dim){
	int i;

	for(i = 0; i < dim; ++i)
		dest[i] = src[i];

	return;
}

void vector_set_zero(float *v, int dim){
	int i;

	for(i=0; i<dim; ++i)
		v[i] = 0;
	
	return;
}

void matvec_mul(float *x, float *y, int d1, int d2, float A[d1][d2]){
  int i, j;

  for(i=0; i<d1; i++){
    y[i] = 0;
    for(j=0; j<d2; j++)
      y[i] += A[i][j]*x[j];
  }
}

void matrix_print(int row, int column, float m[row][column]){

    int i, j;

    for (i = 0; i < row; i++){
        for (j = 0; j < column; j++)
           printf("%g ", m[i][j]);
        printf("\n");
    }
}

void vector_print(int column, float v[column]){
	int i;

	for (i = 0; i < column; i++){
		printf("%g ", v[i]);
		printf("\n");
	}
}

void matrix_set_zero(int row, int column, float m[row][column]){
	int i, j;
	
	for(i=0; i<row; i++)
		for(j=0; j<column; j++)
			m[i][j] = 0;
		
}

void matrix_inverse(float A[2][2], float res[2][2]){
	float det, a, b, c, d;
	a = A[0][0];
	b = A[0][1];
	c = A[1][0];
	d = A[1][1];

	det = a*d - b*c;

	if(det == 0){
		printf("Matrice non  invertibile: il determinante è nullo.\n");
		return;
	}

	res[0][0] = d/det;
	res[0][1] = -b/det;
	res[1][0] = -c/det;
	res[1][1] = a/det;
}

/********************************TEST OK *************************************/
/*int main(){
	int i;

	state rob;
	dot_state dot_rob;
	float inv_test[2][2];

	rob.q1 = 1;
	rob.q2 = 2;
	rob.q3 = 1;
	rob.q4 = 0;
	rob.q5 = 2;
	rob.q6 = 0;
	dot_rob.dq1 = 2;
	dot_rob.dq2 = 0;
	dot_rob.dq3 = 2;
	dot_rob.dq4 = 2;
	dot_rob.dq5 = 1;
	dot_rob.dq6 = 0;

	float Tsee[4][4], M1[2][2], C1[2][2], S2[4][2], M2[2][2], C2[2][2];
	float G1[2], G2[2];

	update_S2(S2, rob);
	update_kyn(Tsee, rob);
	update_M1(M1, rob);
	update_C1(C1, rob, dot_rob);
	update_G1(G1, rob);
	update_S2(S2, rob);
	update_G2(G2, rob);
	update_M2(M2, rob);
	update_C2(C2, rob, dot_rob);

	matrix_print(4, 4, Tsee);
	matrix_inverse(M1, inv_test);
	matrix_print(2, 2, M1);
	matrix_print(2, 2, inv_test);
	/*matrix_print(4, 2, S2);
	matrix_print(2, 2, M1);
	matrix_print(2, 2, C1);
	matrix_print(2, 2, M2);
	matrix_print(2, 2, C2);
	vector_print(2, G1);
	vector_print(2, G2);


	/*for(i = 0; i < 1000000; i++){
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
