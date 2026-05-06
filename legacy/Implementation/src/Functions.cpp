/*
 * Functions.cpp
 *
 *  Created on: 15.10.2016
 *      Author: willy
 */
#include <unordered_map>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include "../headers/Functions.h"
#include "../headers/Warning.h"

using namespace std;

double digamma(const double &x) {
	/*! \brief Implementation of the digamma-function.
	 *
	 *  ALGORITHM AS 103  APPL. STATIST. (1976) VOL.25, NO.3
	 */

	// static needed so it isn't destroyed at end of function
//	static unordered_map<double, double> diGammaHash;
//
//	if(diGammaHash.find(x) != diGammaHash.end()) {
//		return diGammaHash[x];
//	} else {

	//C
	//C     ALGORITHM AS 103  APPL. STATIST. (1976) VOL.25, NO.3
	//C
	//C     Calculates DIGAMMA(X) = D( LOG( GAMMA(X))) / DX
	//C
	//     REAL ZERO, HALF, ONE
		double zero, half, one;
	//C
	//C     Set constants, SN = Nth Stirling coefficient, D1 = DIGAMMA(1.0)
	//C
	//     DATA ZERO/0.0/, HALF/0.5/, ONE/1.0/
	//     DATA S, C, S3, S4, S5, D1 /1.E-05, 8.5, 8.333333333E-02,
	//    *    8.3333333333E-03, 3.96825 3968E-03, -0.57721 56649/
		zero = 0.0;
		half = 0.5;
		one = 1.0;
		double s = 1.e-05;
		double c = 8.5;
		double s3 = 8.333333333e-02;
		double s4 = 8.3333333333e-03;
		double s5 = 3.968253968e-03;
		double d1 = -0.5772156649;


	//C
	//C     Check argument is positive
	//C
	//     DIGAMA = ZERO
	//     Y = X
	//     IFAULT = 1
	//     IF (Y .LE. ZERO) RETURN
	//     IFAULT = 0

		double DIGAMA = zero;
		double y = x;
		if(y <= 0){
			std::cerr << "ERROR::Digamma not defined for x <= 0! Exiting Application." << std::endl;
			exit(0);
		}

	//C
	//C     Use approximation if argument <= S
	//C
	//     IF (Y .LE. S) THEN
	//	DIGAMA = D1 - ONE / Y
	//	RETURN
	//     END IF
		if(y <= s) {
			return d1 -one/x;
		}

	//C
	//C     Reduce to DIGAMA(X + N) where (X + N) >= C
	//C
	//   1 IF (Y .GE. C) GO TO 2
	//     DIGAMA = DIGAMA - ONE/Y
	//     Y = Y + ONE
	//     GO TO 1

		while(y < c) {
			DIGAMA -= one/y;
			y = y + one;
		}

	//C
	//C     Use Stirling's (actually de Moivre's) expansion if argument > C
	//C
	//   2 R = ONE / Y
	//     DIGAMA = DIGAMA + LOG(Y) - HALF*R
	//     R = R * R
	//     DIGAMA = DIGAMA - R*(S3 - R*(S4 - R*S5))
	//     RETURN
	//     END
		double r = one/y;
		DIGAMA = DIGAMA + log(y) -half*r;
		r = r*r;
		DIGAMA = DIGAMA - r*(s3-r*(s4-r*s5));

//		diGammaHash[x] = DIGAMA;

		return DIGAMA;
//	}
}

double trigamma(const double &x) {
	/*! \brief Implementation of the trigamma-function.
	 *
	 *  algorithm as121   Appl. Statist. (1978) vol 27, no. 1
	 */

//double precision function trigam(x, ifault)
//implicit double precision (a-h,o-z)
//c
//c        algorithm as121   Appl. Statist. (1978) vol 27, no. 1
//c
//c        calculates trigamma(x) = d**2(log(gamma(x))) / dx**2
//c
//double precision a, b, one, half, b2, b4, b6,b8, x, y, z, zero
//data a, b, one, half /1.0d-4, 5.0d0, 1.0d0, 0.5d0/
//data zero /0.0d0/



//	double a, b, one, half, b2, b4, b6, b8, y, z, zero;
//	a = 1.0e-04;
//	b = 5.0;
//	one = 1.0;
//	half = 0.5;
//	zero = 0.0;
//
//
////c
////c        b2, b4, b6 and b8 are Bernoulli numbers
////c
////data b2, b4, b6,b8
////*/0.1666666667d0, -0.03333333333d0, 0.02380952381, -0.03333333333/
//

//	b2 = 0.1666666667;
//	b4 = -0.03333333333;
//	b6 = 0.02380952381;
//	b8 = -0.03333333333;


//c
//c        check for positive value of x
//c
//trigam = zero
//ifault = 1
//if (x.le.zero) return
//ifault = 0
//z = x

	double y,z;
	double trigam = zero;
	if(x <= 0){
		std::cerr << "ERROR::Digamma not defined for x <= 0! Exiting Application." << std::endl;
		exit(0);
	}
	z = x;

//c
//c        use small value approximation if x .le. a
//c
//if (z .gt. a) goto 10
//trigam = one / (z * z)
//return

	if(z < a) {
		return one/(z*z);
	}

//c
//c        increase argument to (x+i) .ge. b
//c
//10 if (z .ge. b) goto 20
//trigam = trigam + one / (z * z)
//z = z + one
//goto 10

	while(z < b) {
		trigam = trigam + one/(z*z);
		z = z+one;
	}

//c
//c        apply asymptotic formula if argument .ge. b
//c
//20 y = one / (z * z)
//trigam = trigam + half * y +
//* (one + y * (b2 + y * (b4 + y * (b6 + y * b8)))) / z
//return
//end

	y = one/(z*z);
	trigam = trigam + half*y + (one + y*(b2 + y*(b6 + y*b8))) / z;

	return trigam;
}

double Emean(UUmap &c) {
	/*! \brief Calculates the artihmetic mean of the values in c.
	 *
	 */

	double y = 0;
	UUmap::iterator j;
	for(j=c.begin();j!=c.end();j++) {
		y += (double)j->second/(double)c.size();
	}
	return y;
}

double Evariance(UUmap &c) {
	/*! \brief Calculates the variance in the values of c.
	 *
	 */

	double y = 0;
	double mean = Emean(c);
	UUmap::iterator j;
	for(j=c.begin();j!=c.end();j++) {
		y += pow((double)j->second - (double)mean,2)/((double)c.size()-1);
	}
	return y;
}

double Vj(UUmap &cj, UUmap &c, const double qj) {
	/*! \brief Below formula (12) in Altschull.
	 *
	 */

	return (double)Evariance(cj)/(double)pow(Emean(c),2) - ((double)Evariance(c)/(double)pow(Emean(c),2))*(double)pow(qj,2);
}

double newtonsMethod(double start, double (*f)(double), double (*df)(double)){
	double x_k;
//	if(start == 0) {
//		WARN("Warning start-value == 0. Setting x_k = 1.0");
//		x_k = 1.0;
////		x_k *=-1;
//	} else if(start < 0){
//		WARN("Warning start-value < 0. Setting x_k *= -1.0");
//		x_k *= -1.0;
//	}
//	else{
//		x_k = start;
//	}
//
//	double x_o = x_k;
//
//	vector<double> values;
//	values.push_back(x_k);
//
//	DEBUG(4,"newtonsMethod(" << start << ", " << maxIterations << ", " << tol << ")");
//
//
//
//	for(unsigned int i = 0; i < maxIterations; i++) {
//		x_o = x_k;
//		DEBUG(1,"x[" << i << "] = " << x_k);
////		x_k = x_k - F1(x_k)/dF1(x_k);
//		if(ddL(c, x_k) == 0){
//			DEBUG(1,"ddL(c, x_k) == 0");
//			DEBUG(1,"ddL: " << ddL(c,x_k));
//			break;
//		}
//
//		x_k = x_k - dL(c, x_k)/ddL(c, x_k);
//		values.push_back(x_k);
//
//		if(fabs((double)(x_k - x_o)) < tol) { break;}
//		if(x_k <= 0) {
//			DEBUG(1,x_k << " <= 0; set x_k := " << x_o*0.5);
//			x_k = 0.5*x_o;
//			values.push_back(x_k);
//		}	// <----------------------------------
//	}
//
////	exportLikelihoodCSV(current_trainingsIteration,compNum,c,0.001,2*x_k,1000,values);
//	DEBUG(4,"ddL: " << ddL(c,x_k));
//	return x_k;
	return x_k;
}


