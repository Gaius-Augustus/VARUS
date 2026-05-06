/*
 * Functions.h
 *
 *  Created on: 15.10.2016
 *      Author: willy
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "TypeConventions.h"
#include "math.h"

double digamma(const double &x);
/*
 * parameters for trigamma()
 */
const static double a = 1.0e-04;
const static double b = 5.0;
const static double one = 1.0;
const static double half = 0.5;
const static double zero = 0.0;


//c
//c        b2, b4, b6 and b8 are Bernoulli numbers
//c
//data b2, b4, b6,b8
//*/0.1666666667d0, -0.03333333333d0, 0.02380952381, -0.03333333333/

const static double b2 = 0.1666666667;
const static double b4 = -0.03333333333;
const static double b6 = 0.02380952381;
const static double b8 = -0.03333333333;

double trigamma(const double &x);

double Emean(UUmap &c);
double Evariance(UUmap &c);

double Vj(UUmap &cj, UUmap &c, const double qj);

double newtonsMethod(double start, double (*f)(double), double (*df)(double));

#endif /* FUNCTIONS_H_ */
