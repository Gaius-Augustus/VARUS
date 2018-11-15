/*
 * DirichletMixture.h
 *
 *  Created on: 20.08.2016
 *      Author: willy
 */

/*! \brief The dirichlet-mixture is used, when it is assumed that the runs
 *  are not independently distributed, but instead a handfull of distributions represent all runs.
 *
 *  Detailed description starts here.
 */

#ifndef DIRICHLETMIXTURE_H_
#define DIRICHLETMIXTURE_H_
#include "DirichletDistribution.h"
#include "Estimator.h"
#include "Run.h"
#include "ParameterHandler.h"
#include <map>

#include <random>
#include <stdio.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "myRandomEngine.h"

class DirichletMixture : public Estimator {


public:
	ParameterHandler *param;
	myRandomEngine *sampler;

	std::vector<DirichletDistribution> components;

	unsigned int current_trainingsIteration;
	unsigned int iterationNumber;

	std::unordered_map<double,double> lgammaHash;
//	std::unordered_map<double,double> diGammaHash;


//	double calculateLikelihood(DirichletDistribution &component, Run *r);
//	double digamma(const double x);		//http://lib.stat.cmu.edu/apstat/103
//	double trigamma(const double x);	//http://lib.stat.cmu.edu/apstat/103

	double logLikelihood(DirichletDistribution &comp);
	double logLikelihood(DirichletDistribution &comp, const double alphaSum);
	double dL(DirichletDistribution &comp, const double x);
	double ddL(DirichletDistribution &comp, const double x);
	double newtonsMethod(unsigned int compNum, DirichletDistribution &c, const double start, unsigned int maxIterations, const double tol);

	double XiLog(std::vector<double> &logBComponents2, Run *r, const U32 i);
	double myLogBeta(UDmap &v);
	double myLogBeta(UDmap &v, UUmap &v2);


	void calculateMixtureParameters(const unsigned int num);
	void train(std::vector<Run*> &runs);
	void calculateLikelihoods(Run* r, std::vector<double> &likelihoods);
	bool checkComponentsValid();


	void estimatePWithNoObs(std::vector<Run*> &runs, unsigned int iterationNumber);
	void estimateP(std::vector<Run*> &runs, const unsigned int IterationNumber);
	void initializeP(std::vector<Run*> &runs);
	void calculateP(std::vector<Run*> &runs);


	void exportLikelihoodCSV(const unsigned int it, const unsigned int compNum, DirichletDistribution &comp, const double a, const double b, const int n, std::vector<double> &newton);
	void exportMixture(const unsigned int iterationNumber, USmap &translate2Str);
	USmap translate2Str;

	void setTranslate2Str(USmap &translate2Str);


	void exportRunNames();


	void calculateSimpleP();

	DirichletMixture(ParameterHandler *p);
	DirichletMixture();
	virtual ~DirichletMixture();

};

#endif /* DIRICHLETMIXTURE_H_ */
