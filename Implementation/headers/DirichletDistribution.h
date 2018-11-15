/*
 * DirichletDistribution.h
 *
 *  Created on: 20.08.2016
 *      Author: willy
 */

#ifndef DIRICHLETDISTRIBUTION_H_
#define DIRICHLETDISTRIBUTION_H_
#include <unordered_map>
#include <vector>
#include "Run.h"
#include "TypeConventions.h"
#include "ParameterHandler.h"
#include "debug.h"

class DirichletDistribution{
	/*! \brief Part of the Dirichlet-Mixture.
	 *
	 * Each Dirichletdistribution can be imagined as a die-factory that produces a certain type of dice.
	 * However the actual dice may vary around this distributio.
	 */

public:
	ParameterHandler *param;

	Run *pseudo1;
	Run *pseudo2;

	DirichletDistribution(ParameterHandler *p, const double w);
	DirichletDistribution();
	~DirichletDistribution();

	void calculateWeight(const unsigned int n);
	void reset();

	void calculateCsum();
	void calculateCj();

	void calculateCCountsAndQij();

	void calculateMixtureParameters(const unsigned int num);

	double EalphaSum();

	double calculateLikelihood(Run *r);

	void resetComplete(unsigned int numOfComps);

	UDmap alpha;
	double weight;			//m
	double alphaSum;
	unsigned int cSum;
//	double weightPseudoCount = 0.1;

	unsigned int size;

	std::vector<Run*> runs;
	UUmap aggregateCCounts;
	UDmap Qij;
	std::unordered_map<U32, UUmap > cj;
	// for each transcript j we save a map, with an entry for each run k.
	// each of these entries is the number of reads in the run k on this transcript j.

	void estimateP(Run *r);

	// experimental!
	double calculateDistance(Run *r);

	void calculateSimpleP();
};

#endif /* DIRICHLETDISTRIBUTION_H_ */
