/*
 * ClusterEstimator.h
 *
 *  Created on: 29.12.2016
 *      Author: willy
 */

#ifndef CLUSTERESTIMATOR_H_
#define CLUSTERESTIMATOR_H_
#include "Estimator.h"
#include <unordered_map>
#include <vector>
#include "ClusterComponent.h"
#include "myRandomEngine.h"

class ClusterEstimator : public Estimator {
	/*! \brief Calculates the estimation for the runs by first clustering the runs based on the
	 * similarity of the observations.
	 *
	 * This estimator is inspired by the Dirichlet-Mixture and aims to yield reasonable results
	 * while still having affordable computation-times.
	 */

public:

	std::vector<ClusterComponent> components;
	myRandomEngine *ran;

	ClusterEstimator(ParameterHandler *p);
	virtual ~ClusterEstimator();

	void calculatePnoObs(std::vector<Run*> &runs);

	void calculateRunQs(std::vector<Run*> &runs);

	void initializeClusters(std::vector<Run*> &runs);

	void kMeans(std::vector<Run*> &runs);

	void estimatePWithNoObs(std::vector<Run*> &runs, unsigned int iterationNumber);

	void estimateP(std::vector<Run*> &runs, unsigned int iterationNumber);

	double distance(Run *r, ClusterComponent &c);

	void exportClusters(unsigned int iterationNumber);

//	bool pred(Run* s);
};

#endif /* CLUSTERESTIMATOR_H_ */
