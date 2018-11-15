/*
 * ClusterComponent.h
 *
 *  Created on: 29.12.2016
 *      Author: willy
 */

#ifndef CLUSTERCOMPONENT_H_
#define CLUSTERCOMPONENT_H_
#include "TypeConventions.h"
#include "Run.h"
#include "ParameterHandler.h"

class ClusterComponent {
	/*! \brief This Class represents a component of the cluster used in the clusterEstimator.
	 *
	 * Each cluster is assigned a part of the runs during clustering.
	 * Afterwards the runs are estimated based on the runs in the same cluster with them.
	 *
	 */


public:
	ParameterHandler *param;

	UDmap alpha;		// cluster-center
	double alphaSum;
	std::vector<Run*> runs;

	void setSeed(UDmap &start);

	void updateCenter();

	void releaseRuns();

	void calculateAlphaSum();

	void calculateP();

	double calculateVariance();

	ClusterComponent(ParameterHandler *p);
	virtual ~ClusterComponent();
};

#endif /* CLUSTERCOMPONENT_H_ */
