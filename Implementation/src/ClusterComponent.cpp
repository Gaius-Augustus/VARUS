/*
 * ClusterComponent.cpp
 *
 *  Created on: 29.12.2016
 *      Author: willy
 */

#include "../headers/ClusterComponent.h"
#include "../headers/debug.h"
#include "../headers/Operators.h"
#include "math.h"

ClusterComponent::ClusterComponent(ParameterHandler *p) {
	param = p;
	alphaSum = 0;
}

ClusterComponent::~ClusterComponent() {
	// TODO Auto-generated destructor stub
}

void ClusterComponent::setSeed(UDmap &start){
	/*! \brief Sets the starting point for this cluster.
	 *
	 * For each run the distance is calculated to this starting point
	 * in order to assign it to the run that is closest situated.
	 */

	alpha = start;
	alphaSum = 1;
}

void ClusterComponent::updateCenter(){
	/*! \brief In each step of the k_means algorithm the center needs to be updated.
	 *
	 * This is done by averaging over the runs assigned to this component.
	 */


	if(runs.size() > 0){
		UDmap::iterator j;
		for(j = alpha.begin(); j != alpha.end(); j++){
			for(unsigned int k = 0; k < runs.size(); k++){
				alpha[j->first] += runs[k]->q[j->first];
			}
			alpha[j->first] /= runs.size();
		}
	}
}

void ClusterComponent::releaseRuns(){
	/*! \brief The runs clustered into this component are released again.
	 *
	 * In the next step of the clustering the runs can then be assigned to a component
	 * again.
	 */

	runs.clear();
}

void ClusterComponent::calculateAlphaSum(){
	/*! \brief Calculates the 'concentration-parameter'.
	 *
	 * A higher variance in the frequencies of the runs leads to a smaller \lambda.
	 * This makes the estimation for runs with little reads more precisse.
	 * If the runs in the component only have very little variance, it is a good guess
	 *  to assume that the reads in the runs with only little reads follow a similar
	 *  distribution as the reads in the runs with more observations.
	 */


	if(runs.size() > 1){
		double var = calculateVariance();
		alphaSum = (1.0/(var))*param->numOfBlocks;
	} else{
		alphaSum = param->lambda*param->numOfBlocks;
	}
	DEBUG(4,"AlphaSum: " << alphaSum);
	scaleMap(alpha,alphaSum);
}

void ClusterComponent::calculateP(){
	/*! \brief The estimation p is calculated for all runs assigned to this component in the
	 * clustering-step.
	 *
	 * The estimation is the same as in the advanced estimator, only for a single component.
	 */


	double psum = 0.0;
	for(unsigned int k = 0; k < runs.size(); k++) {
		// c^k_{n_k}* + a/(c*+1)T + lambda
		double normalizingConst = runs[k]->observationSum
								+ (double)param->pseudoCount*param->numOfBlocks
								+ alphaSum;

		DEBUG(4,normalizingConst);
		UDmap::iterator j;
		for(j=alpha.begin(); j!=alpha.end(); j++) {
			// p_j = c_j + a/(c*+1) + lambda*p_total
			runs[k]->p[j->first] =
					(runs[k]->observations[j->first]
					+ (double)param->pseudoCount
					+ alpha[j->first])
					/normalizingConst;

		DEBUG(4,runs[k]->observations[j->first]);
		DEBUG(4,param->pseudoCount);
		DEBUG(4,alpha[j->first]);


		DEBUG(4,runs[k]->p[j->first]);
		}
	}
}

double ClusterComponent::calculateVariance(){
	/*! \brief Calculates the variance in the observations of the runs assigned to this component.
	 *
	 * The euclidean norm is used to get a scalar out of the vector.
	 */


	double var = 1.0;
	UDmap v;
	for(unsigned int k = 0; k < runs.size(); k++){
		UDmap::iterator j;
		for(j = alpha.begin(); j != alpha.end(); j++){
			v[j->first] += (pow(runs[k]->q[j->first] - j->second,2))/runs.size();
		}
	}

	var = euklidNorm(v);

	return var;
}


