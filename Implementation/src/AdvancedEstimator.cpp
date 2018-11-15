/*
 * AdvancedEstimator.cpp
 *
 *  Created on: 31.08.2016
 *      Author: willy
 */

#include "../headers/AdvancedEstimator.h"
#include "../headers/debug.h"
#include "../headers/Operators.h"
#include "assert.h"
#include "math.h"
#include "../headers/Warning.h"

using namespace std;


AdvancedEstimator::AdvancedEstimator(ParameterHandler *p) : Estimator(p) {
}

AdvancedEstimator::~AdvancedEstimator() {
}

void AdvancedEstimator::estimateP(std::vector<Run*> &runs, unsigned int iterationNumber) {
/* 		p_i = (c_i + a + lambda*p_total)normalizingConst
 *
 * p_i: 	Estimator for the i-th Run
 * c_i: 	Observations in the i-th Run
 * a: 		base-pseudoCount
 * lambda: 	Coefficient
 * p_total: Estimator for the all observations combined
 * 	p_total_{i} = observations_total_{i} / sum(observations_total);
 */


//	DEBUG(3,"runs.size() " << runs.size());

	UUmap obs_total;	// c_j
	unsigned int observationSum = 0;							// c^*
	// calculate observations_total and sum(observations)
	for(unsigned int k = 0; k < runs.size(); k++) {

		observationSum 	+= runs[k]->observationSum;
		obs_total 		+= runs[k]->observations;

//		DEBUG(3,"runs[k]->observationSum " << runs[k]->observationSum );
//		DEBUG(3,"runs[k]->observations " << runs[k]->observations );
	}

//	DEBUG(3,"observationSum: " << observationSum);
//	DEBUG(3,"obs_total: " << obs_total);

	UDmap p_total;
	UUmap::iterator j;

	double pC = 1.0;

	unsigned int totalReadsDownloaded = iterationNumber*param->batchSize;

	for(j=obs_total.begin(); j!=obs_total.end(); j++) {
			// p_j = c_j/c*
//			p_total[j->first] = ((double)obs_total[j->first] + pC)/((double)observationSum + pC*obs_total.size());	// pre 7.12.17
			p_total[j->first] = ((double)obs_total[j->first] + pC)/((double)totalReadsDownloaded + pC*obs_total.size()); // after 7.12.17

//			DEBUG(3,"p_total[j->first] " << p_total[j->first]);
	}

	UDmap pnoObs;
	// Make estimation for runs with 0 reads only once
	{
		double normalizingConst = 	+ (double)param->pseudoCount*param->numOfBlocks
									+ param->lambda*param->numOfBlocks;

		UDmap::iterator j;
		// runs[k]->p.size() will be equal to totalObs.size() afterwards even though
		// runs[k]->getObservationsSize may be smaller
		for(j=p_total.begin(); j!=p_total.end(); j++) {
			// p_j = c_j + a/(c*+1) + lambda*p_total
			pnoObs[j->first] =
					((double)param->pseudoCount
					+ param->lambda*p_total[j->first]*param->numOfBlocks)
					/normalizingConst;
		}
	}

	double psum = 0.0;
	bool faultyCount = 0;
	for(unsigned int k = 0; k < runs.size(); k++) {
		// c^k_{n_k}* + a/(c*+1)T + lambda
//		double normalizingConst = runs[k]->observationSum
//								+ (double)param->pseudoCount*param->numOfBlocks
//								+ param->lambda*param->numOfBlocks;								// pre 7.12.17

		if(runs[k]->timesDownloaded == 0){
			runs[k]->p = pnoObs;
		} else {

			double normalizingConst = runs[k]->timesDownloaded*param->batchSize
									+ (double)param->pseudoCount*param->numOfBlocks
									+ param->lambda*param->numOfBlocks;								// after 7.12.17

	//		DEBUG(1,normalizingConst);
			UDmap::iterator j;
			// runs[k]->p.size() will be equal to totalObs.size() afterwards even though
			// runs[k]->getObservationsSize may be smaller
			for(j=p_total.begin(); j!=p_total.end(); j++) {
				// p_j = c_j + a/(c*+1) + lambda*p_total
				runs[k]->p[j->first] =
						(runs[k]->observations[j->first]
						+ (double)param->pseudoCount
						+ param->lambda*p_total[j->first]*param->numOfBlocks)
						/normalizingConst;
				psum += runs[k]->p[j->first];
			}

			// pNoObs = a/ c^k* + Ta + lambda
	//		runs[k]->setPNoObs(basePseudoCount/
	//										(runs[k]->getObservationSum()
	//										+ runs[k]->getSize()*basePseudoCount
	//										+ lambda*sum(p_total)));
	//		psum += runs[k]->getPNoObs()*(runs[k]->getSize() - runs[k]->getPSize());
	//		DEBUG(1, psum);
	//		assert(fabs(1.0 - psum) <= precision);



	//		if(fabs(1.0 - psum) > precision) {
	////			DEBUG((fabs(1.0 - psum) << " > " << precision);
	//			faultyCount++;
	//		}
			psum = 0.0;
		}

		DEBUG(2,"p_total sums to " << sum(p_total));
	}

//	if(faultyCount != 0) {
//		DEBUG(0,"WARNING: " << faultyCount << " runs probably have wrong estimations! P does not sum to 1!");
//	}

}

//void AdvancedEstimator::initializeP(std::vector<Experiment*> &runs) {
//	for(unsigned int k = 0; k < runs.size(); k++) {
//		runs[k]->setAllPValues(1.0/(double)runs[k]->getSize());
//	}
//}

void AdvancedEstimator::estimatePWithNoObs(std::vector<Run*> &runs, unsigned int iterationNumber){

}
