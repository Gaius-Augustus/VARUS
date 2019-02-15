/*
 * AdvancedEstimator.cpp
 *
 *  Created on: 31.08.2016
 *      Authors: willy, Mario Stanke
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
/* 		p_i = (c_i + a + lambda*p_total_{i})*normalizingConst
 *
 * p_i: 	Estimator for the i-th Run
 * c_i: 	Observations in the i-th Run
 * a: 		base-pseudoCount
 * lambda: 	Coefficient
 * p_total: Estimator for the all observations combined
 * 	p_total_{i} = observations_total_{i} / sum(observations_total);
 */


    DEBUG(2, "runs.size() " << runs.size());

    UUmap obs_total;	// c_j
    // calculate observations_total and sum(observations)
    for (unsigned int k = 0; k < runs.size(); k++) {
	obs_total += runs[k]->observations; // insert into map
    }
    // DEBUG(2,"obs_total: " << obs_total);

    UDmap p_total;
    double pC = 0.0;
    
    unsigned int totalAliCount = 0;
    for (UUmap::iterator j = obs_total.begin(); j != obs_total.end(); j++)
	totalAliCount += j->second;
    
    double normalizer = 1.0 / (totalAliCount + pC * obs_total.size());
    for (UUmap::iterator j = obs_total.begin(); j != obs_total.end(); j++) 
	p_total[j->first] = normalizer * (pC + obs_total[j->first]);
    
    DEBUG(1, "total alignment count=" << totalAliCount);

    // Make only once the estimation for runs with 0 reads
    UDmap pnoObs;
    double normalizingConst = (param->pseudoCount + param->lambda) * param->numOfBlocks;	    

    DEBUG(2, "Making estimation for runs not downloaded yet...");
    
    normalizer = 0.0;
    for (UDmap::iterator i = p_total.begin(); i != p_total.end(); i++)
	normalizer += param->pseudoCount  + param->lambda * p_total[i->first] * param->numOfBlocks;
    normalizer = 1.0 / normalizer;
    for (UDmap::iterator i = p_total.begin(); i != p_total.end(); i++){
	pnoObs[i->first] = normalizer * (param->pseudoCount  + param->lambda * p_total[i->first] * param->numOfBlocks);
	// with lamba=1 the pseudocount from the background distribution from all runs is as if there was on average
	// 1 read from each block
	DEBUG(2, "pseudoCount=" << param->pseudoCount << "\tlambda* p_total[i->first] * param->numOfBlocks= " << param->lambda
	      << " * " <<  p_total[i->first] << " * " << param->numOfBlocks << " = " << param->lambda * p_total[i->first] * param->numOfBlocks
	      << "\tpnoObs=" << pnoObs[i->first]);
    }
    DEBUG(2, "p_total sums to " << sum(p_total));

    Run *firstRunWoObs = nullptr; // this run's will represent all runs p, that have not been downloaded yet
    for (unsigned int k = 0; k < runs.size(); k++) {
	runs[k]->pRep = nullptr; // unless warranted otherwise no run has another representative
	if (runs[k]->timesDownloaded == 0){
	    if (firstRunWoObs == nullptr){ 
		runs[k]->p = pnoObs; // first run encountered that has not been downloaded
		firstRunWoObs = runs[k];
	    } else {
		runs[k]->pRep = firstRunWoObs;
	    }
	} else {
	    double normalizingConst = 0.0;
	    UDmap::iterator j;
	    for (j = p_total.begin(); j != p_total.end(); j++)
		normalizingConst += runs[k]->observations[j->first] + param->pseudoCount
		    + param->lambda*p_total[j->first] * param->numOfBlocks;
	    for (j = p_total.begin(); j != p_total.end(); j++) {
		runs[k]->p[j->first] = (runs[k]->observations[j->first]	+ param->pseudoCount
					+ param->lambda*p_total[j->first]*param->numOfBlocks) / normalizingConst;
	    }
	}
    }
}


void AdvancedEstimator::estimatePWithNoObs(std::vector<Run*> &runs, unsigned int iterationNumber){
}
