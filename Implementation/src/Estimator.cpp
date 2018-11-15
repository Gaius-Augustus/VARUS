/*
 * Estimator.cpp
 *
 *  Created on: 07.12.2016
 *      Author: willy
 */

#include "../headers/Estimator.h"
#include "../headers/debug.h"

Estimator::Estimator() {
	param = NULL;

}

Estimator::Estimator(ParameterHandler *p) {
	param = p;
}

Estimator::~Estimator() {
	// TODO Auto-generated destructor stub
}

void Estimator::initializeRuns(std::vector<Run*> &runs, UUmap &transcriptUnits){
	DEBUG(1,"Initializing Runs");
	// this is actually not necessary
//	UUmap::iterator j;
//	for(unsigned int i = 0; i < runs.size(); i++){
//		for(j = transcriptUnits.begin(); j != transcriptUnits.end(); j++){
//			runs[i]->observations[j->first] = 0;
//			runs[i]->p[j->first] = 1.0/(double)transcriptUnits.size();
//		}
//	}
	DEBUG(1,"Done Initializing Runs");
}
