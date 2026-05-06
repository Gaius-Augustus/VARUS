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
    DEBUG(5,"Initializing Runs");
    DEBUG(5,"Done Initializing Runs");
}
