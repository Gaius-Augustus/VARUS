/*
 * Estimator.h
 *
 *  Created on: 25.08.2016
 *      Author: willy
 */

#ifndef ESTIMATOR_H_
#define ESTIMATOR_H_
#include <unordered_map>
#include "Run.h"
#include "ParameterHandler.h"
#include <vector>
#include <memory>
#include "TypeConventions.h"

class Estimator {
public:
    constexpr static double precision = 0.0001;
    ParameterHandler *param;
    Estimator();
    Estimator(ParameterHandler *p);
    virtual ~Estimator();
    
    virtual void estimatePWithNoObs(std::vector<Run*> &runs, unsigned int iterationNumber) = 0;
    virtual void estimateP(std::vector<Run*> &runs, unsigned int iterationNumber) = 0;
    void initializeRuns(std::vector<Run*> &runs, UUmap &transcriptUnits);
};

#endif /* ESTIMATOR_H_ */
