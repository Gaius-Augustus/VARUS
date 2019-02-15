/*
 * AdvancedEstimator.h
 *
 *  Created on: 31.08.2016
 *     Authors: willy, Mario Stanke
 */

#ifndef ADVANCEDESTIMATOR_H_
#define ADVANCEDESTIMATOR_H_
#include "Estimator.h"
#include <unordered_map>
#include <vector>

class AdvancedEstimator : public Estimator {
    /*! \brief The advanced estimator makes use of a pseudocount and a special pseudocount
     * being distributed like a mixture of all runs combined.
     */

 public:
    AdvancedEstimator(ParameterHandler *p);
    virtual ~AdvancedEstimator();
    void estimatePWithNoObs(std::vector<Run*> &runs, unsigned int iterationNumber);
    void estimateP(std::vector<Run*> &runs, unsigned int iterationNumber);
    //	void initializeP(std::vector<Experiment*> &experiments);
};

#endif /* ADVANCEDESTIMATOR_H_ */
