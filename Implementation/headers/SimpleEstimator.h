#ifndef SIMPLEESTIMATOR_H_
#define SIMPLEESTIMATOR_H_
#include "Estimator.h"
#include "Run.h"
#include <unordered_map>

class SimpleEstimator : public Estimator {
public:
//	double pseudoCount;	// value added to each observations for a less strict estimator

	SimpleEstimator(ParameterHandler *p);
	virtual ~SimpleEstimator();

	void estimateP(std::vector<Run*> &runs, unsigned int IterationNumber);
	void estimatePWithNoObs(std::vector<Run*> &runs, unsigned int iterationNumber);
};

#endif /* SIMPLEESTIMATOR_H_ */
