#include "../headers/SimpleEstimator.h"
#include "../headers/debug.h"
#include "../headers/Operators.h"

using namespace std;

SimpleEstimator::SimpleEstimator(ParameterHandler *p) : Estimator(p) {
    // TODO Auto-generated constructor stub
}

SimpleEstimator::~SimpleEstimator() {
    // TODO Auto-generated destructor stub
}

void SimpleEstimator::estimateP(std::vector<Run*> &runs, unsigned int IterationNumber) {
    /*	the estimator is calculated as follows:
     * 	p_j := (c_j + a)/(c* + aT),
     * 	a: pseudoCount
     * 	c: observations in this run
     * 	c*: sum of observations in this run, 	experiments[k]->getObservationSum()
     * 	T: number of transcript-units, 			experiments[k]->getSize()
     */

    for(unsigned int k = 0; k < runs.size(); k++) {
	// experiments[k]->getSize() holds the actual number of possible outcomes of an Experiment
	//		double totalCount = (double)pseudoCount*runs[k]->size
	//										+ runs[k]->getObservationSum();
	
	double totalCount = (double)param->pseudoCount*param->numOfBlocks + sum(runs[k]->observations);


//		TODO: which one correct?
//		unordered_map<unsigned int, double>::iterator j;
//		for(j=experiments[k]->getPBegin(); j!=experiments[k]->getPEnd(); j++) {
//			experiments[k]->setP(j->first,
//					(double)((double)experiments[k]->getObservation(j->first) + pseudoCount)/totalCount );
//		}

	UUmap::iterator j;
	for(j=runs[k]->observations.begin(); j!=runs[k]->observations.end(); j++) {
	    //			runs[k]->setP(j->first,
//					(double)((double)experiments[k]->getObservation(j->first) + pseudoCount)/totalCount );
	    
	    runs[k]->p[j->first] = (double)((double)runs[k]->observations[j->first] + param->pseudoCount)/(double)totalCount;
	}
	

	//		experiments[k]->setPNoObs(pseudoCount/(double)experiments[k]->getSize());
	//		experiments[k]->setPNoObs(pseudoCount/totalCount);
    }
}

void SimpleEstimator::estimatePWithNoObs(std::vector<Run*> &runs, unsigned int iterationNumber){

}
