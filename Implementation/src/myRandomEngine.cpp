/*
 * myRandomEngine.cpp
 *
 *  Created on: 06.12.2016
 *      Author: willy
 */

#include "../headers/myRandomEngine.h"

myRandomEngine::myRandomEngine(ParameterHandler *p) {
    /*! \brief Seed is set randomly depending on the time or according to the value
     * randomSeed in ParameterHandler param.
     */
    
    param = p;
    if (p->randomSeed > 0){
	gen.seed(p->randomSeed);
    } else {
	gen.seed(time(0));
    }
}

myRandomEngine::myRandomEngine() {
    /*! \brief Seed is set randomly depending on the time.
     */
    
    param = nullptr;
    gen.seed(time(0));
}

myRandomEngine::~myRandomEngine() {
    // TODO Auto-generated destructor stub
}

unsigned int myRandomEngine::select_randomly(unsigned int max){
    std::uniform_int_distribution<> dis(0, max);
    int index = dis(gen);
    return index;
}

unsigned int myRandomEngine::select_randomly_p(std::vector<double> &vec){
    std::discrete_distribution<> dis(vec.begin(), vec.end());
    int index = dis(gen);
    return index;
}


void myRandomEngine::seed(unsigned int t){
    /*! \brief Set the seed of the engine. Only needed for unit-testing to make
     * things reproduce-able.
     */
    gen.seed(t);
}


unsigned myRandomEngine::biasSelect(const std::vector<double> &scores, 
				    unsigned bestOfK){
    unsigned choice = 0;
    double bestScore = std::numeric_limits<double>::min();
    for (unsigned i=0; i < bestOfK && i< scores.size(); ++i){
	unsigned idx = select_randomly(scores.size()-1);
	if (scores[idx] > bestScore){
	    bestScore = scores[idx];
	    choice = idx;
	}
    }
    return choice;
}
