/*
 * myRandomEngine.h
 *
 *  Created on: 06.12.2016
 *     Authors: willy, mario
 */

#ifndef MYRANDOMENGINE_H_
#define MYRANDOMENGINE_H_

#include <algorithm>    // std::random_shuffle
#include <cstdlib>
#include <ctime>        // std::time
#include <random>
#include <algorithm>

#include "ParameterHandler.h"
#include "TypeConventions.h"

class myRandomEngine {
public:
    std::mt19937 gen; // Mersenne Twister

    ParameterHandler *param;
    myRandomEngine();
    myRandomEngine(ParameterHandler *p);
    virtual ~myRandomEngine();

    template<typename V>
	V select_randomly(std::vector<V> &vec){
	return vec[select_randomly(vec.size()-1)];
    }

    template<typename V>
	void shuffleExceptLast(std::vector<V> &vec){
	if (vec.size() > 1){
	    typename std::vector<V>::iterator i = vec.end();
	    i--;
	    std::shuffle(vec.begin(), i, gen);
	}
    }

    template<typename V>
	void shuffle(std::vector<V> &vec){
	std::shuffle(vec.begin(), vec.end(), gen);
    }

    unsigned int select_randomly(unsigned int max);
    unsigned int select_randomly_p(std::vector<double> &vec);
    void seed(unsigned int t);

    /* 
     * Randomly select, but introduce a bias for a high score.
     * Chose an item with highest score from bestOfK items
     * chosen uniformly.
     */
    unsigned biasSelect(const std::vector<double> &scores,
			unsigned bestOfK = 4);
};


#endif /* MYRANDOMENGINE_H_ */
