/*
 * Simulator.h
 *
 *  Created on: 08.12.2016
 *      Author: willy
 */

#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include "TypeConventions.h"
#include "ParameterHandler.h"
#include "Run.h"
#include "myRandomEngine.h"

#include "debug.h"
#include <string>
#include <cstdlib>

/*! \brief This class simulates the download and allignment-steps.
 *
 * The simulation basically consists of producing vectors according to
 * a random-distribution of the runs. These vectors are then interpreted
 * as actual observations of reads mapping to transcripts.
 *
 */

class Simulator {
public:
	ParameterHandler *param;
	myRandomEngine	 *ran;

	std::vector<std::string> inputRuns;
	std::unordered_map<std::string,std::vector<double>> p_hidden; 	// store a p_hidden for each run


    SUmap translate2int;

    USmap translate2str;

    std::vector<U32> order;


	Simulator(ParameterHandler *p);
	~Simulator();

	void readInputRuns();
	void initializeRuns(std::vector<Run*> &runs);
	void simulateObservations(Run* r, UUmap &totalObservations);
};

#endif /* SIMULATOR_H_ */
