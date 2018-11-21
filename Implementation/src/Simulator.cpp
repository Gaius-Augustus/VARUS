/*
 * Simulator.cpp
 *
 *  Created on: 08.12.2016
 *      Author: willy
 */

#include "../headers/Simulator.h"
#include "../headers/debug.h"
#include <math.h>
#include <memory>
#include <algorithm>    // std::random_shuffle
#include <cstdlib>
#include <ctime>        // std::time
#include "../headers/Operators.h"

#include <random>
#include <stdio.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


using namespace std;

Simulator::Simulator(ParameterHandler *p) {
    param = p;
    ran = new myRandomEngine(p);
}

Simulator::~Simulator() {
    // TODO Auto-generated destructor stub
}

void Simulator::readInputRuns(){
    /*! \brief Reads in the accession-ids of the Runs to download from.
     *
     * It is read from "Runlist.txt"
     *
     * TODO: redundant with readInputRuns() in ChromosomeInitializer
     */

    stringstream st;
    st << param->pathToRuns << "/Runlist.txt";
    std::string line;

    ifstream myfile (st.str().c_str());

    if (!myfile.is_open()){
	DEBUG(0,"Unable to open Runfile:" << st.str());
	param->exit_text();
    }

    while (getline (myfile,line)){
	if(!line.empty() && line[0] != '@'){
	    std::string tmp;
	    std::stringstream ssin(line);
	    std::vector<std::string> wrds;
	    while (ssin >> tmp)	{
		wrds.push_back(tmp);
	    }
	    inputRuns.push_back(wrds[0]);
	}
    }
    myfile.close();

    if(inputRuns.size() < 1) {
	DEBUG(0,"Runlist.txt is empty: " << st.str());
	param->exit_text();
    }
}

void Simulator::initializeRuns(std::vector<Run*> &runs){
    /* example:
     * 	a;10
     * 	b;0
     * 	c;0
     * 	d;5
     *
     * 	names: 	a	b	c	d
     * 	num:	0	1	2	3
     * 	value:	10	0	0	5
     *
     * 			pHidden[0] = 10
     * 			pHidden[1] = 5
     *
     * 			pHiddenTranslate[0] = 0
     * 			pHiddenTranslate[1] = 3
     *
     * 			translate2str[0] = a
     * 			translate2str[1] = b
     * 				...
     * 			translate2str[3] = d
     *
     * 			translate2int[a] = 0
     * 			translate2int[b] = 1
     * 				...
     * 			translate2int[d] = 3
     */

    for(unsigned int j= 0; j < inputRuns.size(); j++) {
	// DEBUG(1,input->InputRuns[j]);
	
	string runName = inputRuns[j];
	string name = param->pathToDice;
	name += "/"+runName;
	name += ".csv";

	const char * c = name.c_str();
	string line;
	ifstream file (c);

	if (!file.is_open()) {
	    DEBUG(0,"Can't read dice: Unable to open file:" << name);
	    DEBUG(0,"You can either try to delete " << name <<
		  " in the Runlist.txt or add " << name <<
		  " to the correct folder. Exiting programm");
	    exit(-1);
	}
	bool first = true;
	
	vector<double> pHidden;
	unsigned int numOfReads = 0;
	
	unsigned int num = 0;
	while (getline (file,line)){
	    if(line.size() > 1)	{
		string tmp;
		stringstream ssin(line);
		
		vector<string> wrds;
		
		while (ssin){
		    string s;
		    if (!getline( ssin, s, ';' )) break;
		    wrds.push_back( s );
		}

		// in a line with an '@' we have the number of reads
		if(wrds[0][0] == '@'){
		    numOfReads = std::stoi(wrds[1]);
		} else {
		    //blockname;observations
		    translate2str[num] = wrds[0];
		    translate2int[wrds[0]] = num;
		    
		    pHidden.push_back(atoll(wrds[1].c_str()));
		    
		    num++;
		}
	    }
	}
	
	// makes pHidden sum to 1
	normVector(pHidden);
	DEBUG(3,pHidden);
	DEBUG(3,pHidden.size());
	
	Run *r = new Run(runName,param->numOfBlocks,numOfReads,param->batchSize);
	runs.push_back(r);
	
	//		p_hidden.push_back(pHidden);
	p_hidden[runName] = pHidden;
	
    }

	return;
}

void Simulator::simulateObservations(Run* r, UUmap &totalObservations){

    std::discrete_distribution<> d(p_hidden[r->accesionId].begin(),p_hidden[r->accesionId].end());

    for(int n=0; n < param->batchSize; ++n) {
    	unsigned int ind = d(ran->gen);
        ++r->observations[ind];
        ++r->observationSum;
        ++totalObservations[ind];
    }
    r->timesDownloaded++;
}
