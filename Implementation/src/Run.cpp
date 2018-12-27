/*
 * Run.cpp
 *
 *  Created on: 28.11.2016
 *      Author: willy
 */

#include "../headers/Run.h"
#include "../headers/Operators.h"
#include "../headers/debug.h"
#include <algorithm>    // std::random_shuffle
#include <cstdlib>
#include <ctime>        // std::time
#include <assert.h>

using namespace std;

extern ParameterHandler *param;
//const int verbosity_debug_level =-3;
//const int verbosity_out_level =3;

Run::Run() {
    paired = false;
    N = 0;
    X = 0;
    numOfSpots = 0;
    maxNumOfBatches = 0;
    sigmaIndex = 0;
    batchSize = 0;
    timesDownloaded = 0;
    transcriptBlocks = 0;
    observationSum = 0;
    pNoObs = 0.0;
    expectedProfit = 0.0;
    accesionId = "noId";
    badQuality = false;
    avgUmrPercent = 0.0;
    avgSpliced = 0.0;
    pRep = nullptr;
}

Run::Run(std::string accesionId, const unsigned int transcriptBlocks, const unsigned int numOfSpots,
	 const unsigned int batchSize) {
    this->accesionId = accesionId;
    paired = false;
    N = 0;
    X = 0;
    this->numOfSpots = numOfSpots;
    //	assert(batchSize != 0);
    DEBUG(2, "Creating run with numOfSpots: " << numOfSpots << ", batchSize: " << batchSize);
    maxNumOfBatches = ceil((double)numOfSpots/(double)batchSize);
    sigmaIndex = 0;
    this->batchSize = batchSize;
    timesDownloaded = 0;
    this->transcriptBlocks = transcriptBlocks;
    observationSum = 0;
    pNoObs = 0.0;
    expectedProfit = 0.0;
    badQuality = false;
    avgUmrPercent = 0.0;
    avgSpliced = 0.0;
    pRep = nullptr;
}

Run::~Run() {
    // TODO Auto-generated destructor stub
}
