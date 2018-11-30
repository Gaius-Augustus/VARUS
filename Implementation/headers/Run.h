/*
 * Run.h
 *
 *  Created on: 28.11.2016
 *      Author: willy
 */

#ifndef RUN_H_
#define RUN_H_

#include <unordered_map>
#include "ParameterHandler.h"
#include "TypeConventions.h"

class Run {
 public:
    Run();
    Run(std::string accesionId,
	const unsigned int transcriptBlocks,
	const unsigned int numOfSpots,
	const unsigned int batchSize);
    ~Run();

    UUmap observations;
    UDmap p;   // estimator for all entries with an observation
    double pNoObs; // the estimator for all entries with no observations

    unsigned int N; // min-Index for fastq-dump
    unsigned int X; // max-Index for fastq-dump
    std::string accesionId; // <-- name
    bool paired;   // indicates if the run has paired reads

    unsigned int numOfSpots; // number of reads in the run that can be downloaded
    unsigned int maxNumOfBatches; //the number of batches of size batchSize that can be downloaded
    unsigned int sigmaIndex; // Index pointing to the last downloaded batch in sigma
    unsigned int batchSize; // number of reads to be downloaded at once
    std::vector<int> sigma; // holds the order in which the batches should be downloaded from this run
    unsigned int timesDownloaded;

    unsigned int observationSum; // sum of observations of all possible outcomes
    U32 transcriptBlocks;  // number of possible outcomes
    double expectedProfit;
    bool badQuality;
    double avgUmrPercent;
    double avgSpliced;
    double avglen;
    // only for clusterEstimator
    UDmap q;

    //	friend std::ostream& operator<<(std::ostream& os, Run &r);
};

#endif /* RUN_H_ */
