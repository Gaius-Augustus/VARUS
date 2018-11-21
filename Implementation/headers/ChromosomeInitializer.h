/*
 * ChromosomInitializer.h
 *
 *  Created on: 28.11.2016
 *      Author: willy
 */

#ifndef CHROMOSOMEINITIALIZER_H_
#define CHROMOSOMEINITIALIZER_H_
#include "Run.h"
#include "ParameterHandler.h"
#include <memory>

#include "myRandomEngine.h"

class ChromosomeInitializer {
    /*! \brief This class initializes the runs.
     *
     * - Sigma order is set with this class. That is the order in which the batches should be downloaded.
     * - A translation-convention from UUmap to SUmap is created.
     * UUmaps are used most of the time in this implementation, because they use less memory.
     *
     */

 public:
    ParameterHandler *param;
    myRandomEngine * ran;

    SUmap chromosomLengths;		// holds the lengths of the chromosomes

    SUmap transcriptUnitsByNames;	// has size ~= sum(chromosomLengths)/blocksize

    USmap transcriptUnits;

    SUmap translate2int;

    USmap translate2str;

    std::vector<U32> order;

    std::vector<std::string> inputRuns; // names of the given runs


    ChromosomeInitializer(ParameterHandler *p);
    ~ChromosomeInitializer();

    void readInputRuns();	// reads in the list of runs that should be downloaded from

    void getChromosomeLengths();

    void create_transcript_units(UUmap &totalObservations);		// initializes the total observations in
    // controler as well

    //	void translateTranscriptNames2ints(SUmap &transcriptUnitsByNames,
    //										USmap &transcriptUnits);

    void initializeRuns(std::vector<Run*> &runs);
    
    //	void initializeSigma(Run* r, int t = -1);
    void initializeSigma(Run* r);
};

#endif /* CHROMOSOMINITIALIZER_H_ */
