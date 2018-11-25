/*
 * Aligner.h
 *
 *  Created on: 30.11.2016
 *      Author: willy
 */

#ifndef ALIGNER_H_
#define ALIGNER_H_
#include "ParameterHandler.h"
#include "Run.h"
#include "RNAread.h"
#include "ChromosomeInitializer.h"

class Aligner {
    /*! \brief This class controlls the alignment-steps for the runs.
     *
     * The allignemnt is done with STAR.
     */

 public:
    ParameterHandler *param;
    
    Aligner(ParameterHandler *p);
    virtual ~Aligner();

    std::string shellCommand(Run *r);

    void mapReads(Run *r);

    void getAlignedReads(std::unordered_map<std::string, RNAread> &reads, Run *r);

    void updateObservations(Run *r, UUmap &totalObservations,
			    std::unordered_map<std::string, RNAread> &reads, ChromosomeInitializer *c);

    void update(Run *r, UUmap &totalObservations, ChromosomeInitializer *c);

    void checkQuality(const std::string &s, Run *r);
};

#endif /* ALIGNER_H_ */
