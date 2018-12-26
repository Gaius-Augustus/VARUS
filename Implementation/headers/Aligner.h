/*
 * Aligner.h
 *
 *  Created on: 30.11.2016
 *      Authors: willy, Mario
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
     * The allignemnt is done with STAR or HISAT.
     */

 public:
    ParameterHandler *param;
    
    Aligner(ParameterHandler *p);
    ~Aligner();

    virtual std::string shellCommand(Run *r) = 0;

    void mapReads(Run *r);

    virtual void getAlignedReads(std::unordered_map<std::string, RNAread> &reads, Run *r) = 0;

    void updateObservations(Run *r, UUmap &totalObservations,
			    std::unordered_map<std::string, RNAread> &reads,
			    ChromosomeInitializer *c);

    void update(Run *r, UUmap &totalObservations, ChromosomeInitializer *c);

};

#endif /* ALIGNER_H_ */
