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
    /*! \brief This class controlls the alignment steps for the runs.
     *
     * The allignemnt is done with STAR or HISAT.
     */

 public:
    ParameterHandler *param;
    
    Aligner(ParameterHandler *p);
    ~Aligner();

    virtual std::string shellCommand(Run *r) = 0;

    int mapReads(Run *r); // returns 0 iff successful

    virtual void getAlignedReads(std::unordered_map<std::string,
				 RNAread> &reads,
				 Run *r,
				 int batchNr) = 0;
    void sam2transcriptUnits(std::unordered_map<std::string, RNAread> &reads,
			     std::string samFname);
    
    void updateObservations(Run *r, UUmap &totalObservations,
			    std::unordered_map<std::string, RNAread> &reads,
			    ChromosomeInitializer *c);

    void update(Run *r, UUmap &totalObservations, ChromosomeInitializer *c, int batchNr);

    std::string batchDir(Run *r);

    // filter significantly frequent introns and make intron 'database' for STAR or HISAT
    // return number of (intron) hints, including multiplicity
    int filterGTF(std::string gtfInFileName, std::string intronDBFname="", unsigned thresh=0);

    void cleanupAfterAlignment(Run *r);
};

#endif /* ALIGNER_H_ */
