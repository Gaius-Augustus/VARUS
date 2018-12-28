/*
 * HISAT_Aligner.h
 *
 *  Created on: 30.11.2016
 *     Authors: willy, Mario
 */

#ifndef HISAT_ALIGNER_H_
#define HISAT_ALIGNER_H_
#include "Aligner.h"
#include "ParameterHandler.h"
#include "Run.h"
#include "RNAread.h"
#include "ChromosomeInitializer.h"

class HISAT_Aligner : public Aligner {
    /*! \brief 
     * HISAT aligner
     */
 public:
    HISAT_Aligner(ParameterHandler *p) : Aligner(p) {}
    std::string shellCommand(Run *r);
    void getAlignedReads(std::unordered_map<std::string, RNAread> &reads, Run *r, int batchNr);
    void checkQuality(const std::string &sam,
		      const std::string &log,
		      Run *r,
		      int batchNr);
};

#endif /* HISAT_ALIGNER_H_ */
