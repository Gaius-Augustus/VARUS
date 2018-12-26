/*
 * Aligner.h
 *
 *  Created on: 30.11.2016
 *      Author: willy
 */

#ifndef STAR_ALIGNER_H_
#define STAR_ALIGNER_H_
#include "Aligner.h"
#include "ParameterHandler.h"
#include "Run.h"
#include "RNAread.h"
#include "ChromosomeInitializer.h"

class STAR_Aligner : public Aligner {
    /*! \brief 
     *
     * STAR aligner
     */
 public:
    STAR_Aligner(ParameterHandler *p) : Aligner(p) {}
    std::string shellCommand(Run *r);
    void getAlignedReads(std::unordered_map<std::string, RNAread> &reads, Run *r);
    void checkQuality(const std::string &sam,
		      const std::string &log,
		      Run *r);
};

#endif /* STAR_ALIGNER_H_ */
