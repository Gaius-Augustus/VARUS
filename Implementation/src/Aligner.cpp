/*
 * Aligner.cpp
 *
 *  Created on: 30.11.2016
 *     Authors: willy, Mario
 */

#include "../headers/Aligner.h"
#include "../headers/debug.h"
#include "math.h"

#include "../headers/systemFunctions.h"


using namespace std;

Aligner::Aligner(ParameterHandler *p) {
    param = p;
}

Aligner::~Aligner() {
}


void Aligner::mapReads(Run *r) {
    /*! \brief Calls aligner through the shell. The reads of the last downloaded batch 
     * are mapped to the transcript-units.
     */

    // build the command
    std::string s = shellCommand(r);

    const char * c = s.c_str();
    int status = system(c);
    if (0 != status) {
	DEBUG(0, "Failed to run aligner properly!");
    }
}

void Aligner::updateObservations(Run *r, UUmap &totalObservations,
				 unordered_map<string, RNAread> &reads, ChromosomeInitializer * c){
    /* !\brief if a read is mapped to only one chromosomblock and only once, the number
     * of reads mapped to that chromosomblock is increased by one
     * Note: this frequency.size() != total_freq_Dist.size()
     *       that means unmapped chromosomblocks will not be initialized
     *
     */

    string chromosom_block;
    U32 UMRcount = 0;
    unordered_map<string, RNAread>::iterator it;
    for (it = reads.begin(); it != reads.end(); ++it) {
	// evaluates if the read is an UMR to a block not UMR to a Chromosome

	if (it->second.UMR()) {
	    UMRcount++;
	    chromosom_block = it->second.transcriptUnits.begin()->first;

	    DEBUG(4, chromosom_block << " is UMR");
	    //			r->observationSum++;

	    U32 key = c->translate2int[chromosom_block];

	    // tmp_run

	    //			r->addObservation(key);
	    if(r->observations.find(key) != r->observations.end())
		{
		    r->observations[key] += 1;
		}
	    else
		{
		    r->observations[key] = 1;
		}

	    ++r->observationSum;

	    // totalObservations
	    if (totalObservations.find(key) != totalObservations.end()) {
		totalObservations[key] += 1;
	    } else {
		totalObservations[key] = 1;
	    }
	} else {
	    DEBUG(4, chromosom_block << " is not UMR");
	}
    }
//	NULLreads = index - UMRcount;
}

void Aligner::update(Run *r, UUmap &totalObservations, ChromosomeInitializer *c){
    /* !\brief This function creates an alignment, gets the uniquely mapped reads and updates the run
     * and totalObservation accordingly.
     */

    mapReads(r);
    unordered_map<string,RNAread> reads;
    getAlignedReads(reads,r);
    updateObservations(r,totalObservations,reads,c);
}



std::string Aligner::batchDir(Run *r) {
     /* !\brief The directory of the next batch of the run, e.g. SRR7192719/N6300000X6349999/
     */
    return  param->outFileNamePrefix + r->accesionId + "/"
	+ "N" + to_string(r->N)
	+ "X" + to_string(r->X) + "/";
}

