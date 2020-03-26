/*
 * Downloader.cpp
 *
 *  Created on: 29.11.2016
 *      Author: willy
 */

#include "../headers/Downloader.h"
#include "../headers/debug.h"


Downloader::Downloader(ParameterHandler *p) {
    param = p;
}

Downloader::~Downloader() {
}

std::string Downloader::shellCommand(Run *r){
    std::ostringstream n, x;
    
    n << r->N;
    x << r->X;
    
    std::string s = param->fastqDumpCall + " -N " + n.str() + " -X " + x.str() + " -O " + param->outFileNamePrefix
        + param->aliDirName + "/" + r->accesionId + "/"
	+ "N" + n.str() + "X" + x.str() + "/ --fasta 120 ";
    if (r->paired)
	s += "--split-files ";
    
    s += r->accesionId;
    return s;
}

void Downloader::nextBatchIndices(Run *r) {
    /*! \brief	Calculates the batch-indices N and X.
     */
    //	assert(sigmaIndex < sigma.size());
    r->N = r->sigma[r->sigmaIndex] * r->batchSize;
    r->X = r->N + r->batchSize - 1;
    
    if (r->X > r->numOfSpots)
	r->X = r->numOfSpots;

    r->sigmaIndex++;
}

bool Downloader::getBatch(Run *r, bool all){
    /*! \brief	The next batch will be downloaded according to the sigma-vector.
     * 	The command for the download is invoked through the shell and calls fastq-dump.
     *
     *	Sometimes batches are not loaded correctly. This is a problem with fastq-dump.
     *	In this case the reads can not be aligned. This function returns "false", if the
     *	download failed.
     */
    unsigned numberOfAttempts = 2;
    unsigned attempt = 0;
    
    if (all == false)
	nextBatchIndices(r);	// chooses the correct N and X
    else {
	r->N = 0;
	r->X = r->numOfSpots;
    }
    
    std::string s = shellCommand(r);
    
    std::ostringstream msg;
    msg << "getBatch(" << r->accesionId << ", X=" << r->X << ", N=" << r->N << "): " << s;
    DEBUG(0,msg.str());
    int status = -1;

    while (status != 0 && ++attempt <= numberOfAttempts){
	status = system(s.c_str()); // run fastq-dump
	if (status != 0) {
	    DEBUG(0,"Failed to save batch "  << r->accesionId << " -N  " << r->N << " -X " << r->X <<
		  (attempt < numberOfAttempts? ". Retrying..." : "") );
	}
    }
    if (status != 0){
	DEBUG(0,"Definitively failed to save batch "  << r->accesionId << " -N  " << r->N << " -X " << r->X);
	return false;
    }

    r->timesDownloaded++;
    
    return true;
}
