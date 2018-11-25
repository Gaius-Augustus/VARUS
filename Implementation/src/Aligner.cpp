/*
 * Aligner.cpp
 *
 *  Created on: 30.11.2016
 *      Author: willy
 */

#include "../headers/Aligner.h"
//#include "../headers/Operators.h"
#include "../headers/debug.h"
#include "math.h"

#include "../headers/systemFunctions.h"


using namespace std;

Aligner::Aligner(ParameterHandler *p) {
    param = p;
}

Aligner::~Aligner() {
}

std::string Aligner::shellCommand(Run *r) {
    /*! \brief Constructs the shell command that calls STAR.
     *
     */

    std::ostringstream n,x,thread;
    n << r->N;
    x << r->X;
    thread << param->runThreadN;

    std::string s;
    if(r->paired == false) {
	std::string t = r->accesionId + "/" + "N" + n.str() + "X" + x.str() + "/";
	s =      param->pathToSTAR + "STAR "
	    + "--runThreadN " + thread.str() + " "
	    + "--genomeDir " + param->genomeDir + " "
	    + "--readFilesIn " + param->outFileNamePrefix + t + r->accesionId + ".fasta "
	    + "--outFileNamePrefix " + param->outFileNamePrefix + t
	    //							+ " --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0"
	    //	                        + " --outSAMtype BAM";
	    ;
    } else { // unpaired reads
	std::string t = param->outFileNamePrefix + "/" + r->accesionId + "/" + "N" + n.str() + "X" + x.str() + "/";

	vector<string> files = filesWithExtInFolder(t,".fasta");

	string f = "";
	//	        	for(unsigned int i = 0; i < files.size(); i++){
	//	        		 f += files[i] + " ";
	//	        	}

	// we use only the first and the last file, because STAR can apparently only
	// handle two files. For Pombe run SRR097898 the middle part is a technical barcode.
	// So our best guess is to just take reads 1 and 3.
	if(!files.empty()){
	    DEBUG(0,"Found " << files.size() << " Fasta-files!");
	    f = files[0];
	    if(files.size() > 1){
		f += " " + files[files.size()-1];
	    }
	} else {
	    DEBUG(0,"No fasta-files found!");
	}

	s = param->pathToSTAR + "STAR "
	    + "--runThreadN " + thread.str() + " "
	    + "--genomeDir " + param->genomeDir + " "
	    + "--readFilesIn "
	    + f
	    + " --outFileNamePrefix " + t
	    // + " --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0"
	    // + " --outSAMtype BAM";
	    ;
    }
    return s;
}

void Aligner::mapReads(Run *r) {
    /*! \brief Calls STAR through the shell. The reads of the last batch downloaded from the passed run are then
     *  mapped to the transcript-units according to STAR.
     *
     */

    // build the command for STAR
    std::string s = shellCommand(r);

    const char * c = s.c_str();
    int status = system(c);
    if (0 != status) {
	DEBUG(0,"Failed to run STAR properly!");
    }
}

void Aligner::getAlignedReads(unordered_map<string, RNAread> &reads, Run *r) {
	/*! \brief opens the SAM-file and retrieves how many times each read mapped to a transcript-block.
	 *
	 */

	std::ostringstream n, x;
	n << r->N;
	x << r->X;
	std::string st = param->outFileNamePrefix + r->accesionId + "/" + "N" + n.str()
			+ "X" + x.str() + "/Aligned.out.sam";
//    cerr << __FILE__ <<':'<<__LINE__<<':'<< "path;" << s << endl;
//    string s = input->outFileNamePrefix + "Aligned.out.sam";

//	string deleteLater = outFileNamePrefix + r->getAccesionId() + "/" + "N"
//			+ n.str() + "X" + x.str();

	// First check in Log.final.out if the qualityThreshold is guaranteed
	std::string st2 = param->outFileNamePrefix + r->accesionId + "/" + "N" + n.str()
					+ "X" + x.str() + "/Log.final.out";
	checkQuality(st2, r);


	std::string deleteLater = param->outFileNamePrefix + r->accesionId;

	const char * c = st.c_str();
	std::string line;
	std::ifstream myfile(c);

	if (!myfile.is_open()) {
		DEBUG(0,"Unable to open alignment " << st);
		// return;
	}

	while (getline(myfile, line)) {
		if (!line.empty()) {
			if (line[0] != '@') {
				string tmp;
				stringstream ssin(line);
				vector<string> wrds;
				while (ssin >> tmp) {
					wrds.push_back(tmp);
				}

				// SRR034413.5	0	FBgn0086917	22037	0	2S89M	*	0	0
				if (wrds.size() > 3) {
					string readName = wrds[0];
					string chromosomeName = wrds[2];
					string mappedAtChromosomPos = wrds[3];

					// determine which block the read mapped to of chromosom "chromosomeName"
					// mapping-pos/blockSize rounded down
					std::ostringstream oss;
					oss
							<< floor(
									atoll(mappedAtChromosomPos.c_str())
											/ param->blockSize);
					string chromosomBlock = chromosomeName + ":" + oss.str();

					// ToDo: cigar-string
					if (reads.find(readName) != reads.end()) {
						// add the transcript_unit-name to the transcript_unit-list of the read
						if (reads[readName].transcriptUnits.find(chromosomBlock)
								!= reads[readName].transcriptUnits.end()) {
							reads[readName].transcriptUnits[chromosomBlock] =
									reads[readName].transcriptUnits[chromosomBlock]
											+ 1;
						} else {
							reads[readName].transcriptUnits[chromosomBlock] = 1;
						}
					} else {
						RNAread mp;
						mp.transcriptUnits[chromosomBlock] = 1;
						reads[readName] = mp;
					}
				}
			}
		}
	}
	myfile.close();

//	if(param->deleteLater == 1){
//		string remove = "rm -rf " + deleteLater;
//		const char * c2 = remove.c_str();
//		int status = system(c2);
//
//		if (0 != status) {
//			DEBUG(0, "\nFailed to delete:" << deleteLater);
//		}
//	}

	return;
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

void Aligner::checkQuality(const std::string &st, Run *r){
    const char * c = st.c_str();
    std::string line;
    std::ifstream myfile(c);

    if (!myfile.is_open()) {
	DEBUG(0,"Unable to open " << st);
	r->badQuality = true;
	// return;
    }

    while (getline(myfile, line)) {
	if (line.find("Uniquely mapped reads % |") != std::string::npos) {
	    string tmp;
	    stringstream ssin(line);
	    vector<string> wrds;
	    while (ssin >> tmp) {
		wrds.push_back(tmp);
	    }

	    std::string numWithPercent = wrds[wrds.size()-1];
	    std::string nWoP = numWithPercent.substr(0,numWithPercent.size()-1);
	    double quality = std::stof(nWoP);

	    DEBUG(0,"Uniquely mapped reads % |	" << nWoP << "%");

	    if("-nan" == nWoP){
		DEBUG(0,"Detected run with no matches.");
		r->badQuality = true;

	    } else if(param->qualityThreshold > quality){
		DEBUG(0,"Detected run with poor alignability.");
		r->badQuality = true;
	    }
	    // compute average unique mappability over all batches from this run
	    if (r->timesDownloaded>0)
		r->avgUmrPercent = (r->avgUmrPercent*(r->timesDownloaded-1)+quality)/r->timesDownloaded;
	    else
		r->avgUmrPercent = 0.0;
	    DEBUG(0, "Computing avgUmrPercent: " <<  r->avgUmrPercent << "\t" << quality << "\t" << r->timesDownloaded);
	}
    }
    myfile.close();
}
