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
    string alignCmd = shellCommand(r);

    int status = system(alignCmd.c_str());
    if (0 != status)
	DEBUG(0, "Failed to run aligner properly: " << alignCmd);
}

void Aligner::updateObservations(Run *r, UUmap &totalObservations,
				 unordered_map<string, RNAread> &reads, ChromosomeInitializer * c){
    /* !\brief if a read is mapped to only one chromosom_block and only once, the number
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
	RNAread &read = it->second;
	if (read.UMR()) {
	    UMRcount++;
	    chromosom_block = read.transcriptUnits.begin()->first;
	    
	    DEBUG(4, chromosom_block << " is UMR");

	    U32 key = c->translate2int[chromosom_block];

	    if (r->observations.find(key) != r->observations.end())
		r->observations[key] += 1;
	    else
		r->observations[key] = 1;

	    ++r->observationSum;

	    // totalObservations
	    if (totalObservations.find(key) != totalObservations.end())
		totalObservations[key] += 1;
	    else
		totalObservations[key] = 1;
	} else {
	    DEBUG(2, chromosom_block << " is not UMR");
	}
    }
    DEBUG(4, "r->observationSum = " << r->observationSum << " UMRcount = " << UMRcount);
}

void Aligner::update(Run *r, UUmap &totalObservations, ChromosomeInitializer *c, int batchNr){
    /* !\brief This function creates an alignment, gets the uniquely mapped reads
     * and updates the run and totalObservation accordingly.
     */

    mapReads(r);
    unordered_map<string,RNAread> reads;
    getAlignedReads(reads, r, batchNr);
    updateObservations(r, totalObservations, reads, c);
    cleanupAfterAlignment(r);
}


void Aligner::sam2transcriptUnits(unordered_map<string, RNAread> &reads,
				  string samFname){
    string line;
    ifstream samFile(samFname.c_str());

    if (!samFile.is_open()) {
	DEBUG(0,"Unable to open alignment " << samFname);
	return;
    }

    while (getline(samFile, line)) {
	if (!line.empty()) {
	    if (line[0] != '@') {
		string tmp;
		stringstream ssin(line);
		vector<string> wrds;
		while (ssin >> tmp) {
		    wrds.push_back(tmp);
		}

		// SRR034413.5	0   FBgn0086917	22037	0  2S89M  *  0	0
		if (wrds.size() > 3) {
		    string readName = wrds[0];
		    string chromosomeName = wrds[2];
		    string mappedAtChromosomPos = wrds[3];

		    // determine which block the read mapped to of "chromosomeName"
		    // mapping-pos / blockSize rounded down
		    std::ostringstream oss;
		    oss << floor(atoll(mappedAtChromosomPos.c_str()) / param->blockSize);
		    string chromosomBlock = chromosomeName + ":" + oss.str();

		    // ToDo: cigar-string
		    if (reads.find(readName) != reads.end()) {
			// add the transcript_unit-name to the transcript_unit-list of the read
			if (reads[readName].transcriptUnits.find(chromosomBlock)
			    != reads[readName].transcriptUnits.end()) {
			    reads[readName].transcriptUnits[chromosomBlock] += 1;
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
    samFile.close();
}

string Aligner::batchDir(Run *r) {
     /* !\brief The directory of the next batch of the run, e.g. SRR7192719/N6300000X6349999/
     */
    return  param->outFileNamePrefix + r->accesionId + "/"
	+ "N" + to_string(r->N)
	+ "X" + to_string(r->X) + "/";
}

int Aligner::filterGTF(string gtfInFileName, string intronDBFname, unsigned thresh){
    bool write = (intronDBFname != "");
    unsigned numintrons = 0; // mult not ignored, return value
    ofstream intronDBFile;
    
    std::ifstream cumIntronFile(gtfInFileName.c_str());
    if (!cumIntronFile.is_open()) {
	DEBUG(0,"Unable to open " << gtfInFileName << endl);
	return -1;
	// "Won't use intron database for aligning next run."
    }
    if (write){
	intronDBFile.open(intronDBFname.c_str());
	if (!intronDBFile.is_open())
	    DEBUG(0,"Unable to open " << intronDBFname << " for writing.");
    }
    
    bool signif;
    unsigned numsignif = 0; // number of significant introns, mult ignored
    char *gff[7];
    string line;
    while (getline(cumIntronFile, line)) {
	signif = false;
	int mult = 1;
	size_t mpos = line.find("mult=");
	if (mpos != string::npos) {
	    mpos += strlen("mult=");
	    mult = stoi(line.substr(mpos));
	}
	numintrons += mult;
	if (mult >= param->mincovthresh){
	    signif = true;
	    numsignif++;
	}
	
	if (signif && write){
	    char *s = new char[line.length() + 1];
	    strcpy(s, line.c_str());
	    unsigned col = 0;
	    gff[col] = strtok(s, "\t");
	    while (gff[col] && col < 7-1)
		gff[++col] = strtok(NULL, "\t");

	    if (col >= 7-1){
		unsigned strand = 0; // unknown
		if (*gff[6] == '+')
		    strand = 1;
		else if (*gff[6] == '-')
		    strand = 2;
		// STAR format: chr <tab> start <tab> end <tab> strand(0,1,2)
		intronDBFile << gff[0] << "\t" << gff[3] << "\t"
			     << gff[4] << "\t" << strand << endl;
	    }
	    delete [] s;
	}
    }

    cumIntronFile.close();
    if (write){
	intronDBFile.close();
	DEBUG(0, "number of significant (coverage >= " << thresh << ") introns is " << numsignif);
    }
    return numintrons;
}

void Aligner::cleanupAfterAlignment(Run *r){
    // zip .fasta files
    string batchDir_ = batchDir(r);
    string zipCmd = "gzip";
    zipCmd += " " + batchDir_ + "*.fasta";
    DEBUG(0, "cleanupAfterAlignmen: " << zipCmd);
    int status = system(zipCmd.c_str());
    if (status != 0)
	DEBUG(0, "Failed to run gzip properly: " << zipCmd);

    // delete .sam files
    string delCmd = "rm";
    delCmd += " " + batchDir_ + "*.sam";
    DEBUG(0, "cleanupAfterAlignmen: " << delCmd);
    status = system(delCmd.c_str());
    if (status != 0)
	DEBUG(0, "Failed to delete .sam file: " << delCmd);
    
}
