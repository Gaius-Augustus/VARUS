/*
 * HISAT_Aligner.cpp
 *
 *  Created on: 25.12.2018
 *     Authors: Mario, Willy
 */

#include "../headers/HISAT_Aligner.h"
#include "../headers/debug.h"
#include "math.h"

#include "../headers/systemFunctions.h"


using namespace std;

std::string HISAT_Aligner::shellCommand(Run *r) {
    /*! \brief
     * Constructs the shell command that calls HISAT.
     */

    string batchDir_ = batchDir(r);

    string cmd = "hisat";
    cmd += " -p " + to_string(param->runThreadN)
	+ " -x " + param->genomeDir + "/hisatidx";
    
    if (r->paired == false) { // unpaired reads
	cmd += " -U " + batchDir_ + r->accesionId + ".fasta";
    } else { // paired reads
	vector<string> files = filesWithExtInFolder(batchDir_, ".fasta");
	// we use only the first and the last file.
	// For Pombe run SRR097898 the middle part is a technical barcode.
	// So our best guess is to just take reads 1 and 3.
	if (files.empty()){
	    DEBUG(0, "No fasta files found!");
	    exit(1);
	}
	if (files.size() < 2){
	    DEBUG(0, "Paired run but found only one fasta files!");
	    exit(1);
	}
	cmd += " -1 " + files[0] + " -2 " + files[files.size()-1];
    }
    cmd += " -S " + batchDir_ + "Alignment.out.sam";

    bool useIntronDB = true;
    string intronDBFname = "intronDB.gff"; // see below
    ifstream f(intronDBFname.c_str());
    if (!f.good()) // test for existence
	useIntronDB = false;
    f.close();

    if (useIntronDB)
	cmd += " --known-splicesite-infile " + intronDBFname;
   
    return cmd;
}
 
void HISAT_Aligner::getAlignedReads(unordered_map<string, RNAread> &reads, Run *r, int batchNr) {
    /*! \brief opens the SAM file and retrieves 
     * how many times each read mapped to a transcript-block.
     */

    string samFname = batchDir(r) + "Aligned.out.sam";
    string logFname = batchDir(r) + "Log.final.out";
    // check log files if the qualityThreshold is guaranteed
    checkQuality(samFname, logFname, r, batchNr);

    sam2transcriptUnits(reads, samFname);
}


void HISAT_Aligner::checkQuality(const string &samfilename,
				const string &logfilename,
				Run *r, int batchNr){
    /*
     * evaluate HISAT log file
     */
    string batchDir_ = batchDir(r);
    string line;
    std::ifstream logfile(logfilename.c_str());
    
    if (!logfile.is_open()) {
	DEBUG(0,"Unable to open " << logfilename);
	r->badQuality = true;
    }
    
    while (getline(logfile, line)) {
	if (line.find("Uniquely mapped reads number |") != string::npos) {
	    string tmp;
	    stringstream ssin(line);
	    vector<string> wrds;
	    while (ssin >> tmp) {
		wrds.push_back(tmp);
	    }
	    
	    string numberstr = wrds[wrds.size()-1];
	    double quality = 100 * std::stof(numberstr) / r->batchSize;
	    
	    DEBUG(0,"Uniquely mapped reads % |	" << quality << "%");
	    
	    if("-nan" == numberstr){
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
	    DEBUG(0, "Computing avgUmrPercent: " <<  r->avgUmrPercent);
	}
    }
    logfile.close();

    /*
     * evaluate alignment SAM file
     */
    std::ifstream samfile(samfilename.c_str());
    
    if (!samfile.is_open()) {
	DEBUG(0,"Unable to open " << samfilename);
	r->badQuality = true;
    } else {
	samfile.close();
    }
    // sort and convert to BAM
    size_t lastindex = samfilename.find_last_of("."); 
    string stem = samfilename.substr(0, lastindex);
    string bamfilename = stem + ".bam";
    string sortcmd = "samtools sort " + samfilename + " -O BAM -o " + bamfilename;
    int status = system(sortcmd.c_str());
    if (status != 0) {
	DEBUG(0, string("Failed to run samtools sort properly:") + sortcmd);
	exit(1);
    }
    
    // make introns with bam2hints
    string hintsfilename = batchDir_ + "introns.gff";
    string hintscmd = "bam2hints --in=" + bamfilename
	+ " --out=" + hintsfilename + " --intronsonly";
    status = system(hintscmd.c_str());
    if (status != 0) {
	DEBUG(0, string("Failed to run bam2hints properly:") + hintscmd);
	exit(1);
    }
    // count number of introns in introns.gff (with multiplicity)
    int numintrons = filterGTF(hintsfilename);
    if (numintrons<0) {
	DEBUG(0,"Could not open " << hintsfilename);
	numintrons = 0; // ignore
    }
	
    double percentSpliced = 100.0 * numintrons / r->batchSize;

    DEBUG(0,"Percent spliced reads   |	" << percentSpliced << "%");
    
    if (r->timesDownloaded>0)
	r->avgSpliced = (r->avgSpliced*(r->timesDownloaded-1) + percentSpliced)/r->timesDownloaded;
    else
	r->avgSpliced = 0.0;
    DEBUG(0, "avgSpliced: " <<  r->avgSpliced);
    
    // accumulate introns
    string cumintronsFname = "cumintrons.gff" ,
	cumintronsTempFname = "cumintrons.temp.gff",
    	cumintronsStrandedFname = "cumintrons.stranded.gff";
    string intronDBFname = "intronDB.gff"; // see above
    string joincmd = "cat " + hintsfilename;
    ifstream f(cumintronsFname.c_str());
    if (f.good()) // test for existence, otherwise cat gives an error for the very first batch
	joincmd += " " + cumintronsFname;
    f.close();
    joincmd += " | join_mult_hints.pl >" + cumintronsTempFname;
    status = system(joincmd.c_str());
    if (status != 0) {
	DEBUG(0, string("Failed to run join_mult_hints.pl properly: ") + joincmd);
	exit(1);
    }
    if (rename(cumintronsTempFname.c_str(), cumintronsFname.c_str()) != 0){
	DEBUG(0, string("Failed to move ") + cumintronsTempFname);
	exit(1);
    }

    // The intron DB update is not done every time as the increments can be small.
    unsigned period = 500000 / param->batchSize; // update every half million reads
    if (period < 1)
	period = 1;
    if (batchNr % period == 0){
	// determine strands of introns
	string strandCmd = "filterIntronsFindStrand.pl";
	strandCmd += " " + param->genomeFaDir + " " + cumintronsFname + " > " + cumintronsStrandedFname;
	status = system(strandCmd.c_str());
	if (status != 0){
	    DEBUG(0, string("Failed to run filterIntronsFindStrand.pl properly: ")
		  + strandCmd);
	    exit(1);
	}
	// make 4 column intron DB for aligers
	filterGTF(cumintronsStrandedFname, intronDBFname, param->mincovthresh);
    }
}
