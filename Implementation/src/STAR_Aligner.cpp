/*
 * STAR_Aligner.cpp
 *
 *  Created on: 25.12.2018
 *     Authors: Mario, Willy
 */

#include "../headers/STAR_Aligner.h"
#include "../headers/debug.h"
#include "math.h"

#include "../headers/systemFunctions.h"


using namespace std;

std::string STAR_Aligner::shellCommand(Run *r) {
    /*! \brief
     * Constructs the shell command that calls STAR.
     */

    string cmd;
    string batchDir_ = batchDir(r);

    cmd = param->pathToSTAR + "STAR"
	+ " --runThreadN " + to_string(param->runThreadN)
	+ " --genomeDir " + param->genomeDir;

    if (r->paired == false) { // unpaired reads
	cmd += " --readFilesIn " + batchDir_ + r->accesionId + ".fasta";
    } else { // paired reads
	vector<string> files = filesWithExtInFolder(batchDir_, ".fasta");
	string f = "";
	// we use only the first and the last file, because STAR can only
	// handle two files. For Pombe run SRR097898 the middle part is a technical barcode.
	// So our best guess is to just take reads 1 and 3.
	if (!files.empty()){
	    DEBUG(0,"Found " << files.size() << " Fasta-files!");
	    f = files[0];
	    if (files.size() > 1){
		f += " " + files[files.size()-1];
	    }
	} else {
	    DEBUG(0,"No fasta-files found!");
	}
	cmd += " --readFilesIn " + f;
    }
    cmd += " --outFileNamePrefix " + batchDir_;

    bool useIntronDB = true;
    string intronDBFname = "intronDB.gff"; // see below
    ifstream f(intronDBFname.c_str());
    if (!f.good()) // test for existence
	useIntronDB = false;
    f.close();

    if (useIntronDB)
	cmd += " --sjdbFileChrStartEnd " + intronDBFname;
   
    // + " --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0"
    // + " --outSAMtype BAM";
    return cmd;
}
 
void STAR_Aligner::getAlignedReads(unordered_map<string, RNAread> &reads, Run *r) {
    /*! \brief opens the SAM-file and retrieves 
     * how many times each read mapped to a transcript-block.
     */

    string samFname = batchDir(r) + "Aligned.out.sam";
    string logFname = batchDir(r) + "Log.final.out";
    // First check in Log.final.out if the qualityThreshold is guaranteed
    checkQuality(samFname, logFname, r);

    // string deleteLater = param->outFileNamePrefix + r->accesionId;

    string line;
    std::ifstream samFile(samFname.c_str());

    if (!samFile.is_open()) {
	DEBUG(0,"Unable to open alignment " << samFname);
	// return;
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

		// SRR034413.5	0	FBgn0086917	22037	0	2S89M	*	0	0
		if (wrds.size() > 3) {
		    string readName = wrds[0];
		    string chromosomeName = wrds[2];
		    string mappedAtChromosomPos = wrds[3];

		    // determine which block the read mapped to of chromosome "chromosomeName"
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


void STAR_Aligner::checkQuality(const string &samfilename,
				const string &logfilename,
				Run *r){
    /*
     * evaluate STAR log file
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
    // parse introns
    // count number of lines of introns.gff
    std::ifstream hintsfile(hintsfilename);
    
    if (!hintsfile.is_open()) {
	DEBUG(0,"Unable to open hints file " << hintsfilename);
	exit(1);
    }
    unsigned numintrons = 0;
    while (getline(hintsfile, line))
	numintrons++;
    hintsfile.close();

    double percentSpliced = 100.0 * numintrons / r->batchSize;

    DEBUG(0,"Percent spliced reads   |	" << percentSpliced << "%");
    
    if (r->timesDownloaded>0)
	r->avgSpliced = (r->avgSpliced*(r->timesDownloaded-1) + percentSpliced)/r->timesDownloaded;
    else
	r->avgSpliced = 0.0;
    DEBUG(0, "avgSpliced: " <<  r->avgSpliced);
    
    // accumulate introns
    string cumintronsFname = "cumintrons.gff" ,
	cumintronsTempFname = "cumintrons.temp.gff";
    string intronDBFname = "intronDB.gff"; // see above
    string joincmd = "cat " + hintsfilename;
    ifstream f(cumintronsFname.c_str());
    if (f.good()) // test for existence, otherwise cat gives an error for the very first batch
	joincmd += " " + cumintronsFname;
    f.close();
    joincmd += " | join_mult_hints.pl >" + cumintronsTempFname;
    status = system(joincmd.c_str());
    if (status != 0) {
	DEBUG(0, string("Failed to run join_mult_hints.pl properly:") + joincmd);
	exit(1);
    }
    if (rename(cumintronsTempFname.c_str(), cumintronsFname.c_str()) != 0){
	DEBUG(0, string("Failed to move ") + cumintronsTempFname);
	exit(1);
    }
    
    // filter significantly frequent introns and make intron 'database' for STAR
    std::ifstream cumIntronFile(cumintronsFname.c_str());
    if (!cumIntronFile.is_open()) {
	DEBUG(0,"Unable to open " << cumintronsFname << endl <<
	      "Won't use intron database for aligning next run.");
    }
    std::ofstream intronDBFile(intronDBFname.c_str());
    if (!intronDBFile.is_open())
	DEBUG(0,"Unable to open " << intronDBFname << " for writing.");
    
    bool signif;
    unsigned numsignif = 0; // number of significant introns
    char *gff[7];
    while (getline(cumIntronFile, line)) {
	signif = false;
	if (param->mincovthresh <= 1)
	    signif = true;
	else {
	    int mult;
	    size_t mpos = line.find("mult=");
	    if (mpos != string::npos) {
		mpos += strlen("mult=");
		mult = stoi(line.substr(mpos));
		if (mult >= param->mincovthresh)
		    signif = true;
	    }
	}
	if (signif){
	    numsignif++;
	    char *s = new char[line.length() + 1];
	    strcpy(s, line.c_str());
	    unsigned col = 0;
	    gff[col] = strtok(s, "\t");
	    while (gff[col] && col < 7-1 ) {
		gff[++col] = strtok(NULL, "\t");
	    }
	    if (col >= 7-1){
		unsigned strand = 0; // unknown
		if (gff[6] == "+")
		    strand = 1;
		else if (gff[6] == "-")
		    strand = 2;
		// STAR format: chr <tab> start <tab> end <tab> strand
		intronDBFile << gff[0] << "\t" << gff[3] << "\t"
			     << gff[4] << "\t" << strand << endl;
	    }
	    delete [] s;
	}
    }

    cumIntronFile.close();
    intronDBFile.close();
    DEBUG(0, "number of significant (coverage >= " << param->mincovthresh << ") introns is "
	  << numsignif);
    // TODO:
    // filterIntronsFindStrand.pl " + intronDBFname > introns.s.f.gff
}
