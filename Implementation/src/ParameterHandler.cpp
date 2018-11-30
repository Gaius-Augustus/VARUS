/*
 * ParameterHandler.cpp
 *
 *  Created on: 25.11.2016
 *      Author: willy
 */


#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "../headers/vname.h"
#include "../headers/debug.h"
#include "../headers/Warning.h"

#include "../headers/ParameterHandler.h"
#include "../headers/Operators.h"
#include <iterator>
#include <string.h>
#include <iomanip>

#include "assert.h"

using namespace std;

int getoptCount = 0;


#ifndef PARAM
#define PARAM(x) {          \
		string name = #x;	\
    add_parameter(name, x); \
}
#endif

#ifndef INITPARAM
#define INITPARAM(x,z,use,dep) {          \
		string name = #x;	\
    add_parameter(name, x, z, use,dep); \
}
#endif

std::string ParameterHandler::printTextWithWidthAndLeftOffset(const std::string s, const unsigned int l, const unsigned int off){

	assert(l > off);

    std::string tmp;
    std::stringstream ssin(s);
    std::vector<std::string> wrds;
    while (ssin >> tmp)
    {
        wrds.push_back(tmp);
    }

    unsigned int whitespaceSize = 1;
    unsigned int ind = 0;

    unsigned int curr_length = off;
    string ret = "";
	for(unsigned int j = 0; j < off; j++){
		ret += " ";
	}

	for(unsigned int i = 0; i < wrds.size(); i++){
		// adding the next words leads to stepping over the border
		// +2 because of "\n"
		if(curr_length + wrds[i].size() + whitespaceSize> l){
			ret += "\n";
			curr_length = off;

			for(unsigned int j = 0; j < off; j++){
				ret += " ";
			}
		}

		ret += wrds[i] + " ";
		curr_length += wrds[i].size() + whitespaceSize;
	}

//	unsigned int times = floor(s.size()/l);
//
//	for(unsigned int i = 1; i <= times; i++){
//		s.insert(l*i,"\n\t\t");
//	}
//
//	return s;

	return ret;
}

string ParameterHandler::lineLength(string s, const unsigned int l, const unsigned int maxS){
    std::string tmp;
    std::stringstream ssin(s);
    std::vector<std::string> wrds;
    while (ssin >> tmp)
    {
        wrds.push_back(tmp);
    }


    unsigned int whitespaceSize = 1;

    unsigned int ind = 0;

    string par = wrds[0];

    string val;
    if(wrds.size() >=2){
    	val = wrds[1];
    }

	string ret = par;

	for(unsigned int j = 0; j < maxS-par.size(); j++){
		ret += " ";
	}

	ret += val + "\n";

	for(unsigned int j = 0; j < maxS-par.size(); j++){
		ret += " ";
	}

    unsigned int curr_length = maxS;

	for(unsigned int i = 2; i < wrds.size(); i++){
		// adding the next words leads to stepping over the border
		// +2 because of "\n"
		if(curr_length + wrds[i].size() + whitespaceSize> l){
			ret += "\n";
			for(unsigned int j = 0; j < maxS; j++){
				ret += " ";
			}
			curr_length = maxS;
		}
		ret += " " + wrds[i];
		curr_length += wrds[i].size() + whitespaceSize;
	}


//	unsigned int times = floor(s.size()/l);
//
//	for(unsigned int i = 1; i <= times; i++){
//		s.insert(l*i,"\n\t\t");
//	}
//
//	return s;

	return ret;
}

void ParameterHandler::print_usage() {
	/*! \brief Prints the usage of the program.
	 *
	 * The map<string,string> usage holds the name and its description for each parameter.
	 */

    string us = "Usage:\n"
    		"Online-algorithm to automatically draw optimal samples from libraries of \n"
    		"RNA-seqdata from the NCBI with the aim to achieve a coverage as high as \n"
    		"possible, with as little data being downloaded as possible.\n"
    		"In a stepwise procedure parts of these libraries are downloaded and aligned \n"
    		"to the Genome of the species with STAR. At each step the run that is expected \n"
    		"to yield the most increase in coverage, is chosen. \n\n"
    		"Input: \n"
    		"- Genome of the target species specified with --genomeDir\n"
    		" The Genome files have to be prepared for STAR. After indexing the Folder should \n"
    		" contain the following files: \n"
    		" chrLenghth.txt, chrName.txt, chrNameLength.txt, chrStart.txt, Genome, \n"
    		" genomeParameters.txt, SA, SAindex \n\n"
    		+STARmanual+ "\n\n"
    		"- Runlist named \"Runlist.txt\". \n"
    		"Runlist-Format: \n"
    		"@Run_acc total_spots total_bases bool:paired #tabulator separated \n"
    		"ERR1328564	75912627	1	5041721263	0 \n"
    		"ERR1328563	88959042	1	7616693616	0 \n"
    		"DRR030358	20853417	1	750723012	0 \n"
    		"					...						\n\n"
    		"An exemplary call of VARUS might look like this:\n\n"
    		"./VARUS --blockSize 5000 --batchSize 10000 --genomeDir <...> --pathToStar <...>\n"
    		" --pathToRuns <...> --outFileNamePrefix <...> \n\n"
    		"We divide the genome in many smaller blocks.This is done because on some larger \n"
    		"transcripts the coverage is also desired to be evenly distributed.\n"
    		"At each step a fixed number of reads is downloaded. We call a set of reads \n"
    		"downloaded at once a batch. The batchSize specifies how many reads are \n "
    		"downloaded at once.";

    cout << us << endl;



    for(unsigned int i = 0; i < lineWidth; i++){
        cout << "-";
    }
    cout << "\n";



//    map<string,string>::iterator i;
    unsigned int maxParLength = 0;
//    for(i = usage.begin(); i != usage.end(); i++) {
//    	if(i->first.size() > maxParLength) maxParLength = i->first.size();
//    }

    map<string,Parameter>::iterator i;

    for(i = parameters.begin(); i != parameters.end(); i++) {
    	if(i->first.size() > maxParLength) maxParLength = i->first.size();
    }

    // for "--" and ":"
    maxParLength += 3;

    // First print the most important parameters
//	cout << "Paths and basic settings:\n";
//	cout << "-----------------------------------------------------\n";

    printParameterCategory(MANDATORY, "",maxParLength);

//    // paths and STAR-parameters
//    printParameterCategory(BASICSETTINGS,  maxParLength);

    printParameterCategory(BASICSETTINGS,
    						"If you choose to not use an estimator, runs are downloaded randomly"
    						"until you have reached maxBatches.",
							maxParLength);
//
//    // estimators
    printParameterCategory(ESTIMATOR,
    						"In each step of the algorithm the distributions of the reads in each run "
    						"are estimated based on the reads that are downloaded already of this run."
    						" We introduced a cost-function that assigns higher scores to runs that "
    						"are expected to contain reads mapping to regions in the genome with only "
    						"little observations at the given step. The download stops if the "
    						"expected profit is lower than the cost for downloading another batch. "
    						"The default of cost = 0 lets VARUS download maxbatches reads. "
    						"With cost > 0 the download can stop prematurely. "
    						"A sensible choice could be cost = 1.",
							maxParLength);

    printParameterCategory(SIMPLEANDADVANCED,
    						"The simple estimator adds pseudocounts to the actual observations. The "
    						"distribution of the reads is then estimated with the maximum-likelihood-method."
    						" The advanced estimator adds the frequencies of the sum of all observations "
    						"in all runs. This sum is multiplied with a weight lambda to control how "
    						"similar the runs are to be assumed.",
							maxParLength);

    printParameterCategory(DIRICHLETMIXTURE, "", maxParLength);
    printParameterCategory(CLUSTERESTIMATOR, "", maxParLength);
//
//    // simulation
    printParameterCategory(SIMULATION, "", maxParLength);

    printParameterCategory(COMMANDLINE,
    						" Instead of typing all commands in the console you can read in and write to a file "
    						" that stores all your settings. This file can also be read in for future program-runs"
    						" to ensure to use the same settings.",
							maxParLength);

//	cout << parameterCategories[BASICSETTINGS];
//	cout << ":\n-----------------------------------------------------\n";
//
//	for(i = parameters.begin(); i != parameters.end(); i++) {
//		if(i->second.category == parameterCategories[BASICSETTINGS]){
//			string out = "--" + i->first + ": " + i->second.category + " " + i->second.description;
//			cout << lineLength(out,80,maxParLength);
//			cout << "\n" << endl;
//		}
//    }

//	cout << "Advanced settings and modelparameters:\n";
//	cout << "-----------------------------------------------------\n";
//
//	for(i = parameters.begin(); i != parameters.end(); i++) {
//		if(i->second.category == "2"){
//			string out = "--" + i->first + ": " + i->second.category + " " + i->second.description;
//			cout << lineLength(out,80,maxParLength);
//			cout << "\n" << endl;
//		}
//    }
}

void ParameterHandler::printParameterCategory(const paramCat cat, const std::string des, const unsigned int maxParLength){
	/*
	 * Prints the usage for the given parameter.
	 */



	bool first = true;
	map<string,Parameter>::iterator i;

	const char separator    = ' ';

	for(i = parameters.begin(); i != parameters.end(); i++) {
		if(i->second.deprecated == false && i->second.category == parameterCategories[cat]){

			// only print the category-description if there are any parameters in it
			if(first){
				cout << parameterCategories[cat] << ":\n";;
				cout << printTextWithWidthAndLeftOffset(des,80,maxParLength) << "\n";


				for(unsigned int i = 0; i < lineWidth; i++){
					cout << "-";
				}
				cout << "\n";

				first = false;
			}

			string out = "--" + i->first + ": " + i->second.value;

			cout << lineLength(out,lineWidth,maxParLength) << "\n";

			string out2 = i->second.description;
			cout << printTextWithWidthAndLeftOffset(out2,lineWidth,maxParLength);
//
////			cout.width(maxParLength + 2 + 1 + 1);
////			cout << "--" << i->first << ":";
////			cout.fill(' ');
////			cout.width(lineWidth);
////			cout << i->second.value << "\n";
////			cout << std::setw(10) << endl;
//
//			string par = "--" + i->first + ":";
//		    cout << left << setw(maxParLength) << setfill(separator) << par;
//		    cout << left << setw(lineWidth-maxParLength) << setfill(separator) << i->second.value << "\n\n";
//
////		    cout << left << setw(maxParLength) << setfill(separator) << " ";
//
//		    // print all lines of the description
//		    for(unsigned int j = 0; j < i->second.description.size(); j += lineWidth-maxParLength){
//		    	cout << left << setw(maxParLength) << setfill(separator) << " ";
//		    	cout << left << setw(lineWidth-maxParLength) << setfill(separator) << i->second.description.substr(j,lineWidth-maxParLength) << "\n";
//
//		    }
//
			cout << "\n" << endl;
		}
    }

    for(unsigned int i = 0; i < lineWidth; i++){
        cout << "-";
    }
    cout << "\n";

}

void ParameterHandler::exit_text(){
	cout << "Try ./VARUS --help for more information." << endl;
	exit(1);
}

ParameterHandler::ParameterHandler() {
	/*! \brief All parameters are initialized with default-values. To change these values call
	 * readArguments().
	 */

	// initialize the parameterCategories first

	parameterCategories[MANDATORY] 			= "Mandatory Inputs";
	parameterCategories[BASICSETTINGS] 		= "Basic Settings";
	parameterCategories[DIRICHLETMIXTURE] 	= "Dirichlet Mixture";
	parameterCategories[SIMULATION] 		= "Simulation";
	parameterCategories[SIMPLEANDADVANCED] 	= "Simple and advanced Estimator";
	parameterCategories[CLUSTERESTIMATOR] 	= "Cluster Estimator";
	parameterCategories[ESTIMATOR] 			= "Estimators";
	parameterCategories[COMMANDLINE]		= "Commandline handling";

	parameterCategories[TEST] 	= "Test";

	STARmanual = "https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf";

//	pathToVARUS = "/home/willy/BRAKER/VARUS/Implementation/";

	varusID = "(V)";

	param = this;

	mergeThreshold = 2;
	INITPARAM(mergeThreshold,
		  parameterCategories[BASICSETTINGS],
		  "Specifies after how many downloads a merge of all the "
		  "existing alignments should be done. Note that after "
		  "the alignments are merged all subfolders with smaller"
		  " alignments are deleted. At the end one final merge is done. "
		  "2 is the minimal value. "
		  "Note: If you set the parameter deleteLater = 0, then only"
		  " one merge at the end of execution is done.",false);

	qualityThreshold = 5.0;
	INITPARAM(qualityThreshold,
		  parameterCategories[BASICSETTINGS],
		  "Specifies how much of uniquely mapping reads you are allowing."
		  " STAR produces a file 'Log.final.out' with information of "
		  "how much percent of the reads did map. Really small values "
		  "like 0.1 can be an indicator of intrinsic dna.",false);

	fastqDumpCall = "fastq-dump";
	INITPARAM(fastqDumpCall,
		  parameterCategories[MANDATORY],
		  "if fastq-dump is not in your $PATH, you can specify "
		  "the path as well like /path/to/fastqdump/./fastq-dump."
		  "Note that this string has to contain './fastq-dump' as well.",false);

	//"/home/willy/Bachelorarbeit/STAR-2.5.1b/source/"
	pathToSTAR    = "/home/willy/Bachelorarbeit/STAR-2.5.1b/source/";
    	INITPARAM(pathToSTAR,
		  parameterCategories[MANDATORY],
		  "specifies the path to the executable of STAR, the software that"
		  " aligns the reads to the transcriptome/genome.",false);

	pathToVARUS = "/home/willy/BRAKER/VARUS/Implementation/";;
	INITPARAM(pathToVARUS,
		  parameterCategories[MANDATORY],
		  "specifies the path to the executable of VARUS",false);

	//	"/home/willy/workspace/UnitTests/ParameterHandler/script/"
	pathToRuns 			= "../Tutorial/Drosophila/";
    	INITPARAM(pathToRuns,
		  parameterCategories[MANDATORY],
		  "specifies the path to Runlist.txt."
		  "You can retrieve a Runlist with the perl-script 'RunListRetriever.pl'.",false);

	//	readFilesIn 		= "default";
	//    	INITPARAM(readFilesIn, "specifies the path to the folder in which the ");

	//    	"/home/willy/Bachelorarbeit/Manager/genome/"
	genomeDir 			= "/home/willy/Bachelorarbeit/Manager/genome/";
    	INITPARAM(genomeDir,
		  parameterCategories[MANDATORY],
		  "specifies the path to the genome/transcriptome"
		  " that you are investigating.",false);

	outFileNamePrefix 	= "../Tutorial/Drosophila/";
    	INITPARAM(outFileNamePrefix,
		  parameterCategories[MANDATORY],
		  "specifies the path in which all output of VARUS should be stored.",false);

	runThreadN = 4;
    	INITPARAM(runThreadN,
		  parameterCategories[BASICSETTINGS],
		  "Number of threads to run STAR parallel with. Read STAR-manual for more information.\n",false);

	blockSize = 5000;		// <---- 5000 default
    	INITPARAM(blockSize,
		  parameterCategories[BASICSETTINGS],
		  "the number of bases one block will have. "
		  "This is done in order to be able"
		  " to compare the coverage of larger chromosomes with smaller ones, since "
		  "larger chromosomes will naturally have more reads mapped to them than smaller ones.",false);

	batchSize = 100000;
    	INITPARAM(batchSize,
		  parameterCategories[BASICSETTINGS],
		  "the number of reads to be downloaded in each step. "
		  "Typically a library contains about 10e+6 to 10e+7 reads.",false);

	pseudoCount = 1;
    	INITPARAM(pseudoCount,
		  parameterCategories[SIMPLEANDADVANCED],
		  "adds a pseudocount to all possible observations. Only relevant for "
		  "estimators 1, 2.",false);

//    coverage = 10;
//    	PARAM(coverage);
//    zValue = 20;
//    	PARAM(zValue);
    lambda = 10.0;
    	INITPARAM(lambda,
    			parameterCategories[SIMPLEANDADVANCED],
				"parameter for estimator 2. The weight with which the sum of all observations"
				" in all runs is added as additional pseudocounts for the estimation.",false);

//    logScore = 1;	// 1 == logScore, 0 == minScore
//    	PARAM(logScore);
    createDice = 0;		// 1 => creating dice, 0 => other options
    	INITPARAM(createDice,
    			parameterCategories[SIMULATION],
				"if set to 1 dice will be created from real runs."
    			"Only needed for testing/simulation-purposes.",true);

	simulation = 0;		// 1 => running simulation, 0 => real downloadprogramm
    	INITPARAM(simulation,
    			parameterCategories[SIMULATION],
				"if set to 1, the program will simulate downloads.",true);

	estimator = 2;	// 1 == simple, 2 == advanced, 3 == DM, else naive
    	INITPARAM(estimator,
		  parameterCategories[BASICSETTINGS],
		  "1 == simple, 2 == advanced, 3 == Dirichlet mixture, 4 == cluster estimator, otherwise downloads will "
		  "be done choosing the runs randomly. Note that if you choose to download randomly, you "
		  "should specify maxBatches in order to let the program end at some point.",false);

//	dieList = 1;	// 1== uses "Dielist.txt", 0 == uses "Runlist.txt
//    	INITPARAM(dieList, "1== uses \"Dielist.txt\", 0 == uses \"Runlist.txt\"");

    lessInfo = 0; 	// 1 == Toy prints less info
    	INITPARAM(lessInfo,
    			parameterCategories[SIMULATION],
				"if set to 1, Simulator prints less info. Only Relevant for simulation.",true);

    pathToDice 	= "/home/willy/Bachelorarbeit/Manager/test138/";
    	INITPARAM(pathToDice,
    			parameterCategories[SIMULATION],
				"specifies the path to the dice. Dice are saved in csv-format.",true);

    cost = 0;
    	INITPARAM(cost,
    			parameterCategories[ESTIMATOR],
				"sets the cost for downloading one read."
    			" The cost-function is a monotonically increasing function with monotonically decreasing derivative.",false);

    components = 1;
    	INITPARAM(components,
    			parameterCategories[DIRICHLETMIXTURE],
				"sets the number of components of the dirichlet mixture or the cluster estimator. "
    			"Only relevant for estimators 3 and 4.",true);

    numOfBlocks = 29919;
    	INITPARAM(numOfBlocks,
    			parameterCategories[BASICSETTINGS],
    			"sets the number of blocks into which the genome should be divided."
				"To calculate numOfBlocks each transcript is divided by blockLength rounded up. "
				"The sum of these numbers is numOfBlocks."
    			"This is done because naturally it is expected that larger transcripts have a higher chance"
    			" that reads are mapping to them. Also the coverage among the larger transcripts is desired "
    			"to be more or less uniformly distributed.",true);

    trainingsIterations = 1;
    	INITPARAM(trainingsIterations,
    			parameterCategories[DIRICHLETMIXTURE],
				"the number of iterations the dirichlet mixture or cluster estimator will be "
    			"trained in each step. Only relevant for estimators 3 and 4.",true);

    loadAllOnce = 0;
    	INITPARAM(loadAllOnce,
    			parameterCategories[ESTIMATOR],
				"if set to 1, a single batch from each run will be downloaded once, before"
    			" the estimation process starts.",false
    			);

	verbosityDebug = 0;
		INITPARAM(verbosityDebug,
				parameterCategories[BASICSETTINGS],
				"sets the debug-verbosity-level. There are 4 verbosity-levels. Higher values mean more output.",false);
//		verbosity_out_level = verbosityDebug;

	newtonIterations = 10;
		INITPARAM(newtonIterations,
				parameterCategories[DIRICHLETMIXTURE],
				"number of times the newtons method is done to find the maximum-likelihood"
				" for the alpha-sum. Only relevant for estimator 3.",true);

	newtonPrecision = 0.001;
		INITPARAM(newtonPrecision,
				parameterCategories[DIRICHLETMIXTURE],
				"threshold at which the newtons-method will aboard, and return the value."
				" Only relevant for estimator 3.",true);

	readParametersFromFile = 1;
		INITPARAM(readParametersFromFile,
				parameterCategories[COMMANDLINE],
				"if set to 1, the program will look for a file specified with "
				"pathToParameters and interpret its content as command-line arguments. "
				"These parameters will then be used to run the program. Note: additional parameters"
				" passed with the command line will overwrite the parameters read from the parameters-file.",false);
			readAllready = false;

	currentWorkingDirectory = get_working_path();

	pathToParameters = currentWorkingDirectory + "/VARUSparameters.txt";
		INITPARAM(pathToParameters,
				parameterCategories[COMMANDLINE],
				"specifies the path and name to the parameters-file that should be read in and written to. "
				"Default is the current working directory and a file called 'VARUSparameters.txt'.",false);

	exportParametersToFile = 0;
		INITPARAM(exportParametersToFile,
				parameterCategories[COMMANDLINE],
				"if set to 1 the parameters used for this execution of the program"
				" will be exported into a parametersfile",false);

	exportObservationsToFile = 0;
		INITPARAM(exportObservationsToFile,
				parameterCategories[BASICSETTINGS],
				"if set to 1 the program will output the observations in all"
				" runs in all steps into CSV-files. This does not refer to the fasta.files"
				" but to the quantitative distributions of the reads among the genome.\n"
				"WARNING: This will overwrite older files! \n"
				"NOTE: Using this option can lead to performance-issues.",false);

	deleteLater = 0;
		INITPARAM(deleteLater,
				parameterCategories[BASICSETTINGS],
				"if set to 1, the fasta-files and alignment-files will be deleted after"
				" they are used to identify the next run to be downloaded from. If you want to use the "
				"reads for your genome-annotation you should not use this option.",false);

	maxBatches = 100;
		INITPARAM(maxBatches,
				parameterCategories[BASICSETTINGS],
				"if  maxBatches > 0, the program will exit as soon as it loaded maxBatches batches.",false);

	randomSeed = -1;
		INITPARAM(randomSeed,
				parameterCategories[BASICSETTINGS],
				"if randomSeed > 0, the program will have deterministic results. Else the "
				"seed will be set according to the current time",false);

	profitCondition = 0;
		INITPARAM(profitCondition,
				parameterCategories[ESTIMATOR],
				"if profitConditon == 1, the program will exit if the"
				" expected profit falls below 0. Note that the expected profit can lead to the program downloading"
				" for a very long time, since some of the estimators tend to be very optimistic "
				"if the parameters are not set adequately.",false);

	ignoreReadNum = 0;
		INITPARAM(ignoreReadNum,
				parameterCategories[SIMULATION],
				"Only important for simulation: if ignoreReadNum == 1, it will be ignored in case of a simulation "
				"if a run has no reads left. Otherwise the run is not an option to download from after the maximum number of reads"
				" is downloaded from this run.",true)

	simpleDM = 0;
		INITPARAM(simpleDM,
				parameterCategories[DIRICHLETMIXTURE],
				"refers to estimator 3. With this estimation procedure the calculation times are a bit better "
				"than with the normal dirichlet mixture. However the estimation is not that accurate.",true);

//	kMeansIterations = 10;
//		INITPARAM(kMeansIterations, "only relevant for estimator 4: Number of iterations for "
//				"the k-means-cluster-algorithm.");

	exportNewtons = 0;
		INITPARAM(exportNewtons,
				parameterCategories[DIRICHLETMIXTURE],
				"exports the steps of the newtons method. Deprecated since it is very expensive."
				" Only relevant for estimator 3.",true);


//	param.reset(this);
}
ParameterHandler::~ParameterHandler() {
	cerr << "parameterHandler deleted" << endl;
}

void ParameterHandler::readArguments(int argc, char *argv[]) {
	/*! \brief Reads the command line arguments.
	 */

    const char* const short_opts = "h";
    const option long_opts[] = {
	{"batchSize", 1, nullptr, 'b'},
	{"blockSize", 1, nullptr, 'B'},
	{"components", 1, nullptr, 'c'},
	{"cost", 1, nullptr, 'C'},
	{"createDice", 1, nullptr, 'D'},
	//			{"dieList", 1, nullptr, 'L'},
	{"estimator", 1, nullptr, 'e'},
	{"genomeDir", 1, nullptr, 'g'},
	{"lambda", 1, nullptr, 'l'},
	{"lessInfo", 1, nullptr, 'I'},
	{"loadAllOnce", 1, nullptr, 'A'},
	{"numOfBlocks", 1, nullptr, 'O'},
	{"outFileNamePrefix", 1, nullptr, 'F'},
	{"pathToDice", 1, nullptr, 'T'},
	{"pathToRuns", 1, nullptr, 'R'},
	{"pathToSTAR", 1, nullptr, 'S'},
	{"fastqDumpCall", 1, nullptr, 'f'},
	{"pseudoCount", 1, nullptr, 'p'},
	{"runThreadN", 1, nullptr, 'N'},
	{"simulation", 1, nullptr, 's'},
	{"trainingsIterations", 1, nullptr, 't'},
	{"verbosityDebug", 1, nullptr, 'v'},
	{"newtonIterations", 1, nullptr, 'i'},
	{"newtonPrecision", 1, nullptr, 'P'},
	{"readParametersFromFile", 1, nullptr, 'q'},
	{"pathToParameters", 1, nullptr, 'Q'},
	{"exportParametersToFile", 1, nullptr, 'Z'},
	{"exportObservationsToFile", 1, nullptr, 'z'},
	{"deleteLater", 1, nullptr, 'd'},
	{"maxBatches", 1, nullptr, 'm'},
	{"randomSeed", 1, nullptr, 'a'},
	{"profitCondition", 1, nullptr, 'E'},
	{"ignoreReadNum", 1, nullptr, 'G'},
	{"simpleDM", 1, nullptr, 'H'},
	//			{"kMeansIterations", 1, nullptr, 'J'},
	{"exportNewtons", 1, nullptr, 'j'},
	{"help", 0, nullptr, 'h'},
	{"mergeThreshold", 1, nullptr, 'w'},
	{"pathToVARUS", 1, nullptr, 'W'},
	{"qualityThreshold", 1, nullptr, 'V'},
	{nullptr, 0, nullptr, 0}
    };

	optind = 0;

    while (true)
    {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt)
        {
        case 'b':
            batchSize = std::stoi(optarg);
            	PARAM(batchSize);
            break;

        case 'B':
        	blockSize = std::stoi(optarg);
        		PARAM(blockSize);
            break;

        case 'c':
        	components = std::stoi(optarg);
        		PARAM(components);
            break;

        case 'C':
        	cost = std::stof(optarg);
        		PARAM(cost);
            break;

        case 'D':
            createDice = std::stoi(optarg);
            	PARAM(createDice);
            break;

//        case 'L':
//            dieList = std::stoi(optarg);
//            	PARAM(dieList);
//            break;

        case 'e':
            estimator = std::stoi(optarg);
            	PARAM(estimator);
            break;

        case 'g':
            genomeDir = std::string(optarg);
            	PARAM(genomeDir);
            break;

        case 'l':
            lambda = std::stof(optarg);
            	PARAM(lambda);
            break;

        case 'I':
            lessInfo = std::stoi(optarg);
            	PARAM(lessInfo);
            break;

        case 'A':
            loadAllOnce = std::stoi(optarg);
            	PARAM(loadAllOnce);
            break;

        case 'O':
            numOfBlocks = std::stoi(optarg);
            	PARAM(numOfBlocks);
            break;

        case 'F':
            outFileNamePrefix = std::string(optarg);
            	PARAM(outFileNamePrefix);
            break;

        case 'T':
            pathToDice = std::string(optarg);
            	PARAM(pathToDice);
            break;

        case 'R':
            pathToRuns = std::string(optarg);
            	PARAM(pathToRuns);
            break;

        case 'S':
        	pathToSTAR = std::string(optarg);
            	PARAM(pathToSTAR);
            break;

        case 'f':
        	fastqDumpCall = std::string(optarg);
            	PARAM(fastqDumpCall);
            break;

        case 'p':
            pseudoCount = std::stof(optarg);
            	PARAM(pseudoCount);
            break;

        case 'N':
            runThreadN = std::stoi(optarg);
            	PARAM(runThreadN);
            break;

        case 's':
            simulation = std::stoi(optarg);
            	PARAM(simulation);
            break;

        case 't':
            trainingsIterations = std::stoi(optarg);
            	PARAM(trainingsIterations);
            break;

        case 'v':
            verbosityDebug = std::stoi(optarg);
            	PARAM(verbosityDebug);
//            	verbosity_out_level = verbosityDebug;
            break;

        case 'i':
            newtonIterations = std::stoi(optarg);
            	PARAM(newtonIterations);
            break;

        case 'P':
        	newtonPrecision = std::stof(optarg);
            	PARAM(newtonPrecision);
            break;

        case 'q':
            readParametersFromFile = std::stoi(optarg);
            	PARAM(readParametersFromFile);
            break;

        case 'Q':
            pathToParameters = std::string(optarg);
            	PARAM(pathToParameters);
            break;

        case 'Z':
            exportParametersToFile = std::stoi(optarg);
            	PARAM(exportParametersToFile);
            break;

        case 'z':
            exportObservationsToFile = std::stoi(optarg);
            	PARAM(exportObservationsToFile);
            break;

        case 'd':
            deleteLater = std::stoi(optarg);
            	PARAM(deleteLater);
            break;

        case 'm':
            maxBatches = std::stoi(optarg);
            	PARAM(maxBatches);
            break;

        case 'a':
            randomSeed = std::stoi(optarg);
            	PARAM(randomSeed);
            break;

        case 'E':
            profitCondition = std::stoi(optarg);
            	PARAM(profitCondition);
            break;

        case 'G':
            ignoreReadNum = std::stoi(optarg);
            	PARAM(ignoreReadNum);
            break;

        case 'H':
            simpleDM = std::stoi(optarg);
            	PARAM(simpleDM);
            break;

//        case 'J':
//            kMeansIterations = std::stoi(optarg);
//            	PARAM(kMeansIterations);
//            break;

        case 'j':
        	exportNewtons = std::stoi(optarg);
            	PARAM(exportNewtons);
            break;

        case 'w':
        	mergeThreshold = std::stoi(optarg);
        		if(mergeThreshold < 2){
        			DEBUG(0,"WARNING: Setting mergeThreshold = 2.")
        			mergeThreshold = 2;
        		}
            	PARAM(mergeThreshold);
            break;

        case 'V':
        	qualityThreshold = std::stof(optarg);
        		if(qualityThreshold < 0.0){
        			DEBUG(0,"WARNING: Setting qualityThreshold = 0.1")
        			mergeThreshold = 0.1;
        		}
        		else if(qualityThreshold > 100.0){
        			DEBUG(0,"WARNING: Setting qualityThreshold = 100.0")
        			mergeThreshold = 100.0;
        		}
            	PARAM(qualityThreshold);
            break;

        case 'W':
        	pathToVARUS = std::string(optarg);
            	PARAM(pathToVARUS);
            break;

        case 'h': // -h or --help
			print_usage();
			exit(-1);
        	break;
        case '?': // Unrecognized option
        	print_usage();
        	exit(1);
        	break;
        default:
            print_usage();
            break;
        }
    }


    /**
     *  if the user selects the simulation, no real runs will be used
     */

//    if(simulation == 0) {
//    	dieList = 0;
//    }

    /**
     *  If the user selected a path readParameters from file, we have to call
     *  read_parameters_from_file(), which in turn will call readArguments.
     *  So we need to set a flag in order to not get stuck in a loop.
     */

    if(readParametersFromFile != 0)
    {
		if(readAllready == false)
		{
			argc_saved = argc;
			argv_saved = new char*[argc_saved];

			for(unsigned int i = 0; i < argc; i++) {
				argv_saved[i] = argv[i];
			}


			readAllready = true;
			read_parameters_from_file(pathToParameters);
		}
    }

}

std::ostream & ParameterHandler::printParameters(std::ostream& os) {
	/*! \brief The parameters are appended to the std::ostream os.
	 *
	 * The format is in such a way, that a call of read_parameters_from_file()
	 *  with an output-file created with this method can be read in again.
	 * For example: Let lambda=1.2
	 *
	 * outFile:
	 * .
	 * .
	 * .
	 * --lambda 1.2
	 * .
	 * .
	 * .
	 */

//	map<string,string>::iterator it;
	map<string,Parameter>::iterator it;
	for(it = parameters.begin(); it != parameters.end(); it++) {
//		os << "--" << it->first << " " << it->second << "\n";
//		os << "--" << it->first << " " << it->second.description << "\n";

		if(it->second.deprecated == false){
			os << "--" << it->first << " " << it->second.value << "\n";
		}
	}
	return os;
}

void ParameterHandler::export_parameters() {
	/*! \brief The parameters are exported into a file.
	 *
	 * The parameters can be read in from this file with read_parameters_from_file().
	 */

	if(exportParametersToFile != 0)
	{
		DEBUG(0,"Exporting parameters to " << pathToParameters);

		ofstream f;
//		f.open(pathToParameters);

		f.open(pathToParameters, ios::out);
		if(!f) { DEBUG(0,"Cant export to file " << pathToParameters);}


		ostringstream ss;

		printParameters(ss);
		f << ss.str();
//		cout << ss.str() << endl;

		f.close();
	}
}

template<typename K>
void ParameterHandler::add_parameter(string name, K a, string category, string use, bool deprecated) {
	/*! \brief This function is called when a new parameter should be added.
	 *
	 * With this function and the map parameters it is possible to
	 * import and export the values of the parameters easily.
	 */

    stringstream ss;
    ss << a;
    string value = ss.str();

//	string value = std::to_string(a);

//    parameters[name] = value;


    parameters[name].value = value;
    parameters[name].name = name;

    if(category != "")
    	parameters[name].category = category;

    if(use != "")
    	parameters[name].description = use;

    parameters[name].deprecated = deprecated;
}

void ParameterHandler::read_parameters_from_file(string path) {
	/*! \brief The parameters are read from a file.
	 *
	 * The parameters can be written to a file with exportParameters().
	 */

//	DEBUG(1,"reading from " << path);
	ifstream filestr;
	filestr.open (path);
	if (!filestr.is_open())
	{
		DEBUG(0,"Unable to open parameters-file " << path);
		if(path.empty()) DEBUG(0,"path is empty!");
	}


	std::vector<std::string> wrds;
	wrds.push_back("./RNAM");

	string line;
	while (getline (filestr,line))
	{
		if(!line.empty() && line.size() > 1)
		{
			if(line[0] == '-' && line[1] == '-') {
				std::string tmp;
				std::stringstream ssin(line);

				while (ssin >> tmp)
				{
					wrds.push_back(tmp);
				}
			} else {
				DEBUG(0,"unable to read parameter \"" << line << "\" in parameter-file!");
			}
		}
	}
	filestr.close();

	int argc_test = wrds.size();
//	cerr << __FUNCTION__ << argc_test << endl;


	char **argv_test = new char*[argc_test];
	for(unsigned int i = 0; i < wrds.size(); i++) {
		char *p = new char[wrds[i].size() + 1];
		strncpy(p, wrds[i].c_str(), wrds[i].size());
		p[wrds[i].size()] = '\0';
		argv_test[i] = p;
	}

	/**
	 * readArguments() is called twice, because in the first call the parameters form
	 * the file are read, in the second call the parameters from the command-line are added.
	 */

	readArguments(argc_test,argv_test);

//	cerr << "argc_saved " << argc_saved << endl;
	readArguments(argc_saved,argv_saved);

//	cerr << "yo" << endl;
}

std::string get_working_path()
{
	const int MAXPATHLEN = 1000;
	char temp[MAXPATHLEN];
	return ( getcwd(temp, MAXPATHLEN) ? std::string( temp ) : std::string("") );
}
