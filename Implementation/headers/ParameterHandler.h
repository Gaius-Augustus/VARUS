/*
 * ParameterHandler.h
 *
 *  Created on: 25.11.2016
 *      Author: willy
 */

#ifndef PARAMETERHANDLER_H_
#define PARAMETERHANDLER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include "TypeConventions.h"
#include <memory>
#include <set>


/*! \brief This class holds all the data, that is read in as command-line arguments.
 *
 * A raw pointer is passed to all other instances, that need to have access to the input data.
 * A std::unique_ptr is created in main.cpp that is responsible for this instance.
 */

/*
 * Holds the description and the category for a given parameter.
 * E.g. categories could be the level of details for each model.
 * More general information should be on a higher level in the manual and
 * should be easier accessible.
 */

struct Parameter{
public:
    std::string name;
    std::string value;
    std::string description;
    std::string category;

    bool deprecated;
};

class ParameterHandler
{


 public:
    ////     paths, commandline
    std::string pathToSTAR;
    std::string pathToRuns;
    //    std::string readFilesIn;
    std::string genomeDir;
    std::string outFileNamePrefix;
    uI runThreadN;
    uI blockSize;
    uI batchSize;
    double pseudoCount;
    double lambda;
    int logScore;
    double cost;
    int createDice;
    int simulation;
    int estimator;
    //	int dieList;
    int lessInfo = 1;
    uI components;
    uI numOfBlocks;
    uI trainingsIterations;
    int loadAllOnce;
    int verbosity;
    int verbosityDebug;
    int randomSeed;
    int profitCondition;

    std::string fastqDumpCall;


    uI newtonIterations;
    double newtonPrecision;

    int exportObservationsToFile;

    int readParametersFromFile;
    std::string pathToParameters;

    //		bool readAllready;
    int exportParametersToFile;

    std::string pathToDice;
    int deleteLater;

    unsigned int maxBatches;

    int ignoreReadNum;

    int simpleDM;

    //	int kMeansIterations;

    int exportNewtons;


    //	std::shared_ptr<Parameters> parameters;


    bool readAllready;

    int mergeThreshold;

    double qualityThreshold;

    std::string currentWorkingDirectory;

    std::string pathToVARUS;

    std::string varusID;

    // for reading parameters from file and afterwards add commandline arguments
    int argc_saved;

    char **argv_saved;

    //	std::map <std::string, std::string> parameters;
    //	std::map<std::string, std::string> usage;


    std::map <std::string, Parameter> parameters;

    enum paramCat{
	QUICKSTART,
	BASICSETTINGS,
	DIRICHLETMIXTURE,
	SIMULATION,
	SIMPLEANDADVANCED,
	CLUSTERESTIMATOR,
	MANDATORY,
	ESTIMATOR,
	TEST,
	COMMANDLINE
    };

    std::string STARmanual;

    std::map<paramCat, std::string> parameterCategories;		//key is name of the category in code, value is name of the category for the usage

    unsigned mincovthresh = 3;  // minimum number of splices to make the db

    void printParameterCategory(const paramCat cat, const std::string des, const unsigned int maxParLength);

    void read_parameters_from_file(std::string path);	// uses parameters stored in the file

    void export_parameters();	// exports the parameters into a file, so that they can later be used

    void readArguments(int argc, char *argv[]);	// reads the arguments from the commandline
    void readInputRuns();	// reads in the list of runs that should be downloaded from
    std::ostream & printParameters(std::ostream& os);	// prints the parameters

    void print_usage();

    const unsigned int lineWidth = 80;

    std::string printTextWithWidthAndLeftOffset(const std::string s, const unsigned int l, const unsigned int off);

    std::string lineLength(const std::string s, const unsigned int l, const unsigned int maxS);


    void exit_text();

    template<typename K>
	void add_parameter(std::string name, K a, std::string category = "", std::string description = "", bool deprecated = false);

    ParameterHandler();
    ~ParameterHandler();

 private:
    ParameterHandler *param;
};

std::string get_working_path();


#endif /* PARAMETERHANDLER_H_ */


