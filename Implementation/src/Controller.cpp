/*
 * Controller.cpp
 *
 *  Created on: 30.11.2016
 *      Author: willy
 */

#include "../headers/Controller.h"
#include <random>
#include "../headers/Warning.h"
#include "../headers/debug.h"
#include <assert.h>

using namespace std;

Controller::Controller(ParameterHandler *p) {
    /*! \brief Different Objects are created based on the settings in the passed ParameterHandler p.
     *
     */

    param 	= p;
    DEBUG(1,"Creating Controller");
    
    if (param->simulation != 1){
	chrom 	= new ChromosomeInitializer(p);
	down 	= new Downloader(p);
	allign 	= new Alligner(p);
	sim 	= nullptr;
    } else {
	sim 	= new Simulator(p);
	chrom 	= nullptr;
	down 	= nullptr;
	allign 	= nullptr;
    }


    ran	= new myRandomEngine(p);
    
    if(p->estimator == 1){
	est = new SimpleEstimator(p);
    } else if(p->estimator == 2){
	est = new AdvancedEstimator(p);
    } else if(p->estimator == 3){
	est = new DirichletMixture(p);
    } else if(p->estimator == 4){
	est = new ClusterEstimator(p);
    } else {
	est = nullptr;
    }

    batchCount = 0;
    lastMerge = 0;
    maxProfit = 1;
    totalScore = 0;
    goodBatchCount = 0;
    
    totalProfit = 0.0;
    
    DEBUG(1,"Done Creating Controller");
}

Controller::~Controller() {
}

void Controller::initialize(){
    /*! \brief Initialization of the chromosomeinitializer.
     *  Also sets the seed for the random-engines.
     */
    
    DEBUG(1,"Initializing Controller");
    if (param->simulation != 1){
	chrom->readInputRuns();
	chrom->getChromosomeLengths();
	chrom->create_transcript_units(totalObservations);
	chrom->initializeRuns(runs);
    } else {
	DEBUG(1,"Initializing simulator")
	    sim->readInputRuns();
	sim->initializeRuns(runs);
	initMap(totalObservations,(unsigned int) 0, param->numOfBlocks);
	DEBUG(1,"Done initializing simulator");
    }

    downloadableRuns = runs;

    if(param->estimator != 0){
	est->initializeRuns(runs,totalObservations);
	if(param->estimator == 3){
	    DirichletMixture *d = (DirichletMixture*) est;
	    if(param->simulation == 1){
		d->setTranslate2Str(sim->translate2str);
	    } else {
		d->setTranslate2Str(chrom->translate2str);
	    }
	    DEBUG(1,"Done initializing translate2str for DM");
	}
    }
    DEBUG(1,"Done Initializing in Controller");
}

void Controller::algorithm(){
	/*! \brief The main part of the program.
	 *
	 * The program will stay in the while-loop until the expected score drops below 0.
	 */

	maxProfit = 1;

	if (param->loadAllOnce == 1){
	    DEBUG(1,"Loading all once");
	    for(unsigned int k = 0; k < runs.size(); k++){
		if(param->simulation != 1){
		    down->getBatch(runs[k]);
		    allign->update(runs[k],totalObservations,chrom);
		} else{
		    sim->simulateObservations(runs[k],totalObservations);
		}
	    }
	    
	    if(param->estimator != 0){
		est->estimateP(runs, 0);
	    }
	    // calculates maxProfit
	    calculateProfit(downloadableRuns);
	    Run *tmp_run = chooseNextRun();
	    DEBUG(1,"Done Loading all once");
	}
	
	updateDownloadableRuns();

	while(continuing()){
	    batchCount++;

	    // if maxbatches < batchCount stop downloading
	    if(!continuing()) break;




	    calculateProfit(downloadableRuns);

	    DEBUG(1,"Choosing next Run ...");
	    Run *tmp_run = chooseNextRun();
	    DEBUG(1,"...done Choosing next Run");

	    // maxProfit set in chooseNextRun()
	    DEBUG(1,"Iteration: " << batchCount << "(" << goodBatchCount << " good batches) ("<< runs.size() - badQualityRuns.size() << " good runs of " << runs.size() << ") | maxProfit: " << maxProfit);

	    // if profit < 0 stop downloading
	    if(!continuing()) break;

	    if(param->simulation != 1){
		assert(tmp_run != nullptr);
		DEBUG(1,"Getting batch number " << tmp_run->timesDownloaded + 1 << " for " << tmp_run->accesionId);
		down->getBatch(tmp_run);
		allign->update(tmp_run,totalObservations,chrom);
		if(tmp_run->badQuality == false){
		    goodBatchCount++;
		}
	    } else{
		sim->simulateObservations(tmp_run,totalObservations);
	    }

	    totalScore = score(totalObservations);
	    totalProfit = totalScore - param->cost*param->batchSize*batchCount;

	    if(param->estimator != 0){
		DEBUG(1,"Estimating...")
		    // we pass all runs since the DM uses the information from all runs
		    est->estimateP(runs, batchCount);
		DEBUG(1,"... done Estimating")
		    }

	    calculateProfit(downloadableRuns);

	    DEBUG(1,"Exporting...")
		if(1 == param->lessInfo){
		    exportTotalObservationCSVlessInfo(tmp_run->accesionId);
		}else{
		    exportTotalObservationCSV(tmp_run->accesionId);
		}

	    exportRunStatistics();

	    DEBUG(1,"...done Exporting")

		if(param->ignoreReadNum != 1){
		    updateDownloadableRuns();
		}

	    if(downloadableRuns.size() == 0){
		DEBUG(0,"No downloadable Runs left. Exiting...");
		exportCoverage();
		finalMerge();
		param->exit_text();
	    }

	    if(param->deleteLater == 1 && batchCount%param->mergeThreshold == 0){
		mergeAlignments(goodBatchCount);
	    }
	}

//	if(batchCount%param->mergeThreshold != 0){
	exportCoverage();
	finalMerge();

//	}
}

Run*  Controller::chooseNextRun(){
    /*! \brief This function returns a pointer to the Run that should be downloaded next.
     *
     * If no estimator is selected, the next run is choosen randomly. If two or more
     * Runs yield equal scores the Run is choosen randomly among the runs with the
     * highest score.
     */


    assert(downloadableRuns.size() > 0);
    unsigned int index = 0;
    vector<Run*> candidates;

    // in this case the next run is choosen randomly
    if(param->estimator == 0){
	// chooses a run randomly
	return ran->select_randomly(downloadableRuns);
    } else {
	if(downloadableRuns.size() > 0){
	    maxProfit = downloadableRuns[0]->expectedProfit;
	    candidates.push_back(downloadableRuns[0]);

	    for(unsigned int i = 1; i < downloadableRuns.size(); i++){
		// found a better score
		if(maxProfit < downloadableRuns[i]->expectedProfit){
		    candidates.clear();
		    candidates.push_back(downloadableRuns[i]);
		    maxProfit = downloadableRuns[i]->expectedProfit;
		} else if(maxProfit == downloadableRuns[i]->expectedProfit){ // score is equal
		    candidates.push_back(downloadableRuns[i]);
		}
	    }

	    return ran->select_randomly(candidates);
	}
    }
    return downloadableRuns[index];
}

bool Controller::continuing(){
	/*! \brief This function implements a criteria for exiting the program.
	 *
	 * Returns true if the program should be continued, false otherwise.
	 */

	if(param->maxBatches > 0 && batchCount > param->maxBatches){
		DEBUG(0,"maxBatches: " << param->maxBatches);
		return false;
	}

	if(1 == param->profitCondition && maxProfit <= 0){
		DEBUG(0,"maxProfit: " << maxProfit);
		return false;
	}

	return true;
}

double Controller::score(const double v) {
    /*! \brief Returns the contribution of a transcript.
     */
    return log(1.0 + v);
}

double Controller::score(UUmap &obs){
    /*! \brief Calculates the score for the observations of a run or the totalObservations.
     *
     */
    
    double sc = 0.0;
    UUmap::iterator i;
    for(i = obs.begin(); i != obs.end(); i++){
	sc += score(i->second);
    }
    return sc;
}

void Controller::profit(Run *d) {
	/*! \brief Calculates the expected profit that will be achieved by downloading from the run d*.
	 */

	long double pr = 0.0;
	long double pv = 0.0;

	DEBUG(4,"calculating profit" << d->p.size());

	UDmap::iterator it;
	for (it=d->p.begin(); it != d->p.end(); ++it) {
	    pr += score((double) it->second * (double) param->batchSize
			+ (double) totalObservations[it->first])
		- score((double) totalObservations[it->first]);
	    
	    DEBUG(4,"pr " << pr);
	    pv += it->second;
	}

	pr -= param->cost * (double)param->batchSize;

	// set the score of this Experiment
	d->expectedProfit = pr;
}

void Controller::calculateProfit(vector<Run*> &runs){
    /*! \brief Calculates the expected profit for all runs.
     */

    // calculate profit for runs with no reads only once
    double noReadsProfit = std::numeric_limits<int>::min();

    for(unsigned int i = 0; i < runs.size(); i++){
	if(runs[i]->timesDownloaded == 0 &&
	   noReadsProfit != std::numeric_limits<int>::min()){
	    runs[i]->expectedProfit = noReadsProfit;
	} else {
	    profit(runs[i]);
	    if (runs[i]->timesDownloaded == 0){
		noReadsProfit = runs[i]->expectedProfit;
	    }
	}
    }
}

void Controller::exportTotalObservationCSV(string name) {
    /*! \brief Exports the observations in csv-format.
     *
     * The passed string "name" is the accession-id of the run that was chosen to be downloaded from.
     * Also the observations and p-vectors of all runs are exported.
     */



    if(param->exportObservationsToFile != 0){
	string fileName = param->outFileNamePrefix;
	fileName += "/Toy" + std::to_string(batchCount);
	fileName += ".csv";
	
	fstream f;
	f.open(fileName, ios::out);
	if(!f) { DEBUG(0,"Cant open file " << fileName << "!");}

	// first line
	f << "blockName;totalObservations;totalScore;totalProfit;best";
	for(unsigned int k = 0; k < runs.size(); k++) {
	    f << ";" << runs[k]->accesionId << ".obs;"
	      << runs[k]->accesionId<< ".p;" << runs[k]->accesionId<< ".score";
	}
	f << endl;

	if(param->simulation != 1){
	    for(unsigned int j=0; j < totalObservations.size(); j++) {
		f << chrom->translate2str[j]<< ";" << totalObservations[j]
		  << ";" << totalScore << ";" << totalProfit << ";" << name;

		for(unsigned int k = 0; k < runs.size(); k++) {
		    f << ";" << runs[k]->observations[j] << ";"
		      << runs[k]->p[j] << ";" << runs[k]->expectedProfit;
		}
		f << endl;
	    }
	} else {
	    for(unsigned int j=0; j < totalObservations.size(); j++) {
		f << sim->translate2str[j]<< ";" << totalObservations[j]
		  << ";" << totalScore << ";" << totalProfit << ";" << name;

		for(unsigned int k = 0; k < runs.size(); k++) {
		    f << ";" << runs[k]->observations[j] << ";"
		      << runs[k]->p[j] << ";" << runs[k]->expectedProfit;
		}
		f << endl;
	    }
	}
	f.close();
    }
}

void Controller::exportTotalObservationCSVlessInfo(string name) {
    /*! \brief Exports the observations in csv-format with less information.
     *
     * Only the total observations and the best score and accession-id of the best run is exported.
     * In exportTotalObservationsCSV() also the observations and p-vectors of all runs are exported.
     */


    if(param->exportObservationsToFile != 0)
	{
	    string fileName = param->outFileNamePrefix;
	    fileName += "Coverage" + std::to_string(batchCount);
	    fileName += ".csv";

	    fstream f;
	    f.open(fileName, ios::out);
	    if(!f) { DEBUG(0,"Cant open file " << fileName << "!");}

	    // first line
	    f << "blockName;totalObservations;totalScore;totalProfit;best" << "\n";

	    if(param->simulation != 1){
		for(unsigned int j=0; j < totalObservations.size(); j++) {
		    f << chrom->translate2str[j]<< ";" << totalObservations[j] << ";" << log(totalScore) << ";" << totalProfit << ";" << name << "\n";

		}
	    }else{
		for(unsigned int j=0; j < totalObservations.size(); j++) {
		    f << sim->translate2str[j]<< ";" << totalObservations[j] << ";" << log(totalScore) << ";" << totalProfit << ";" << name << "\n";

		}
	    }
	    f.close();
	}
}

void Controller::exportCoverage(){
    string fileName = param->outFileNamePrefix;
    fileName += "Coverage";
    fileName += ".csv";

    fstream f;
    f.open(fileName, ios::out);
    if(!f) { DEBUG(0,"Cant open file " << fileName << "!");}

    // first line
    f << "blockName;totalObservations" << "\n";

    if(param->simulation != 1){
	for(unsigned int j=0; j < totalObservations.size(); j++) {
	    f << chrom->translate2str[j]<< ";" << totalObservations[j] << "\n";

	}
    }else{
	for(unsigned int j=0; j < totalObservations.size(); j++) {
	    f << sim->translate2str[j]<< ";" << totalObservations[j] << "\n";
	}
    }
    f.close();
}

void Controller::printRun2File(Run *r) {
    /*! \brief Exports the observations of the run into a csv-file.
     *
     * This function is mainly used for creating the dice from real runs.
     * These dice can later be used for simulations. In order to do that
     * all reads in the run are downloaded and then alligned. Based on
     * the observations it is possible to make probability distributions.
     */

    string csvName = param->pathToDice;
    csvName += r->accesionId;
    csvName += ".csv";

    fstream f2;
    f2.open(csvName, ios::out);
    if (!f2) {
	DEBUG(0, "Cant open file " << csvName << "!");
    }
    for (int i = 0; i < chrom->order.size(); i++) {
	f2 << chrom->translate2str[i] << ";" << r->observations[i] << "\n";
    }
    f2.close();
}

void Controller::createDieFromRun(Run *r) {
    /*! \brief Completely downloads the given run stepwise and alligns it.
     *
     * The resulting observations are then printed into a csv-file, which can later
     * be used for simulations.
     */

    int c = 1;
    while(r->maxNumOfBatches >= c) {
	DEBUG(0,"loading batch " << c << "/" << r->maxNumOfBatches);
	c++;
	down->getBatch(r);
	allign->update(r,totalObservations,chrom);
    }

    printRun2File(r);
}

void Controller::createDiceFromRuns() {
    /*! \brief Creates dice from all runs.
     *
     * It is first looked for a file "completedRUns.txt". In this file the runs that are
     * already created as a die are listed. Then a "Runlist.txt" is read in.
     * All runs except the runs being listed on the "completedRuns.txt"-file are then created as dice.
     */

    //	string com = ;
    string com = param->pathToDice + "completedRuns.txt";

    const char * c = com.c_str();
    string line;
    ifstream file(c);

    if (!file.is_open()) {
	DEBUG(0,"Cant open " << com << ". Will create new file completedRuns.txt");
	//        return;
    }

    // getting all the runs already downloaded out of the file "completedRuns.txt"
    // the names are saved line-wise
    unordered_map<string, int> completedRuns;
    while (getline(file, line)) {
	if (line.size() > 1) {
	    string tmp;
	    stringstream ssin(line);
	    ssin >> tmp;
	    completedRuns[tmp] = 1;
	}
    }
    DEBUG(3, completedRuns.size());
    DEBUG(3, completedRuns);

    // to change the order of the runs downloaded randomly
    ran->shuffle(runs);

    // iterating over the runlist, deleting runs that are downloaded already
    for(unsigned int i = 0; i < runs.size(); i++){
	if (completedRuns.find(runs[i]->accesionId) == completedRuns.end()) {

	    createDieFromRun(runs[i]);
	    completedRuns[runs[i]->accesionId] = 1;

	    // update number
	    fstream f;
	    f.open(com, std::ios_base::app);
	    if (!f) {
		DEBUG(0, "Cant open file " << com << "!");
	    }
	    f << runs[i]->accesionId << endl;
	    DEBUG(3,runs[i]->accesionId);
	    f.close();
	    runs.erase(runs.begin()+i);
	} else {
	    runs.erase(runs.begin()+i);
	}
    }
}

void Controller::updateDownloadableRuns(){
    vector< Run* >::iterator it = downloadableRuns.begin();
    
    while(it != downloadableRuns.end()) {
	if((*it)->badQuality == true){
	    DEBUG(0,"Deleting run " << (*it)->accesionId << " because of bad quality.");
	    badQualityRuns.push_back(*it);
	    it = downloadableRuns.erase(it);
	}
	else if((*it)->timesDownloaded*param->batchSize >= (*it)->numOfSpots){
	    DEBUG(0,"Deleting run " << (*it)->accesionId << " because no reads are left to download.");
	    it = downloadableRuns.erase(it);
	}

	else ++it;
    }
}

void Controller::mergeAlignments(const int count){
	/*
	 * Merges all alignments that are in subfolders of the current working directory.
	 * This is done in two steps:
	 * First the existing Aligned.out.sam -files are converted to bam-files using samtools.
	 * In a second step these resulting bam-files are merged using samtools.
	 *
	 * Uses the script /scripts/mergeAlignments.sh
	 */

	/*
	 * Merges certain bam-files that are merged alignments themselves.
	 * Example with mergeThreshold = 10:
	 *
	 * 					...						|	...
	 * 				/			\				|
	 * 		  (1_100)	...		(901_1000)		|	level 2
	 * 		/		\							|
	 * (1_10) ... (91_100)		/ ...	\		|	level 1
	 *
	 *	The required number of levels is log_(mergeThreshold)(batchCount).
	 *	That means over time the number of files is increasing.
	 */

	unsigned int requiredLevels = ceil(log(batchCount) / log(param->mergeThreshold));
	unsigned int levelOfCurrentFile = ceil(log(batchCount-lastMerge) / log(param->mergeThreshold));	// always 1, don't need to calculate here

	/*
	 * For each subfolder search for a sam-file called Aligned.out.sam and convert it to a file
	 * Aligned.out.bam.
	 */


	string s = param->pathToVARUS + "/scripts/mergeAlignments.sh merged" + std::to_string(levelOfCurrentFile)+ "." + std::to_string(lastMerge + 1) + "_" + std::to_string(batchCount) + ".bam " + std::to_string(param->mergeThreshold);

	if(param->deleteLater){
		s = s + " delete";
	}
	s = s + " >nul 2>&1";

	DEBUG(0, "Running '" << s <<"'");
	const char * c = s.c_str();
	int status = system(c);
	if(0 != status)
	{
		DEBUG(0,"Failed to run 'mergeAlignments.sh' properly!");
		param->exit_text();
	}

	lastMerge = batchCount;

}

void Controller::finalMerge(){
	/*
	 * Merges certain bam-files that are merged alignments themselves.
	 * Example with mergeThreshold = 10:
	 *
	 * 					...						|	...
	 * 				/			\				|
	 * 		  (1_100)	...		(901_1000)		|	level 2
	 * 		/		\							|
	 * (1_10) ... (91_100)		/ ...	\		|	level 1
	 *
	 *	The required number of levels is log_(mergeThreshold)(batchCount).
	 *	That means over time the number of files is increasing.
	 */

	unsigned int requiredLevels = ceil(log(batchCount) / log(param->mergeThreshold));
	string s = param->pathToVARUS + "/scripts/finalMerge.sh";

	if(param->deleteLater){
		s = s + " delete";
	}
	s = s + " >nul 2>&1";

	string sCall = "Running '" + s + "'";
	DEBUG(0, sCall);
	const char * c = s.c_str();
	int status = system(c);
	if(0 != status)
	{
		DEBUG(0,"Failed to run 'finalMerge.sh' properly!");
		param->exit_text();
	}
}

bool runCount (Run *r1, Run *r2) {
    if (r1->timesDownloaded > r2->timesDownloaded){
	return true;
    } else if (r1->timesDownloaded == r2->timesDownloaded){
	return (r1->avgUmrPercent > r2->avgUmrPercent);
    }

    return false;
}

void Controller::exportRunStatistics(){
    string fileName = param->outFileNamePrefix;
    fileName += "RunStatistics";
    fileName += ".csv";

    std::sort (runs.begin(), runs.end(), runCount);

    fstream f;
    f.open(fileName, ios::out);
    if(!f) { DEBUG(0,"Cant open file " << fileName << "!");}

    // first line
    f << "accesion-id;timesDownloaded;avgUmrPercent;badQuality;maxNumOfBatches" << "\n";

    for(unsigned int j=0; j < runs.size(); j++) {
	f << runs[j]->accesionId<< ";" << runs[j]->timesDownloaded<< ";" <<runs[j]->avgUmrPercent << "%" << ";" << runs[j]->badQuality << ";" << runs[j]->maxNumOfBatches << "\n";

    }
    f.close();
}
