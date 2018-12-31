/*
 * Controller.h
 *
 *  Created on: 30.11.2016
 *      Author: willy
 */

#ifndef CONTROLLER_H_
#define CONTROLLER_H_
#include "Run.h"
#include "TypeConventions.h"
#include "ChromosomeInitializer.h"
#include "Aligner.h"
#include "STAR_Aligner.h"
#include "HISAT_Aligner.h"
#include "Downloader.h"
#include "ParameterHandler.h"
#include "debug.h"
#include "Operators.h"
#include "myRandomEngine.h"
#include "Estimator.h"
#include "SimpleEstimator.h"
#include "AdvancedEstimator.h"
#include "DirichletMixture.h"
#include "Simulator.h"
#include "ClusterEstimator.h"
#include <memory>

class Controller {
    /*! \brief Controlls the other objects.
     *
     * The main-algorithm is implemented in algorithm().
     */

 public:
    ParameterHandler *param;

    ChromosomeInitializer *chrom;
    Aligner *align;
    Downloader *down;
    Simulator *sim;
    myRandomEngine *ran;
    Estimator *est;

    unsigned int batchCount;
    unsigned int goodBatchCount;
    unsigned int lastMerge;

    double maxProfit;

    double totalProfit;		// profit = gain - costs

    std::vector<Run*> runs;			// all runs
    std::vector<Run*> downloadableRuns;		// runs that still have reads left to be downloaded
    std::vector<Run*> badQualityRuns;		// runs that at some point were found to have bad quality during
    // execution of the program

    UUmap totalObservations;
    double totalScore;
    double avgUniq;
    double avgSpliced;
    
    Controller(ParameterHandler *p);
    virtual ~Controller();

    void initialize();

    void algorithm();

    Run* chooseNextRun();

    bool continuing();

    double score(const double v);

    void profit(Run *d);

    void calculateProfit(std::vector<Run*> &runs);

    double score(UUmap &obs);

    void exportTotalObservationCSV(std::string name);

    void exportTotalObservationCSVlessInfo(Run *r); //std::string name);

    void exportCoverage();

    void printRun2File(Run *r);

    void createDieFromRun(Run *r);

    void createDiceFromRuns();

    void updateDownloadableRuns();

    void mergeAlignments(const int count);

    void finalMerge();

    void exportRunStatistics();
};

#endif /* CONTROLLER_H_ */
