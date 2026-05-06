/*
 * ClusterEstimator.cpp
 *
 *  Created on: 29.12.2016
 *      Author: willy
 */

#include "../headers/ClusterEstimator.h"
#include <assert.h>
#include "../headers/Operators.h"
#include <fstream>

ClusterEstimator::ClusterEstimator(ParameterHandler *p) : Estimator(p) {
	/*! \brief p->components ClusterComponents are created.
	 */

	DEBUG(1,"Creating cluster-components");
	for(unsigned int i = 0; i < param->components; i++){
		ClusterComponent c(p);
		components.push_back(c);
	}

	ran = new myRandomEngine(p);

}

ClusterEstimator::~ClusterEstimator() {
	// TODO Auto-generated destructor stub
}

void ClusterEstimator::initializeClusters(std::vector<Run*> &runs){
	/* Sets the start-points for the clusters.
	 *
	 * Each component is assigned a Runs observation. These observations are used to
	 * build the starting point.
	 *
	 */

	DEBUG(1,"Initializing clusters");
	std::vector<Run*> startPoints = runs;

	int num = startPoints.size();

	for(unsigned int i = 0; i < components.size() && i < num; i++){
		int ind = ran->select_randomly(startPoints.size()-1);
		components[i].setSeed(startPoints[ind]->q);
		startPoints.erase(startPoints.begin() + ind);
	}

	if(components.size() > runs.size()){
		for(unsigned int i = runs.size(); i < components.size(); i++){
			components[i].setSeed(runs[0]->q);
		}
	}
}

void ClusterEstimator::calculateRunQs(std::vector<Run*> &runs){
	/* Calculates the q map of all runs.
	 *
	 * The q-map holds the frequencys of all observations. It is needed for the clustering-step.
	 */

	DEBUG(1,"calculating Run qs");
	for(unsigned int k = 0; k < runs.size(); k++){
		UUmap::iterator j;
		for(j = runs[k]->observations.begin(); j != runs[k]->observations.end(); j++){
			runs[k]->q[j->first] = (double)j->second/(double)runs[k]->observationSum;
		}
	}
}

void ClusterEstimator::kMeans(std::vector<Run*> &runs){
	/*! \brief kMeans is a clustering algorithm.
	 *
	 * Each point is a frequency of observations. Among all frequencies the best cluster given n clusterpoints is
	 * searched.
	 */

	DEBUG(1,"k-means");
	for(unsigned int l = 0; l < param->trainingsIterations; l++){
		for(unsigned int i = 0; i < components.size(); i++){
			components[i].releaseRuns();
		}

		for(unsigned int k = 0; k < runs.size(); k++){
			double minDist = distance(runs[k],components[0]);
			int compInd = 0;

			for(unsigned int i = 1; i < components.size(); i++){
				double dist = distance(runs[k], components[i]);
				DEBUG(4,"distance to cluster " << i << ";" << dist);
				if(minDist > dist){
					minDist = dist;
					compInd = i;
				}
			}
			components[compInd].runs.push_back(runs[k]);
		}

		for(unsigned int i = 0; i < components.size(); i++){
			components[i].updateCenter();
		}
	}
}

double ClusterEstimator::distance(Run *r, ClusterComponent &c){
	/* The Kullback-Leibler metric.
	 *
	 * \begin{equation}
	 * 	\text{dist}(\text{run},\text{cluster}) := \sum_{j=1}^T \log(\frac{q_j}{\alpha_j}) \cdot q_j
	 * \end{equation}
	 *
	 */

	double dist = 0.0;
	UDmap::iterator j;
	for(j = c.alpha.begin(); j != c.alpha.end(); j++){
		dist += log(r->q[j->first]/j->second)*r->q[j->first];
	}
	dist = pow(dist,0.5);
	return dist;
}

bool pred(Run* s) {
	/*! \brief Checks if the observationSum of a given run is equal to 0.
	 */

	if(0 == s->observationSum) return true;
	else return false;
}

void ClusterEstimator::estimateP(std::vector<Run*> &runs, unsigned int iterationNumber){
	/*! \brief Calculates the estimation for the runs by first clustering the runs based on the
	 * similarity of the observations.
	 */


	std::vector<Run*> runsCopy = runs;

	// calculate p for runs with 0 observations
	calculatePnoObs(runsCopy);

	// removes all runs with 0 observations
	runsCopy.erase( std::remove_if( runsCopy.begin(), runsCopy.end(), pred ), runsCopy.end() );

	calculateRunQs(runsCopy);

	initializeClusters(runsCopy);

	kMeans(runsCopy);

	exportClusters(iterationNumber);

	for(unsigned int i = 0; i < components.size(); i++){
		components[i].calculateAlphaSum();
		components[i].calculateP();
	}
}

void ClusterEstimator::exportClusters(unsigned int iterationNumber){
	/*! \brief Exports the cluster into a file.
	 *
	 * After a number indicating the cluster-number runs are written.
	 * All runs written below a number are assigned to the same cluster.
	 */

	std::string fileName = param->outFileNamePrefix;
	fileName += "/ClusterEstimator" + std::to_string(iterationNumber);
	fileName += ".csv";


	std::fstream f;
	f.open(fileName, std::ios::out);
	if(!f) { DEBUG(0,"Cant open file " << fileName << "!");}

	for(unsigned int i = 0; i < components.size(); i++){
		f << "component " << i << "\n";
		for(unsigned int k = 0; k < components[i].runs.size(); k++){
			f << components[i].runs[k]->accesionId << "\n";
		}
		f << "----------------------\n";
	}
}

void ClusterEstimator::calculatePnoObs(std::vector<Run*> &runs){
	/*! \brief Calculates the estimation for runs with no observations.
	 *
	 * The estimation is done by using the observations from all runs with observations.
	 */

	UUmap obs_total;	// c_j
	unsigned int observationSum = 0;							// c^*
	// calculate observations_total and sum(observations)
	for(unsigned int k = 0; k < runs.size(); k++) {
		if(runs[k]->observationSum > 0){
			observationSum 	+= runs[k]->observationSum;
			obs_total 		+= runs[k]->observations;
		}
	}

	UDmap p_total;
	UUmap::iterator j;
	for(j=obs_total.begin(); j!=obs_total.end(); j++) {
		// p_j = c_j/c*
		p_total[j->first] = (double)obs_total[j->first]/(double)observationSum;
	}

	double a_T = param->pseudoCount*param->numOfBlocks;
	double lambda_T = param->lambda*param->numOfBlocks;

	for(unsigned int k = 0; k < runs.size(); k++) {
		if(0 == runs[k]->observationSum){
			// c^k_{n_k}* + a/(c*+1)T + lambda
			double normalizingConst = a_T + lambda_T;

			UUmap::iterator j;
			for(j=runs[k]->observations.begin(); j!=runs[k]->observations.end(); j++) {
				// p_j = c_j + a/(c*+1) + lambda*p_total
				runs[k]->p[j->first] =
						((double)param->pseudoCount
						+ lambda_T*p_total[j->first])
						/normalizingConst;
			}
		}
	}

}

void ClusterEstimator::estimatePWithNoObs(std::vector<Run*> &runs, unsigned int iterationNumber){

}



