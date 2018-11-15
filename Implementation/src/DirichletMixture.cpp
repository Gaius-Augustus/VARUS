/*
 * DirichletMixture.cpp
 *
 *  Created on: 20.08.2016
 *      Author: willy
 */

#include "../headers/DirichletMixture.h"
#include <math.h>
#include "../headers/debug.h"
#include "../headers/Operators.h"
#include <unordered_map>

#include <fstream>
#include <random>
#include <stdio.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <assert.h>
#include "../headers/debug.h"
#include <time.h>
#include "../headers/Functions.h"
#include "../headers/Warning.h"

using namespace std;

DirichletMixture::DirichletMixture() : Estimator() {
	//Default
	param = nullptr;
	current_trainingsIteration = 0;
	iterationNumber = 0;
	sampler = nullptr;
}

DirichletMixture::DirichletMixture(ParameterHandler *p) : Estimator(p){
	/*! \brief Creates the components and initializes other values.
	 *
	 * The components are containers, into which the runs are sampled.
	 * Based on the samplings the parameters are then estimated.
	 */

	param = p;
	DEBUG(1,"Creating DM");
	for(int i = 0; i < param->components; i++) {
		DirichletDistribution c(param,1.0/(double)param->components);
		components.push_back(c);
	}
	current_trainingsIteration = 0;
	iterationNumber = 0;
	sampler = new myRandomEngine(p);

	DEBUG(1,"Done Creating DM");
}

DirichletMixture::~DirichletMixture() {
}

//double DirichletMixture::calculateLikelihood(DirichletDistribution &component, Run *r) {
//	/*! \brief Calculates the likelihood for a DirichletDistribution that it created the
//	 * data in an Experiment.
//	 *
//	 * Formula for calculation can be found in Altschull...
//	 */
//
////	double psC = 1.0;
////	double psCSum = psC*(double)param->numOfBlocks;
//	double pIk = 0.0;
//	UDmap::iterator j;
//
//	pIk = lgamma(component.alphaSum) - lgamma(component.alphaSum + r->observationSum);
//	for(j = component.alpha.begin(); j != component.alpha.end(); j++){
//    		pIk += lgamma(j->second + (double)r->observations[j->first])
//    										- (double) lgamma(j->second);
//    		DEBUG(4,"lgamma(" << j->second << ") = " << lgamma(j->second));
//	}
//    pIk = exp(pIk);
//
//    DEBUG(1,component.weight*pIk);
//	return component.weight*pIk;
//}

double DirichletMixture::logLikelihood(DirichletDistribution &comp) {
	return logLikelihood(comp, comp.alphaSum);
}

double DirichletMixture::logLikelihood(DirichletDistribution &comp, const double alphaSum) {
	/*! \brief Calculates the likelihood for a DirichletDistribution that it created the
	 * data in an Experiment.
	 *
	 * Formula (8) in Altschull-paper.
	 */

	double y = 0;
	double secondSum;
	for(int k = 0; k < comp.runs.size(); k++) {
		y += lgamma(alphaSum) - lgamma(alphaSum + comp.runs[k]->observationSum);

		DEBUG(4,lgamma(alphaSum) << ", " << lgamma(comp.runs[k]->observationSum));

		secondSum = 0;
		UDmap::iterator j;
		for(j = comp.alpha.begin(); j != comp.alpha.end(); j++){
			if(comp.runs[k]->observations[j->first] > 0) {
				DEBUG(4,alphaSum*comp.Qij[j->first]
								+ comp.runs[k]->observations[j->first]);

				secondSum += lgamma(alphaSum*comp.Qij[j->first]
						+ comp.runs[k]->observations[j->first])
						- lgamma(alphaSum*comp.Qij[j->first]);
			}
		}
		y += secondSum;
	}
	return y;
}

double DirichletMixture::dL(DirichletDistribution &comp, const double x) {
	/*! \brief The derivate of the likelihood needed for newtons-method.
	 *
	 * Formula (10) in Altschull-paper.
	 */

	double y = 0;
	double secondSum;
	for(int k = 0; k < comp.runs.size(); k++) {
		DEBUG(4,"x " << x);
		DEBUG(4,"comp.experiments[" << k << "].cCount: " << comp.runs[k]->observationSum);
		y += digamma(x) - digamma(x + comp.runs[k]->observationSum);

		secondSum = 0;
		UDmap::iterator j;
		for(j = comp.alpha.begin(); j != comp.alpha.end(); j++){
//			DEBUG(4,"digamma(" << x*comp.Qij[it->first] << ")");
//			DEBUG(4,"digamma(" << x*comp.dice[k].observations[it->first]<< ")");

			if(comp.Qij[j->first] > 0.) {
				secondSum = secondSum
							+ comp.Qij[j->first]*(digamma(x*comp.Qij[j->first]
													  + comp.runs[k]->observations[j->first])
													- digamma(x*comp.Qij[j->first]));
			}

		}
		y += secondSum;
	}
	return y;
}

double DirichletMixture::ddL(DirichletDistribution &comp, const double x) {
	/*! \brief The second derivate of the likelihood needed for newtons-method.
	 *
	 * Formula (11) in Altschull-paper.
	 */

	double y = 0;
	double secondSum;
	for(int k = 0; k < comp.runs.size(); k++) {
		y += trigamma(x) - trigamma(x + comp.runs[k]->observationSum);

		secondSum = 0;
		UDmap::iterator it;
		for(it = comp.alpha.begin(); it != comp.alpha.end(); it++) {
			if(comp.Qij[it->first] > 0.) {
				secondSum = secondSum
							+ comp.Qij[it->first]*comp.Qij[it->first]*(trigamma(x*comp.Qij[it->first]
																		   + comp.runs[k]->observations[it->first])
																		- trigamma(x*comp.Qij[it->first]));
			}
		}
		y += secondSum;
	}
	return y;
}

void DirichletMixture::exportLikelihoodCSV(unsigned int it, unsigned int compNum, DirichletDistribution &comp, const double a, const double b, const int n
						, vector<double> &newton) {
	/*! \brief Exports the likelihood as csv-file in the range from a to b with n equaly long steps.
	 *
	 *  A folder is created if necessary in which all the data for this component will go.
	 *  The csv-file has four columns "steps;L;dL;newtons". L is the likelihood, dL is its derivate and newtons
	 *  holds the values that have been produced by the newtons-method. The last value is the return-value of the newtons-method.
	 */


	double step = (double)(b-a)/(double)n;

	string command1 = "mkdir -p " + param->outFileNamePrefix + "/comp" + to_string(compNum);
	int status1 = system(command1.c_str());

	if (0 != status1) {
//		WARN("\nFailed to invoke:" << command1);
	}

	string command = "mkdir -p " + param->outFileNamePrefix + "/comp" + to_string(compNum) + "/alphaSum" + to_string(iterationNumber) + "/";
	int status = system(command.c_str());

	if (0 != status) {
//		WARN( "\nFailed to invoke:" << command);
	}


	string fileName = param->outFileNamePrefix + "/comp" + to_string(compNum) + "/alphaSum" + to_string(iterationNumber) + "/";
	fileName += "newton" + to_string(it);
	fileName += ".csv";

	fstream f;
	f.open(fileName, ios::out);
	if(!f) { DEBUG(0,"Cant open file " << fileName << "!");}

	f << "steps;L;dL;newtons" << "\n";
	unsigned int count = 0;
	for(double i = a; i < b; i += step) {
	    f << (double)i << ";" << logLikelihood(comp,i) << ";" << dL(comp,i);

	    if(count < newton.size()) {
	    	f << ";" << newton[count];
	    }

	    count++;
	    f << "\n";
	}
    f.close();

}

void DirichletMixture::exportRunNames() {
	/*! \brief Exports the names of the runs sampled into the bins as txt-files.
	 *
	 * Run-names are written to file line-wise.
	 * */


	for(unsigned int m = 0; m < components.size(); m++) {
		DEBUG(2,current_trainingsIteration << ":components[" << m << "]=" <<  components[m].runs.size());
		string command1 = "mkdir -p " + param->outFileNamePrefix + "/comp" + to_string(m);
		int status1 = system(command1.c_str());

		if (0 != status1) {
//			WARN("\nFailed to invoke:" << command1);
		}

		string command = "mkdir -p " + param->outFileNamePrefix + "/comp" + to_string(m) + "/alphaSum" + to_string(iterationNumber) + "/";
		int status = system(command.c_str());

		if (0 != status) {
//			WARN("\nFailed to invoke:" << command);
		}

		// save dice in this component
		string fileName2 = param->outFileNamePrefix + "/comp" + to_string(m) + "/alphaSum" + to_string(iterationNumber) + "/";
	//	fileName += toStr(comp.Qij);
		fileName2 += "runs" + to_string(current_trainingsIteration);
		fileName2 += ".csv";

		fstream f2;
		f2.open(fileName2, ios::out);
		if(!f2) { DEBUG(0,"Cant open file " << fileName2 << "!");}

		f2 << "runs" << "\n";

		for(unsigned int k = 0; k < components[m].runs.size(); k++) {
			f2 << components[m].runs[k]->accesionId << "\n";
		}

		f2.close();
	}
}

double DirichletMixture::newtonsMethod(unsigned int compNum, DirichletDistribution &c, const double start, const unsigned int maxIterations, const double tol) {

	/*! \brief Newtons method calculates the expected alphasum.
	 *
	 *
	 * The likelihood depending on alphasum is maximized. This is done by calculating the arguments at which
	 * the first derivate is zero.
	 */

	double x_k;
	if(start <= 0) {
//		WARN("Warning start <= 0. Setting x_k = 1.0");
		x_k = 0.01;
//	} else if(start < 0){
//		WARN("Warning start < 0. Setting x_k = -start");
//		x_k = -start;
//	}
	}else{
		x_k = start;
	}

	double x_o = x_k;

	vector<double> values;
	values.push_back(x_k);

	DEBUG(4,"newtonsMethod(" << start << ", " << maxIterations << ", " << tol << ")");



	for(unsigned int i = 0; i < maxIterations; i++) {
		x_o = x_k;
		DEBUG(4,"x[" << i << "] = " << x_k);
//		x_k = x_k - F1(x_k)/dF1(x_k);
		if(ddL(c, x_k) == 0){
			DEBUG(4,"ddL(c, x_k) == 0");
			DEBUG(4,"ddL: " << ddL(c,x_k));
			break;
		}

		x_k = x_k - dL(c, x_k)/ddL(c, x_k);
		values.push_back(x_k);

		if(fabs((double)(x_k - x_o)) < tol) { break;}
		if(x_k <= 0) {
			DEBUG(4,x_k << " <= 0; set x_k := " << x_o*0.5);
			x_k = 0.5*x_o;
			values.push_back(x_k);
		}	// <----------------------------------
	}

	if(param->exportNewtons == 1){
		exportLikelihoodCSV(current_trainingsIteration,compNum,c,0.001,2*x_k,1000,values);
	}
	DEBUG(4,"ddL: " << ddL(c,x_k));
	return x_k;
}


void DirichletMixture::calculateMixtureParameters(const unsigned int num) {
	/*! \brief Calculates the parameters according to the runs sampled in the bins.
	 *
	 */

	for(int i = 0; i < components.size(); i++) {
		components[i].calculateMixtureParameters(num);

		DEBUG(4,"EalphaSum(component[" << i << "]): " << components[i].alphaSum);

		components[i].alphaSum = newtonsMethod(i,components[i],fabs(components[i].alphaSum),param->newtonIterations,param->newtonPrecision);

		components[i].alpha = components[i].Qij;

		scaleMap(components[i].alpha,components[i].alphaSum);

		DEBUG(4, components[i].alpha );
		for(unsigned int k = 0; k < components[i].runs.size(); k++){
			DEBUG(4,components[i].runs[k]->accesionId);
		}
	}
}


void DirichletMixture::calculateLikelihoods(Run* r, vector<double> &likelihoods){
	for(unsigned int m = 0; m < components.size(); m++) {
		double l = components[m].calculateLikelihood(r);
//		double l = components[m].calculateDistance(r);
		likelihoods.push_back(l);
		DEBUG(4,"likelihood["<< m << "]:" << l);
	}

	DEBUG(4,sum(likelihoods));
	if(sum(likelihoods) == 0){
		for(unsigned int m = 0; m < components.size(); m++){
			likelihoods[m] = 1.0;
		}
	}
}

bool DirichletMixture::checkComponentsValid(){
	bool ret = true;
	for(unsigned int m = 0; m < components.size(); m++) {
		if(components[m].runs.size() < 2)
		{
//			WARN("Warning at least one component has fewer than 2 runs!");
			ret = false;
			break;
		}
	}
	return ret;
}

void DirichletMixture::train(std::vector<Run*> &runs) {
	/*! \brief The components are trained with gibbs-samplin.
	 *
	 */

	// resets the components, so that previous trainings will not affect this training
	for(unsigned int m = 0; m < components.size(); m++) {
		components[m].resetComplete(components.size());
	}

	// training
	for(current_trainingsIteration = 0; current_trainingsIteration < param->trainingsIterations; current_trainingsIteration++) {
		// reset everything except for the parameters q,alpha,alphasum
		for(unsigned int m = 0; m < components.size(); m++) {
			components[m].reset();
		}

		DEBUG(4,"_____b)_____");
		// if there is only one component all runs go in the same bin, else calculate as follows
		if(components.size() > 1) {
			// Sampling
			for(unsigned int k = 0; k < runs.size(); k++){

				vector<double> likelihoods;
				calculateLikelihoods(runs[k],likelihoods);

				DEBUG(4,"Sample run based on likelihoods");
				// samples the k-th dice into the component corresponding to the number generated by the discrete distribution
				components[sampler->select_randomly_p(likelihoods)].runs.push_back(runs[k]);
			}

			// throws a warning if a component has less than 2 runs sampled
			checkComponentsValid();

		} else {
			// if there is only one component all runs go in the same bin
			for(unsigned int k = 0; k < runs.size(); k++) {
				components[0].runs.push_back(runs[k]);
			}
		}

		//	c)
		DEBUG(4,"_____c)_____");
		calculateMixtureParameters(runs.size());

		exportMixture(current_trainingsIteration, translate2Str);
		exportRunNames();
	}
}

void DirichletMixture::estimateP(std::vector<Run*> &runs, unsigned int IterationNumber) {
	/*! \brief Estimates the p for all Experiments.
	 *
	 */

	this->iterationNumber = IterationNumber;

	DEBUG(1,"Training");
	train(runs);

	DEBUG(1,"calculateP");
	if(param->simpleDM == 1){
		calculateSimpleP();
	} else{
		calculateP(runs);
	}
}

void DirichletMixture::calculateP(std::vector<Run*> &runs) {
	/*! \brief calculates the estimators p for all runs after training.
	 *
	 */

	vector<double> logBComponents;
	// calculate denominator
	for(unsigned int u = 0; u < components.size(); u++) {
		logBComponents.push_back(myLogBeta(components[u].alpha));
	}

	UDmap::iterator it;
	for(unsigned int k = 0; k < runs.size(); k++) {

		vector<double> BComponents;
		for(unsigned int u = 0; u < components.size(); u++) {
			UDmap m;
//			m = components[u].alpha;
//			m += runs[k]->observations;
//			BComponents.push_back(myLogBeta(m) - logBComponents[u]);
			BComponents.push_back(myLogBeta(components[u].alpha, runs[k]->observations) - logBComponents[u]);
		}

		// calculate max
		double max = BComponents[0];
		for(unsigned int u = 1; u < components.size(); u++) {
			double t = BComponents[u];
			if(t > max) { max = t;}
		}

		for(unsigned int u = 0; u < components.size(); u++) {
			BComponents[u] = exp(BComponents[u] - max);
		}

		double xSum = 0.0;
		UDmap::iterator it;
		for(it=components[0].alpha.begin(); it!=components[0].alpha.end(); it++) {
			DEBUG(4,it->first);

			xSum += XiLog(BComponents, runs[k],it->first);

			DEBUG(4,XiLog(BComponents, runs[k],it->first));

			DEBUG(4,lgammaHash.size());
		}

		assert(xSum != 0);
		assert(components[0].alpha.size()>0);
		for(it=components[0].alpha.begin(); it!=components[0].alpha.end(); it++) {
			DEBUG(4,it->first);
			runs[k]->p[it->first] = (double)XiLog(BComponents, runs[k],it->first)/(double)xSum;
		}
	}
}

void DirichletMixture::initializeP(std::vector<Run*> &experiments) {
	// Do nothing.
}

//double DirichletMixture::myBeta(unordered_map<unsigned int, double> &v) {
//	/*! \brief Calculates the value of the beta-function.
//	 *
//	 */
//
//
//	double denominator = tgamma((float)sum(v));
//	DEBUG(1,sum(v));
//	DEBUG(1,denominator);
//	double ret = 1.0;
//	unordered_map<unsigned int, double>::iterator it;
//	for(it=v.begin(); it!=v.end(); it++) {
//		DEBUG(1,it->second);
//		ret = ret*tgamma(it->second)/denominator;
//	}
//	return ret;
//}

double DirichletMixture::myLogBeta(unordered_map<unsigned int, double> &v) {
	/*! \brief Calculates the logarithm of the beta-function of the vector v.
	 *
	 */

	double ret = 0.0;
	unordered_map<unsigned int, double>::iterator it;
	for(it=v.begin(); it!=v.end(); it++) {
//		if(lgammaHash.find(it->second) == lgammaHash.end()) {
//			lgammaHash[it->second] = lgamma(it->second);
//		}
//
//		ret += lgammaHash[it->second];

		ret += lgamma(it->second);
	}
	double su = sum(v);

//	if(lgammaHash.find(su) == lgammaHash.end()) {
//		lgammaHash[su] = lgamma(su);
//	}
//	ret = ret - lgammaHash[su];

	ret = ret - lgamma(su);

	return ret;
}

double DirichletMixture::myLogBeta(UDmap &v, UUmap &v2) {
	/*! \brief Calculates the logarithm of the beta-function of the vector v.
	 *
	 */

	double ret = 0.0;
	UDmap::iterator it;
	for(it=v.begin(); it!=v.end(); it++) {
//		if(lgammaHash.find(it->second) == lgammaHash.end()) {
//			lgammaHash[it->second] = lgamma(it->second);
//		}
//
//		ret += lgammaHash[it->second];

		ret += lgamma(it->second + v2[it->first]);
	}
	double su = sum(v);

//	if(lgammaHash.find(su) == lgammaHash.end()) {
//		lgammaHash[su] = lgamma(su);
//	}
//	ret = ret - lgammaHash[su];

	ret = ret - lgamma(su);

	return ret;
}

//double DirichletMixture::Xi(Experiment *e, const unsigned int i) {
//	/*! \brief Karplus formula (39).
//	 *
//	 */
//
//
//	double ret = 0.0;
//	for(unsigned int u = 0; u < components.size(); u++) {
//
//		// m = nVec + alpha_jVec
//		unordered_map<unsigned int, double> m;
//		m = components[u].getAlpha();
//		m += e->getObservations();
//
//		ret += components[u].getWeight()
//				*(myBeta(m)
//						/myBeta(components[u].getAlpha()))
//				*((e->getObservation(i) + components[u].getAlpha(i))
//						/(e->getObservationSum() + components[u].getAlphaSum()));
//	}
//	return ret;
//}

double DirichletMixture::XiLog(std::vector<double> &BComponents, Run *e, const U32 i) {
	/*! \brief Karplus formula (39)
	 *
	 */

	double ret = 0.0;
	for(unsigned int u = 0; u < components.size(); u++) {
		ret += components[u].weight
				*BComponents[u]
				*((e->observations[i] + components[u].alpha[i])
						/(e->observationSum + components[u].alphaSum));
	}

	return ret;
}

void DirichletMixture::exportMixture(const unsigned int iterationNumber,
										USmap &translate2Str) {
	/*! \brief Exports the parameters of the mixture as a csv-file.
	 *
	 */

	string fileName = param->outFileNamePrefix;
	fileName += "Dm" + std::to_string(iterationNumber);
	fileName += ".csv";


	fstream f;
	f.open(fileName, ios::out);
	if(!f) { DEBUG(0,"Cant open file " << fileName << "!");}


	DEBUG(4,components.size());
	// first line
	// out of for-loop because of ";"
	f << "blockName;"
	  << "c" << 0 << ".weight;"
	  << "c" << 0 << ".alpha;"
	  << "c" << 0 << ".alphaSum";

	for(unsigned int u = 1; u < components.size(); u++) {
		f << ";c" << u << ".weight;"
		  << "c" << u << ".alpha;"
		  << "c" << u << ".alphaSum";
	}

	f << "\n";
	// first line end

	UDmap::iterator j;
	for(j=components[0].alpha.begin(); j != components[0].alpha.end(); j++) {

		f << translate2Str[j->first] << ";";
		f << components[0].weight << ";"
			<< components[0].alpha[j->first] << ";"
			<< components[0].alphaSum;

		for(unsigned int k = 1; k < components.size(); k++) {
			f << ";" << components[k].weight
				<< ";" << components[k].alpha[j->first]
				<< ";" << components[k].alphaSum;
		}
		f << "\n";
	}
	f.close();
}

void DirichletMixture::setTranslate2Str(USmap &tr){
	// only needed here for exportMixture
	translate2Str = tr;
}

void DirichletMixture::calculateSimpleP(){
	for(unsigned int i = 0; i < components.size(); i++){
		components[i].calculateSimpleP();
	}
}

void DirichletMixture::estimatePWithNoObs(std::vector<Run*> &runs, unsigned int iterationNumber){

}
