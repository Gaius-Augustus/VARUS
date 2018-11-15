/*
 * DirichletDistribution.cpp
 *
 *  Created on: 20.08.2016
 *      Author: willy
 */

#include "../headers/DirichletDistribution.h"
#include "../headers/Operators.h"
#include "../headers/debug.h"
#include "../headers/Functions.h"
#include "../headers/Warning.h"

using namespace std;

DirichletDistribution::DirichletDistribution(ParameterHandler *p, const double w) {

	param = p;
	DEBUG(1,"Creating DirichletDistribution");
	size = param->numOfBlocks;
	weight = w;
	initMap(alpha,1.0,size);
	initMap(aggregateCCounts,(unsigned int)0, size);
	initMap(Qij,0.0, size);

	alphaSum = sum(alpha);
	cSum = 0;

	//------------------------------------------------------------------------
	// experimental!!!

	pseudo1 = new Run("pseudoCountRun1", p->numOfBlocks, 0,0);
		initMap(pseudo1->observations,(U32) 1,size);
		pseudo1->observationSum = size;

	pseudo2 = new Run("pseudoCountRun2", p->numOfBlocks, 0,0);
		initMap(pseudo2->observations,(U32) 1,size);
		pseudo2->observationSum = size;

	DEBUG(1,"Done Creating DirichletDistribution");
}

DirichletDistribution::DirichletDistribution() {
	weight = 0.0;
	alphaSum = 0.0;
	cSum = 0;
	size = 0;
	param = nullptr;

	// experimental
	pseudo1 = nullptr;
	pseudo2 = nullptr;
}

DirichletDistribution::~DirichletDistribution() {

}

void DirichletDistribution::reset() {
	/*! \brief At each gibbs-sampling-step this function is called.
	 *
	 * The runs are released from this component.
	 * The aggregateCCounts are set to 0.
	 */

	runs.clear();
	set0(aggregateCCounts);
	cj.clear();
	cSum = 0;

	// experimental
	runs.push_back(pseudo1);
	runs.push_back(pseudo2);
}

void DirichletDistribution::resetComplete(unsigned int numOfComps) {
	/*! \brief This function is called at the beginning of a new training.
	 */

	runs.clear();
	set0(aggregateCCounts);

	weight = 1.0/(double)numOfComps;

	initMap(alpha,1.0,size);
	initMap(aggregateCCounts,(unsigned int)0, size);
	initMap(Qij,0.0, size);

	alphaSum = sum(alpha);
	cSum = 0;

	// experimental
	runs.push_back(pseudo1);
	runs.push_back(pseudo2);
}

void DirichletDistribution::calculateWeight(const unsigned int n) {
	/*! \brief Calculates the weight of the component based on how many runs have been sampled into this component
	 * in relation to how many runs exist in total.
	 */

	weight = runs.size() / (double)(n + 2*param->components);
}

void DirichletDistribution::calculateCsum() {
	/*! \brief Calculates the total number of reads of all runs in this component.
	 */

	cSum = 0;
	for(unsigned int k = 0; k < runs.size(); k++) {
		cSum += runs[k]->observationSum;
	}
}

void DirichletDistribution::calculateCCountsAndQij() {
	/*! \brief Calculates the CCounts and Qijs.
	 */

	double pseudoCount = 0.0;
	UUmap::iterator j;
	for(j = aggregateCCounts.begin(); j != aggregateCCounts.end(); j++) {
		for(unsigned int k = 0; k < runs.size(); k++) {
			aggregateCCounts[j->first] += runs[k]->observations[j->first];
		}

//		if(cSum > 0) {
//			Qij[j->first] = (double)( aggregateCCounts[j->first] + 1.0)/(double)( cSum + 1.0*aggregateCCounts.size());	// see also Formula (9) Altschull

//			Qij[j->first] = (double)( aggregateCCounts[j->first])/(double)( cSum);

			Qij[j->first] = (double)( aggregateCCounts[j->first] + pseudoCount)/(double)( cSum + pseudoCount*aggregateCCounts.size());

			DEBUG(4,"Qij[" << j->first << "] = " << aggregateCCounts[j->first] << "/" << cSum);
//		}
//		else {
//			Qij[j->first] = 1.0;
//		}
	}
}

void DirichletDistribution::calculateCj() {
	/*! \brief The j-th element of Cj gives the number of reads for all
	 * runs on that specific transcript-unit.
	 *
	 */


	UDmap::iterator j;
	for(j = alpha.begin(); j != alpha.end(); j++) {
		for(unsigned int k = 0; k < runs.size(); k++) {
			cj[j->first][k] = runs[k]->observations[j->first];
			DEBUG(4,"c" << k << ":" << cj[k]);
		}
	}
}

void DirichletDistribution::calculateMixtureParameters(unsigned int const num){
	/*! \brief Calculates the parameters according to the runs sampled in the bins.
	 *
	 */

	//	i)
	calculateWeight(num);
	//	ii)
	calculateCsum();

	calculateCCountsAndQij();

//	iii)
	calculateCj();

	alphaSum = EalphaSum();
}

double DirichletDistribution::calculateLikelihood(Run *r) {
	/*! \brief Calculates the likelihood for a DirichletDistribution that it created the
	 * data in an Experiment. We use the logarithm of the gamma-function because tgamma is
	 * growing too fast. gamma(100) ~= 10^155.
	 * We later transform back with exp().
	 *
	 * \f$\begin{align*}
		& L = \sum_{k=1}^{n_i} \log[\frac{\Gamma(\alpha^{*}_i)}{\Gamma(\alpha^{*}_i + c^{*k}_i)} \cdot \prod_{j=1}^T \frac{\Gamma(\alpha^{*}_i + {c_j^{k}}_i)}{\Gamma(\alpha^{*}_i)}] \label{eq:logLikelihood} \\
		& = \sum_{k=1}^{n_i} \{\log(\Gamma(\alpha^*_i)) - \log(\Gamma(\alpha^*_i + c^{*k}_i)) + \sum_{j=1}^T [ \log(\Gamma(\alpha^*_i{q_j}_i + {c_j^k}_i)) - \log(\Gamma(\alpha^*_i{q_j}_i))] \}
		\end{align*}\f$
	 */

//	double psC = 1.0;
//	double psCSum = psC*(double)param->numOfBlocks;
	double pIk = 0.0;
	UDmap::iterator j;

	pIk = lgamma(alphaSum) - lgamma(alphaSum + (double)r->observationSum);
	for(j = alpha.begin(); j != alpha.end(); j++){
    		pIk += lgamma(j->second + (double)r->observations[j->first])
    										- (double) lgamma(j->second);
	}

//	pIk = exp(pIk);
    pIk = exp(pIk);

//	pIk = gamma(alphaSum)/gamma(alphaSum + r->observationSum);
//	for(j = alpha.begin(); j != alpha.end(); j++){
//		pIk *= gamma(j->second + r->observations[j->first])/gamma(j->second);
//	}

    // is weight really usefull?
//    return pIk;
	return weight*pIk;
}

double DirichletDistribution::EalphaSum() {
	/*! \brief Calculates the start-value for the newtons-method.
	 *
	 */

	double y = 0;
	double counter;
	double denominator;

	// cStar has dimension of number of Runs
	UUmap cStar;
	UDmap::iterator j;
	for(j =alpha.begin(); j != alpha.end(); j++) {
		cStar += cj[j->first];
	}

	for(j =alpha.begin(); j != alpha.end(); j++) {
		if(Qij[j->first]*(1.0-Qij[j->first]) > 0) {
			double vj = Vj(cj[j->first], cStar, Qij[j->first]);
			DEBUG(3,"vj:" << vj);

			UUmap cStar2;
			cStar2 = cStar;
			dotProduct(cStar2);

			counter = Qij[j->first]*(1.0-Qij[j->first])
							*(double)Emean(cStar2)/(double)pow(Emean(cStar),2)
					- vj;
			denominator = vj -(double)Qij[j->first]*(1.0-Qij[j->first])/(double)Emean(cStar);


			y += (Qij[j->first]*(double)counter/(double)denominator);

		}
	}
	return y;
}

double DirichletDistribution::calculateDistance(Run *r){
	double diff = 0.0;
	UDmap::iterator j;
	for(j = alpha.begin(); j != alpha.end(); j++){
		diff += pow(j->second-r->observations[j->first],2);
	}

	return (1/(diff+1))*weight;
}

void DirichletDistribution::calculateSimpleP(){
	/*! \brief Calculates estimation for a one-component dirichlet-mixture.
	 *
	 * We can use this insted of the more complex calculation because we sampled the runs into bins.
	 * The sampling reduced the problem to this very easy estimation.
	 */


	for(unsigned int k = 0; k < runs.size(); k++) {
		double totalCount =  alphaSum + sum(runs[k]->observations);

		UDmap::iterator j;
		for(j=alpha.begin(); j!=alpha.end(); j++) {
			runs[k]->p[j->first] = (double)((double)runs[k]->observations[j->first] + j->second)/(double)totalCount;
		}

	}
}



