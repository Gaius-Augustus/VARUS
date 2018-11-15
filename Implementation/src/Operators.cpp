/*
 * Operators.cpp
 *
 *  Created on: 28.11.2016
 *      Author: willy
 */

#include "../headers/Operators.h"

void remove_substr(const std::string &sub, std::string &str){
//	std::string t = "Banana Republic";
//	std::string s = "nana";

	std::string::size_type i = str.find(sub);

	if (i != std::string::npos)
	   str.erase(i, sub.length());

}

//unsigned int select_randomly(unsigned int max){
//
//	std::default_random_engine generator((unsigned int)time(0) );
//////	std::uniform_int_distribution<unsigned int> distribution(0,size-1);
////	std::uniform_int_distribution<int> distribution(0,9);
//
////    std::random_device rd;
////    std::mt19937 gen(rd());
//
//	std::mt19937 gen;
//	gen.seed(11212);
////	gen.seed(time(nullptr));
//
//    std::uniform_int_distribution<> dis(0, max);
//
//	int index = dis(gen);
//
//	return index;
//}

void scaleMap(UDmap &m,  const double scaleTo) {
	double scalar = scaleTo/sum(m);
	m *= scalar;
}
