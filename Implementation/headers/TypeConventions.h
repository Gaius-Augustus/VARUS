/*
 * TypeConventions.h
 *
 *  Created on: 11.11.2016
 *      Author: willy
 */
#include <unordered_map>
#include <string>

#ifndef TYPECONVENTIONS_H_
#define TYPECONVENTIONS_H_

/* ! \brief: Some typedeffs. We will be using a lot of unsigned ints.
 * Defining these types here once, makes it easier to change the type
 * if it seems necessary.
 *
 * E.g. Consider short as an alternative datatype
 * for int with less memory-usage.
 */

typedef unsigned int U32;	// this is used for the maps
typedef unsigned int uI;	// this is used for non-map related unsigned ints

typedef std::unordered_map<U32,U32> UUmap;

typedef std::unordered_map<std::string,U32> SUmap;

typedef std::unordered_map<U32,std::string> USmap;

typedef std::unordered_map<U32,double> UDmap;

typedef std::unordered_map<std::string,double> SDmap;

#endif /* TYPECONVENTIONS_H_ */
