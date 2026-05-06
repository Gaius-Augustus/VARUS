/*
 * systemFunctions.h
 *
 *  Created on: 31.12.2017
 *      Author: willy
 */

#ifndef HEADERS_SYSTEMFUNCTIONS_H_
#define HEADERS_SYSTEMFUNCTIONS_H_

#include <iostream>
#include <algorithm>
#include <string>
#include <dirent.h>


bool has_suffix(const std::string& s, const std::string& suffix);

std::vector<std::string> filesWithExtInFolder(const std::string& path, const std::string &ext);


#endif /* HEADERS_SYSTEMFUNCTIONS_H_ */
