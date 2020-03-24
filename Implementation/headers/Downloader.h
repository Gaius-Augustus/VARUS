/*
 * Downloader.h
 *
 *  Created on: 29.11.2016
 *      Author: willy
 */

#ifndef DOWNLOADER_H_
#define DOWNLOADER_H_
#include "Run.h"
#include <string.h>
#include "ParameterHandler.h"

class Downloader {
public:
    ParameterHandler *param;

    Downloader(ParameterHandler *p);
    virtual ~Downloader();

    std::string shellCommand(Run *r);
    char const ** shellCommand2(Run *r);
    void nextBatchIndices(Run *r);
    bool getBatch(Run *r, bool all = false);
};

#endif /* DOWNLOADER_H_ */
