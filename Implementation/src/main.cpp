//============================================================================
// Name        : RNAMv2.cpp
// Author      : Willy Bruhn
// Version     :
// Copyright   : Your copyright notice
// Description : RNAMv2 24.11.16
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include "../headers/Run.h"
#include "../headers/TypeConventions.h"
#include "../headers/ChromosomeInitializer.h"
#include "../headers/Aligner.h"
#include "../headers/Downloader.h"
#include "../headers/ParameterHandler.h"
#include "../headers/debug.h"
#include "../headers/Operators.h"
#include <memory>
#include "../headers/Controller.h"

using namespace std;
ParameterHandler *param;

int main(int argc, char *argv[]) {
    param = new ParameterHandler();
    param->readArguments(argc,argv);
    param->export_parameters();

    Controller *c = new Controller(param);
    c->initialize();

    if (1 == param->createDice){
	c->createDiceFromRuns();
    } else {
	c->algorithm();
    }
    delete c;
    return EXIT_SUCCESS;
}
