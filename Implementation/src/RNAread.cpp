#include "../headers/RNAread.h"

RNAread::RNAread(){

}

RNAread::~RNAread(){

}

bool RNAread::UMR(){
    // in order to be an uniquely mapping read (UMR) the read has to map to only one transcript unit
    // Paired reads are given the same name and currently treated as identical.
    // Therefore, we allow 2 matches to the same transcript unit.
    // TODO: The treatment of read pairs could be more precise.
    if (1 == transcriptUnits.size()                 // at most one transcript unit
	&& transcriptUnits.begin()->second <= 2)   // ... at most twice
	return true;
    return false;
}
