#include "../headers/RNAread.h"

RNAread::RNAread(){
}

RNAread::~RNAread(){
}

bool RNAread::UMR(){
    // in order to be an UMR the read has to map to only one transcription unit ...
    if (1 == transcriptUnits.size()) {
        // ...and only once!
        if (1 == transcriptUnits.begin()->second){
            return true;
        }
    }
    return false;
}
