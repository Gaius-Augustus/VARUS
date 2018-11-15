#include "../headers/RNAread.h"

RNAread::RNAread()
{

}

RNAread::~RNAread()
{

}

bool RNAread::UMR()
{
    // in order to be an UMR the read has to mapp to only one chromosom...
    if(1 == transcriptUnits.size())
    {
        // ...and only once!
        if(1 == transcriptUnits.begin()->second)
        {
            return true;
        }
    }
    return false;
}
