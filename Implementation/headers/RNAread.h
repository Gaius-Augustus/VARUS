#ifndef RNAread_H
#define RNAread_H
#include <unordered_map>
#include <string>

class RNAread {
public:
    std::unordered_map<std::string,unsigned int> transcriptUnits;
    // key = referenceGen, value = number of times matched to the referenceGen

    bool UMR();

    RNAread();
    ~RNAread();
};

#endif // RNAread_H
