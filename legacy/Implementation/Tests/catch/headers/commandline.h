#ifndef COMMANDLINE_H_
#define COMMANDLINE_H_


#include <iostream>
using namespace std;

class commandline {
private: 
    int status; 
    std::string value;
    std::string call;
    
    

public:
    
    void makeCall(const std::string &call, const bool verbose = false);    
    void print();
    
    int getStatus();
    std::string getValue();
    std::string getCall();
    
    bool valContains(const std::string &search);
    
};

#endif /* COMMANDLINE_H_ */
