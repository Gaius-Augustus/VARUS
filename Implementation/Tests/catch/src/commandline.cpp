#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string>

#include <sys/wait.h>


#include "../headers/commandline.h"

using namespace std;

std::string exec(const char* cmd) {
    char buffer[128];
    std::string result = "";
    FILE* pipe = popen(cmd, "r");
    if (!pipe) throw std::runtime_error("popen() failed!");
    try {
        while (!feof(pipe)) {
            if (fgets(buffer, 128, pipe) != NULL)
                result += buffer;
        }
    } catch (...) {
        pclose(pipe);
        throw;
    }
    pclose(pipe);
    return result;
}

void commandline::makeCall(const std::string &call, const bool verbose){
 
    this->call = call;
    string s = call + " &> ./Tests/tmp.txt";
    int stat = system(s.c_str());
    
    
    status = WEXITSTATUS(stat);
    
    value = "";
    string line;
    ifstream myfile ("./Tests/tmp.txt");
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
        value += line + '\n';
        }
        myfile.close();
    }

    else cerr << "Unable to open file ./Tests/tmp.txt";
    
    system("rm ./Tests/tmp.txt &> /dev/null");
    
    if(verbose){
        print();
    }
}

void commandline::print(){
    cout << "Call: " << call << " | Status: " << status 
    << "\n---------------------------------------------------------\n"
    << value
    << "\n---------------------------------------------------------\n";
}

int commandline::getStatus(){
    return status;
}

std::string commandline::getValue(){
    return value;
}

std::string commandline::getCall(){
    return call;
}

bool commandline::valContains(const std::string &search){
    return (value.find(search) != std::string::npos);
}
    
