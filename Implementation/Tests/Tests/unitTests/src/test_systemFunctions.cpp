#include "../../../../headers/systemFunctions.h"

#include "../../../catch/catch.hpp"
#include <iostream>
//#include "../../../../headers/Operators.h"

using namespace std;

TEST_CASE( "filesWithExtInFolder", "" )
{
    vector<string> files = filesWithExtInFolder("/home/willy/BRAKER/VARUS/Implementation/Tests/Tests/unitTests/systemFunctions/fileExtensions", ".fasta");
    
    for(unsigned int i = 0; i < files.size(); i++){
        //cout << files[i] << endl;
    }

    REQUIRE(files.size() == 3);
    REQUIRE(files[0] == "/home/willy/BRAKER/VARUS/Implementation/Tests/Tests/unitTests/systemFunctions/fileExtensions/bla_1.fasta");

}



