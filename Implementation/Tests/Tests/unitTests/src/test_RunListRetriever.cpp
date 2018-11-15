/*
 * test_RunListRetriever.cpp
 *
 *  Created on: 18.12.2017
 *      Author: willy
 */

 #include "../../../catch/catch.hpp"
 #include <stdlib.h>

 #include "../../../catch/src/commandline.cpp"

 using namespace std;

 extern commandline COM;

TEST_CASE( "RunListRetriever - html Format correct", "" )
{
    //string s = exec("./Tests/./test.sh");

    // relative to the commandline.cpp
    // COM.makeCall("pwd");
    // /home/willy/BRAKER/VARUS/Implementation
    
    COM.makeCall("../RunListRetriever/./RunListRetriever.pl --help");

    REQUIRE(COM.getStatus() == 0);
    REQUIRE(COM.valContains("Usage:"));

    //COM.print();

}
