/*
 * main_Test.cpp
 *
 *  Created on: 25.11.2016
 *      Author: willy
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <fstream>

#include "../headers/ChromosomeInitializer.h"
#include "../headers/TypeConventions.h"
#include "../headers/debug.h"
#include "../headers/Operators.h"
#include "../headers/ParameterHandler.h"
#include "../headers/Warning.h"

#include "../headers/Run.h"
#include "../headers/Downloader.h"
#include "../headers/Aligner.h"
#include "../headers/RNAread.h"

#include "../headers/Controller.h"

#include "../headers/myRandomEngine.h"

#include <algorithm>

#include "../headers/Estimator.h"

#include "../headers/Simulator.h"
#include "../headers/Functions.h"
#include "../headers/DirichletDistribution.h"
#include "../headers/DirichletMixture.h"

#include "../headers/Functions.h"

#include "../headers/ClusterComponent.h"
#include "../headers/ClusterEstimator.h"

#include "sys/types.h"
#include "sys/sysinfo.h"

ParameterHandler *param = new ParameterHandler();

using namespace std;

unsigned int error_count = 0;
#define ASSRT(val) if (val == false) {					\
	cerr << "->\terror in " << __FUNCTION__ << "()"<< endl;	\
	error_count++;										\
}														\

void printToFile(string path, string content){
	ofstream f;
	f.open(path, ios::out);
//	if(!f) { WARN("Can't open" << path);}
	f << content;
	f.close();
}


string read_file(string path) {

	string file;

	ifstream filestr;
	filestr.open (path);
	if (!filestr.is_open())
	{
//		WARN("Can't open " << path);
	}


	string line;
	while (getline (filestr,line))
	{
		if(!line.empty())
		{
			file += line + "\n";
		}
	}
	filestr.close();

	return file;
}

void test_map_plus() {
	UUmap m;
	m[0] = 1;
	m[1] = 12;
	m[14] = 7;

	UUmap m2;
	m2[0] = 100000;
	m2[2] = 10;
	m2[14] = 3;

	UUmap m3 = m + m2;

	ASSRT((m3[0] == 100001 && m3[2] == 10 && m3[1] == 12 && m3[14] == 10));
}

void test_map_plus_equal() {
	UUmap m;
	m[0] = 1;
	m[1] = 12;
	m[14] = 7;

	UUmap m2;
	m2[0] = 100000;
	m2[2] = 10;
	m2[14] = 3;

	m += m2;

	ASSRT(m[0] == 100001 && m[1] == 12 && m[14] == 10);
}

void test_initMap() {
	unordered_map<unsigned int, unsigned int> a;
	unsigned int size = 10;
	unsigned int num = 1;
	initMap(a,num,size);

	ASSRT(a.size() == 10);
	ASSRT(a[0] == 1);

	{
	unordered_map<unsigned int, double> b;
	unsigned int size = 10;
	double num = 1;
	initMap(b,num,size);

	ASSRT(b.size() == 10);
	ASSRT(b[0] == 1);
	}

	{
	unordered_map<unsigned int, double> c;
	unsigned int size = 10;
	double num = 4.3;
	initMap(c,num,size);

	ASSRT(c.size() == 10);
	ASSRT(c[0] == 4.3);
	}
}

//void test_select_randomly(){
////	for(unsigned int i = 0; i < 100; i++){
////		cerr << select_randomly(10) << endl;
////	}
//
//	std::mt19937 gen;
//	gen.seed(11212);
////	gen.seed(time(nullptr));
//
//	int max = 1;
//    std::uniform_int_distribution<> dis(0, max);
//
//    stringstream comp;
//    for(int i = 0; i < 100; i++){
//    	int index = dis(gen);
//    	cerr << index;
//    	comp << index;
//    }
//    cerr << endl;
//    cerr << comp.str() << endl;
//    string shouldBe = "0100101011110001010110001010101111001111011110101011010110001000111101000100000010111111000101000011";
//    ASSRT(shouldBe==comp.str());
//}

//void test_Vec_select_randomly(){
////	vector<int> vec;
////	for(int i = 0; i < 10; i++) vec.push_back(i);
////
////	cerr << select_randomly(vec) << endl;
//}

void test_DotProduct(){
	UUmap m;
	m[0] = 10;
	m[1] = 1;
	m[2] = 4;
	dotProduct(m);

	stringstream s;
	s << m;

	ASSRT(s.str() == "2:16;1:1;0:100;");

}

void test_ScaleMap(){
	{
		UDmap m;
		m[0] = 10;
		m[1] = 2;

		scaleMap(m,6);

		stringstream s;
		s << m;

		ASSRT("1:1;0:5;" == s.str());
	}
	{
		UDmap m;
		m[0] = 0.9;
		m[1] = 0.1;

		scaleMap(m,85.4);

		stringstream s;
		s << m;

		ASSRT("1:8.54;0:76.86;" == s.str());
	}

}

void test_operators(){
	cerr << "testing operators" << endl;

	test_map_plus();
	test_map_plus_equal();
	test_initMap();
	test_DotProduct();
	test_ScaleMap();
//	test_select_randomly();
//	test_Vec_select_randomly();
}



// myRandomEngine

void test_myRandomEnginge_vectorSelect(){
	vector<int> vec;

	for(int i = 0; i < 10; i++) vec.push_back(i);

	myRandomEngine *eng = new myRandomEngine();
	eng->seed(1241);

    stringstream comp;

    for(int i = 0; i < 100; i++){
    	comp << eng->select_randomly(vec);
    }
//	cerr << comp.str() << endl;
    string shouldBe = "91131326894743745094822007601508195"
    		"824280484309840864733951574673251898412867940"
    		"70937267208348066578";
    ASSRT(shouldBe==comp.str());
}

void test_myRandomEngine_shuffle(){
	myRandomEngine *eng = new myRandomEngine();

	vector<int> vec;
	for(int i = 0; i < 30; i++) vec.push_back(i);

	eng->seed(121313);
	eng->shuffle(vec);

	stringstream comp;
	comp << vec;
//	cerr << vec << endl;

    string shouldBe = "0:22;1:28;2:23;3:24;4:4;5:2;"
    		"6:10;7:27;8:8;9:3;10:6;11:5;12:17;13:29"
    		";14:14;15:1;16:18;17:21;18:25;19:12;20:0;"
    		"21:16;22:26;23:20;24:11;25:19;26:15;27:13;28:9;29:7;";
    ASSRT(shouldBe==comp.str());

}

void test_myRandomEngine_selectRandomly_p(){
	{
		vector<double> probabilities;
		probabilities.resize(2);

		probabilities[0] = 0.1;
		probabilities[1] = 0.9;

		myRandomEngine *r = new myRandomEngine();
		r->seed(19);


		vector<unsigned int> outcomes;
		outcomes.resize(2);
		outcomes[0] = 0;
		outcomes[1] = 0;

		for(unsigned int i = 0; i < 1000; i++){
			outcomes[r->select_randomly_p(probabilities)]++;
		}

		stringstream s;

		s << outcomes;

		ASSRT("0:89;1:911;" == s.str());
	}
	//------------------------------------------------------------
	{
		vector<double> probabilities;
		probabilities.resize(2);

		probabilities[0] = 0.05;
		probabilities[1] = 0.45;

		myRandomEngine *r = new myRandomEngine();
		r->seed(19);


		vector<unsigned int> outcomes;
		outcomes.resize(2);
		outcomes[0] = 0;
		outcomes[1] = 0;

		for(unsigned int i = 0; i < 1000; i++){
			outcomes[r->select_randomly_p(probabilities)]++;
		}

		stringstream s;

		s << outcomes;
		ASSRT("0:89;1:911;" == s.str());
	}
	//------------------------------------------------------------
	{
		vector<double> probabilities;
		probabilities.resize(2);

		probabilities[0] = 0.000000001;
		probabilities[1] = 0.000000009;

		myRandomEngine *r = new myRandomEngine();
		r->seed(19);


		vector<unsigned int> outcomes;
		outcomes.resize(2);
		outcomes[0] = 0;
		outcomes[1] = 0;

		for(unsigned int i = 0; i < 1000; i++){
			outcomes[r->select_randomly_p(probabilities)]++;
		}

		stringstream s;

		s << outcomes;
		ASSRT("0:89;1:911;" == s.str());
	}

}

void test_myRandomEngine(){
	cerr << "testing myRandomEngine" << endl;
//
//	myRandomEngine *eng = new myRandomEngine();
//
//	eng->seed(0);
//
//    stringstream comp;
//
//    for(int i = 0; i < 100; i++){
//    	int index = eng->select_randomly(100);
//    	comp << index;
//    }
////    cerr << comp.str() << endl;
//    string shouldBe = "5559728560865585426265384430905972"
//    		"738487982534857399384734865237849678148787984"
//    		"780804652786811726458145495765210414726187874"
//    		"462157131326215612262399591684536614491701069"
//    		"7676567172136137531613632";
//    ASSRT(shouldBe==comp.str());
//
//
//    eng->seed(12113);
//
//    stringstream comp2;
//
//    for(int i = 0; i < 100; i++){
//    	int index = eng->select_randomly(100);
//    	comp2 << index;
//    }
////    cerr << comp2.str() << endl;
//    string shouldBe2 = "34802029731841035158281427474"
//    		"8434462918469534100449371419148965529169"
//    		"7844799636449621418131140121209633124763"
//    		"8036241888536147767661217029821328893355"
//    		"837982826655257445142299901765500617789";
//    ASSRT(shouldBe2==comp2.str());

//    test_myRandomEnginge_vectorSelect();
//    test_myRandomEngine_shuffle();

	test_myRandomEngine_selectRandomly_p();
}

// --------------------------------------------------------------------------------
// ParameterHandler

void test_ParameterHandler_read_parameters() {

	// somehow the other instances of ParameterHandler crash later
	// if i use a unique_ptr here once.

////	std::shared_ptr<Parameters> parameters = std::make_shared<Parameters>();
//	std::unique_ptr<ParameterHandler> input(new ParameterHandler());
////	ParameterHandler *input = new ParameterHandler();
////	input->read_parameters_from_file("/home/willy/workspace/UnitTests/ParameterHandler");
//	// NEEDS IMPROVEMENT --> if you only add a foldername this thing doesnt work!!!
//
//	input->read_parameters_from_file("/home/willy/workspace/RNAMv2/UnitTests/ParameterHandler/parameters");
//
//	ostringstream ss;
//	input->printParameters(ss);
//
//	string str = ss.str();
//	string t = read_file("/home/willy/workspace/RNAMv2/UnitTests/ParameterHandler/parameters");
//	input->export_parameters();
//
////	cout << str << endl;
////	cout << "second" << endl;
////	cout << t << endl;
//
////	assert(files_equal("/home/willy/workspace/UnitTests/ParameterHandler/parameters1.txt",
////					"/home/willy/workspace/UnitTests/ParameterHandler/parameters"));
//	ASSRT(t == str);

	//-----------------------------------------------------------------------------------------
	ParameterHandler *inp = new ParameterHandler();



	inp->read_parameters_from_file("/home/willy/workspace/RNAMv2/UnitTests/ParameterHandler/parameters");

	delete inp;

	//
//	ostringstream ss;
//	inp->printParameters(ss);
//
//	string str = ss.str();
//	string t = read_file("/home/willy/workspace/RNAMv2/UnitTests/ParameterHandler/parameters");
//	inp->export_parameters();
//
////	cout << str << endl;
////	cout << "second" << endl;
////	cout << t << endl;
//
////	assert(files_equal("/home/willy/workspace/UnitTests/ParameterHandler/parameters1.txt",
////					"/home/willy/workspace/UnitTests/ParameterHandler/parameters"));
//	ASSRT(t == str);



//	ParameterHandler *par1 = new ParameterHandler();

	ParameterHandler par1;
//	par1.read_parameters_from_file("/home/willy/workspace/RNAMv2/UnitTests/ParameterHandler/parameters");


//	ParameterHandler par2;
//	par2.read_parameters_from_file("/home/willy/workspace/RNAMv2/UnitTests/ParameterHandler/parameters");


//
//	cerr << "par" << endl;
//
//	// found this will crash later


//	ParameterHandler *par = new ParameterHandler();
//	par->read_parameters_from_file(
//			"/home/willy/workspace/RNAMv2/UnitTests/Controller/allgorithm/parameters");
}

void test_ParameterHandler_commandPrompt(string commandPrompt) {

	string tmp;
	stringstream ssin(commandPrompt);

	vector<string> wrds;
	unsigned int num = 0;
	while (ssin >> tmp)
	{
		wrds.push_back(tmp);
	}

	int argc_test = wrds.size();

	char *argv_test[argc_test];

	for(unsigned int i = 0; i < wrds.size(); i++) {
		char *p = new char[wrds[i].size() + 1];
		strncpy(p, wrds[i].c_str(), wrds[i].size());
		p[wrds[i].size()] = '\0';
		argv_test[i] = p;
	}

	for(unsigned int i = 0; i < argc_test; i++) {
		cout << "argvInMain :" << argv_test[i] << endl;
	}

//	ParameterHandler *input = new ParameterHandler();
//	input->readArguments(argc_test,argv_test);
//	input->printParameters();

}

void test_ParameterHandler_start_script(std::string path)
{
//	DEBUG(1,"calling script " << path);
    const char * c = path.c_str();
    int status = system(c);
    ASSRT(status==0);
}

void test_ParameterHandler() {
	cerr << "testing ParameterHandler" << endl;

	test_ParameterHandler_read_parameters();

	cerr << "finnished tests on ParameterHandler" << endl;

//	// tests normal execution
//	test_ParameterHandler_start_script("/home/willy/workspace/RNAMv2/UnitTests/ParameterHandler/script/./prog "
//									"--readParametersFromFile 1 "
//						" --pathToParameters /home/willy/workspace/RNAMv2/UnitTests/ParameterHandler/script/"
//						" --verbosityDebug 2");

//	// tests what happens if a parameter is missing in the parametersfile
//	test_ParameterHandler_start_script("/home/willy/workspace/RNAMv2/UnitTests/ParameterHandler/script2/./prog "
//									"--readParametersFromFile 1 "
//									"--exportParametersToFile 0 "
//						" --pathToParameters /home/willy/workspace/RNAMv2/UnitTests/ParameterHandler/script2/"
//						" --verbosityDebug 2");
}

// Run

void test_Run() {
	cerr << "testing Run" << endl;
	unique_ptr<Run> r(new Run());
}

// ChromosomInitializer

void test_ChromosomeInitializer_readInputRuns(){
	cerr << "testing ChromosomeInitializer" << endl;

	ParameterHandler *p = new ParameterHandler();
	p->pathToRuns = "/home/willy/workspace/RNAMv2/UnitTests/ChromosomeInitializer/dielist/";
	//p->dieList = 1;
	unique_ptr<ChromosomeInitializer> chr(new ChromosomeInitializer(p));
	chr->readInputRuns();

	ostringstream ss;
	ss << chr->inputRuns;
	string str = ss.str();

	ASSRT("0:die1;1:die2;2:die36363;" == str);

	p->pathToRuns = "/home/willy/workspace/RNAMv2/UnitTests/ChromosomeInitializer/runlist/";
	//p->dieList = 0;
	unique_ptr<ChromosomeInitializer> chr2(new ChromosomeInitializer(p));
	chr2->readInputRuns();

	ostringstream ss2;
	ss2 << chr2->inputRuns;
	string str2 = ss2.str();


	ASSRT("0:ERR198932;1:ERR198987;2:ERR198889;3:ERR198955;4:ERR199005;5:SRR3720840;6:ERR198948;"
			"7:SRR3720841;8:ERR198980;9:ERR198936;10:SRR3720856;11:SRR3720842;12:SRR3720857;13:SRR3729591;"
			"14:SRR3720851;15:ERR198861;16:SRR3720846;17:SRR3720847;18:SRR3720845;19:SRR3729590;20:SRR3720844;"
			"21:SRR3720843;22:SRR2504046;23:SRR3622109;24:SRR3622110;25:SRR3622111;26:SRR3622112;27:SRR3622114;"
			"28:SRR3622115;29:SRR3622117;30:SRR3622118;31:SRR3622122;32:SRR3622123;33:SRR4044395;34:SRR4044396;"
			"35:SRR4044397;36:SRR001813;37:SRR005108;38:SRR005109;39:SRR006164;40:SRR001337;41:SRR001346;"
			"42:SRR001344;43:SRR2504046;44:SRR3622109;45:SRR3622110;46:SRR3622111;47:SRR3622112;48:SRR3622114;"
			"49:SRR3622115;50:SRR3622117;51:SRR3622118;52:SRR3622122;53:SRR3622123;54:SRR4044395;55:SRR4044396;"
			"56:SRR4044397;57:SRR001813;58:SRR005108;59:SRR005109;60:SRR006164;61:SRR001337;62:SRR001346;"
			"63:SRR001344;64:SRR001348;65:SRR001341;66:SRR001343;67:SRR001664;68:SRR001339;69:SRR014358;"
			"70:SRR014395;71:SRR014396;72:SRR015372;73:SRR016821;74:SRR016822;75:SRR016825;76:SRR016828;"
			"77:SRR016832;78:SRR016834;79:SRR016837;80:SRR019720;81:SRR015074;82:SRR015104;83:SRR037091;"
			"84:SRR037092;85:SRR039383;86:SRR061685;87:SRR064507;88:SRR073282;89:SRR073284;90:SRR073285;"
			"91:SRR073286;92:SRR073289;93:SRR073290;94:SRR073291;95:SRR073292;96:SRR073279;97:SRR077420;"
			"98:SRR023540;99:SRR023549;100:SRR023595;101:SRR023505;102:SRR023609;103:SRR023602;" == str2);
}

void test_ChromosomeInitializer_getChromosomeLengths(){
	ParameterHandler *p = new ParameterHandler();
	p->genomeDir = "/home/willy/workspace/RNAMv2/UnitTests/ChromosomeInitializer/chromosomeLengths/";
	unique_ptr<ChromosomeInitializer> chr(new ChromosomeInitializer(p));
	chr->getChromosomeLengths();

	ostringstream ss;
	ss << chr->chromosomLengths;
	string str = ss.str();

	ASSRT("211000022279446:13234;211000022278908:2389;211000022279696:1464;211000022278306:1710;"
			"211000022278608:2201;211000022278356:1272;211000022278406:1926;211000022280269:1962;"
			"211000022280763:12001;211000022278049:27456;211000022278396:3553;211000022278446:1945;"
			"211000022279406:1186;211000022278416:1157;211000022278526:1067;211000022278256:1380;"
			"211000022280476:6900;211000022278296:1163;211000022280596:13501;211000022278578:2412;"
			"211000022278786:1504;211000022278032:6936;211000022278166:1178;"
			"3Cen_mapped_Scaffold_31_D1643_D1653_D1791:87365;211000022279676:12187;"
			"211000022280446:3970;211000022278888:2200;211000022279114:14983;211000022278276:1063;"
			"211000022279016:1001;211000022279336:1094;211000022280437:5968;211000022280187:13079;"
			"3Cen_mapped_Scaffold_36_D1605:36913;211000022278493:2442;211000022278196:1214;"
			"211000022280046:6062;211000022278343:2113;211000022278136:1232;211000022278286:1446;"
			"211000022278096:1148;211000022278426:2636;211000022280742:12513;211000022278216:2076;"
			"211000022279735:7722;211000022278415:25840;211000022279516:1323;211000022279156:1227"
			";211000022280126:1225;211000022278676:1378;211000022280176:2129;211000022278206:2090"
			";3Cen_mapped_Scaffold_50_D1686:23238;Unmapped_Scaffold_54_D1776:12632;211000022278926"
			":4465;211000022280270:13108;211000022278498:13940;211000022279576:1117;2110000222797"
			"56:2117;211000022278366:1416;211000022278376:1594;211000022278456:1475;2110000222792"
			"46:1538;211000022278536:2304;211000022278546:2668;211000022278386:1029;2110000222798"
			"98:2516;211000022278996:1043;211000022280008:3209;211000022280655:5952;2110000222798"
			"56:1261;211000022278946:1158;211000022278106:1058;211000022278576:1001;2110000222781"
			"85:12681;211000022278686:11126;211000022278716:2383;211000022278776:1076;21100002227"
			"8853:2166;211000022278966:4225;211000022279056:2547;211000022279126:2487;21100002227"
			"9726:1154;211000022279356:1031;211000022279634:5281;211000022278176:1571;21100002227"
			"9396:4370;211000022279416:2795;211000022279436:2421;211000022279456:12095;2110000222"
			"79656:1274;3Cen_mapped_Scaffold_41_D1641:22604;Unmapped_Scaffold_4_D1555_D1692:86267"
			";211000022278317:1970;211000022279476:2187;211000022279506:1053;211000022279526:1023;"
			"211000022279606:3284;211000022279666:1115;211000022279626:1016;" == str);
}

void test_ChromosomeInitializer_create_transcript_units(){
	ParameterHandler *p = new ParameterHandler();
	p->blockSize = 4;
	unique_ptr<ChromosomeInitializer> chr(new ChromosomeInitializer(p));

	chr->chromosomLengths["chr1"] 	= 12;	// 3
	chr->chromosomLengths["X"] 		= 17;	// 5
	chr->chromosomLengths["Y"] 		= 3;	// 1
	chr->chromosomLengths["R"] 		= 14;	// 4


	UUmap totalObservations;
	chr->create_transcript_units(totalObservations);

//	cerr << totalObservations;
//	cerr << endl;
//
//	cerr << chr->transcriptUnitsByNames;
//	cerr << endl;

	ostringstream ss;
	ss << chr->transcriptUnitsByNames;
	string str = ss.str();

	ASSRT(totalObservations.size() == 13
	      && "chr1:0:4;R:3:2;R:1:4;chr1:1:4;R:2:4;Y:0:3;chr1:2:4;X:0:4;"
	      "R:0:4;X:2:4;X:1:4;X:3:4;X:4:1;" == str);
}

void test_ChromosomeInitializer_initializeSigma() {
	ParameterHandler *p = new ParameterHandler();
	unique_ptr<ChromosomeInitializer> chr(new ChromosomeInitializer(p));

	Run* run = new Run("testRun", 30000, atoi("1600000"), 100000);

//	chr->initializeSigma(run);

//	cerr << run->sigma << endl;

//	ASSRT(run->sigma.size() == 16);

	Run* run2 = new Run("testRun2", 30000, atoi("1300100"), 100000);

	vector<Run*> runs;
	runs.push_back(run);
	runs.push_back(run2);

	chr->ran->seed(125421);
	for(unsigned int i = 0; i < runs.size(); i++){
		chr->initializeSigma(runs[i]);
	}

	ASSRT(runs[0]->sigma.size() == 16);
	ASSRT(runs[1]->sigma.size() == 14);

//	cerr << runs[0]->sigma << endl;

	stringstream comp1;
	comp1 << runs[0]->sigma;

    string shouldBe1 = "0:1;1:4;2:3;3:15;4:13;5:11;6:8;"
    		"7:9;8:5;9:10;10:2;11:12;12:6;13:0;14:14;15:7;";
    ASSRT(shouldBe1==comp1.str());


	stringstream comp2;
	comp2 << runs[1]->sigma;
//	cerr << runs[1]->sigma << endl;

    string shouldBe2 = "0:4;1:8;2:10;3:6;4:12;5:0;6:1;"
    		"7:7;8:5;9:9;10:13;11:3;12:11;13:2;";
    ASSRT(shouldBe2==comp2.str());


}

void test_ChromosomeInitializer_initializeRuns(){
	ParameterHandler *p = new ParameterHandler();
	p->pathToRuns = "/home/willy/workspace/RNAMv2/UnitTests/ChromosomeInitializer/initializeRuns";
	unique_ptr<ChromosomeInitializer> chr(new ChromosomeInitializer(p));

	vector<Run*> runs;

	{
		chr->initializeRuns(runs);

		ostringstream ss;
		for(unsigned int i = 0; i < runs.size(); i++){
			ss << i << ":" << runs[i]->accesionId << ";";
		}
		string str = ss.str();

		ASSRT(runs.size() == 10);
		ASSRT("0:DRR016435;1:DRR016443;2:DRR016442;3:DRR016441;4:DRR016440;"
				"5:DRR016439;6:DRR016438;7:DRR016437;8:DRR016436;9:DRR016444;" == str);
	}

}

void test_ChromosomeInitializer() {
	test_ChromosomeInitializer_readInputRuns();
	test_ChromosomeInitializer_getChromosomeLengths();
	test_ChromosomeInitializer_create_transcript_units();

	test_ChromosomeInitializer_initializeSigma();
	test_ChromosomeInitializer_initializeRuns();
}

// Downloader

void test_Downloader_nextBatchIndices(){
	Run* run = new Run("testRun", 30000, 1600000, 100000);
	ParameterHandler *p = new ParameterHandler();
	unique_ptr<ChromosomeInitializer> chr(new ChromosomeInitializer(p));
	chr->initializeSigma(run);


	Downloader *d = new Downloader(p);
	d->nextBatchIndices(run);

//	cerr << run->N << ":" << run->X << endl;
}

void test_Downloader_shellCommand(){
	ParameterHandler *p = new ParameterHandler();
	Downloader *d = new Downloader(p);

	Run *r = new Run();
	string s = d->shellCommand(r);


	ASSRT(s == "fastq-dump -N 0 -X 0 -O defaultnoId/N0X0/ --fasta noId");
}

void test_Downloader_getBatch(){
	Run* run = new Run("DRR016435", 30000, 1600000, 100);
	ParameterHandler *p = new ParameterHandler();
	unique_ptr<ChromosomeInitializer> chr(new ChromosomeInitializer(p));
	chr->initializeSigma(run);

	p->outFileNamePrefix = "/home/willy/workspace/RNAMv2/UnitTests/Downloader/getBatch/";


	Downloader *d = new Downloader(p);
	d->nextBatchIndices(run);

	d->getBatch(run);

//	d->getBatch(run,true);
}

void test_Downloader(){
	cerr << "testing Downloader" << endl;

	test_Downloader_nextBatchIndices();
	test_Downloader_shellCommand();
	test_Downloader_getBatch();
}

// RNAread

void test_RNAread_UMR(){
	RNAread r1;

	r1.transcriptUnits["X1"] = 1;
	r1.transcriptUnits["X2"] = 2;

	ASSRT(r1.UMR() == false);

	RNAread r2;

	r2.transcriptUnits["X1"] = 3;

	ASSRT(r2.UMR() == false);

	RNAread r3;

	r3.transcriptUnits["XY"] = 1;

	ASSRT(r3.UMR() == true);

}

void test_RNAread(){
	cerr << "testing RNAread" << endl;

	test_RNAread_UMR();
}

// Aligner

void test_Aligner_shellCommand(){
	ParameterHandler *p = new ParameterHandler();
	Aligner *a = new Aligner(p);

	Run *r = new Run();
	string s = a->shellCommand(r);



	ASSRT(s == "defaultSTAR --runThreadN 4 --genomeDir default "
			"--readFilesIn defaultnoId/N0X0/noId.fasta --outFileNamePrefix defaultnoId/N0X0/");
}

void test_Aligner_mapReads(){
	ParameterHandler *p = new ParameterHandler();
	Aligner *a = new Aligner(p);

//	N315100X315199
	Run *r = new Run();
	r->N = 315100;
	r->X = 315199;
	r->accesionId = "DRR016435";

	p->outFileNamePrefix = "/home/willy/workspace/RNAMv2/UnitTests/Aligner/getAlignedReads/";
	p->genomeDir = "/home/willy/Bachelorarbeit/Manager/genome/";
	p->pathToSTAR = "/home/willy/Bachelorarbeit/STAR-2.5.1b/source/";
	p->deleteLater = 0;

	a->mapReads(r);
}

void test_Aligner_getAlignedReads(){
	ParameterHandler *p = new ParameterHandler();
	Aligner *a = new Aligner(p);

	Run *r = new Run();

	r->N = 315100;
	r->X = 315199;
	r->accesionId = "DRR016435";

	p->outFileNamePrefix = "/home/willy/workspace/RNAMv2/UnitTests/Aligner/getAlignedReads/";
	p->genomeDir = "/home/willy/Bachelorarbeit/Manager/genome/";
	p->pathToSTAR = "/home/willy/Bachelorarbeit/STAR-2.5.1b/source/";
	p->deleteLater = 0;

	unordered_map<string,RNAread> reads;
	a->getAlignedReads(reads,r);

	stringstream s;
	unordered_map<string,RNAread>::iterator i;
	for(i = reads.begin(); i != reads.end(); i++)
	{
		s << i->second.UMR() << ":";
	}

//	cerr << s.str() << endl;

	ASSRT(s.str() == "1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1"
	":1:1:1:1:1:1:1:1:1:1:1:0:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1"
	":1:1:1:1:1:1:1:1:1:1:1:1:0:1:1:1:1:1:1:1:1:1:1:1:1:1:");

}

void test_Aligner_updateObservations(){
	ParameterHandler *p = new ParameterHandler();
	ChromosomeInitializer *c = new ChromosomeInitializer(p);
	Aligner *a = new Aligner(p);

	unique_ptr<Run> r(new Run());

	unordered_map<string,RNAread> reads;

	UUmap totalObservations;

//	a->updateObservations(r,totalObservations,reads,c);

}

void test_Aligner_update(){

}

void test_Aligner_checkQuality(){

	ParameterHandler *p = new ParameterHandler();
	ChromosomeInitializer *c = new ChromosomeInitializer(p);
	Aligner *a = new Aligner(p);

	a->checkQuality("/home/willy/BRAKER/VARUS/Implementation/unitTests/Aligner/Log.final.out");
}

void test_Aligner(){
	cerr << "testing Aligner" << endl;


	if(false){
		test_Aligner_shellCommand();
	//	test_Aligner_mapReads();
		test_Aligner_getAlignedReads();
		test_Aligner_updateObservations();	// TODO
		test_Aligner_update();				// TODO
	}

	test_Aligner_checkQuality();
}

// Controller
void test_Controller_allgorithm(){
	ParameterHandler *par = new ParameterHandler();
		par->read_parameters_from_file(
				"/home/willy/workspace/RNAMv2/UnitTests/Controller/allgorithm/parameters");
		par->verbosityDebug = 1;



		par->maxBatches = 1;

//		par->export_parameters();

	Controller *c = new Controller(par);

	// sets the seed for the random-choice of the runs
	c->ran->seed(0);

	// sets the seed for the sigma-order
	c->chrom->ran->seed(0);
	c->initialize();


	c->algorithm();

	stringstream s;
	s << c->totalObservations << endl;

	// used this to create reference
//		printToFile("/home/willy/workspace/RNAMv2/UnitTests/Controller/allgorithm/test1.txt",s.str());



	string shouldBe = read_file("/home/willy/workspace/RNAMv2/UnitTests/Controller/allgorithm/test1.txt");

	ASSRT(shouldBe == s.str());

	// num of uniquely mapped reads
	ASSRT(854 == sum(c->totalObservations));
}

void test_Controller_allgorithm2(){
	ParameterHandler *par = new ParameterHandler();
		par->read_parameters_from_file(
				"/home/willy/workspace/RNAMv2/UnitTests/Controller/allgorithm/parameters");
		par->verbosityDebug = 1;

		par->maxBatches = 1;
		par->estimator = 1;

		Controller *c = new Controller(par);

	// sets the seed for the random-choice of the runs
	c->ran->seed(0);

	// sets the seed for the sigma-order
	c->chrom->ran->seed(0);
	c->initialize();


	c->algorithm();

	stringstream s;
	s << c->totalObservations << endl;

	// used this to create reference
//		printToFile("/home/willy/workspace/RNAMv2/UnitTests/Controller/allgorithm/test1.txt",s.str());



	string shouldBe = read_file("/home/willy/workspace/RNAMv2/UnitTests/Controller/allgorithm/test1.txt");

	ASSRT(shouldBe == s.str());

	// num of uniquely mapped reads
	ASSRT(854 == sum(c->totalObservations));
}

void test_Controller_allgorithmSimulation(){
	ParameterHandler *par = new ParameterHandler();
		par->read_parameters_from_file("/home/willy/workspace/RNAMv2/UnitTests/Controller/simulated/parametersCopy");
		par->verbosityDebug = -1;

		par->maxBatches = 0;
		par->estimator = 1;
		par->simulation = 1;
		par->batchSize = 10;
		par->cost = 0.01;
		par->numOfBlocks = 2;
		par->exportObservationsToFile = 1;
		par->createDice = 0;



	Controller *c = new Controller(par);
	// sets the seed for the random-choice of the runs
		c->ran->seed(1);

//		cerr << "hi" << endl;

		c->initialize();

//		cerr << c->runs[0]->observations << endl;

		c->algorithm();

	stringstream s;
	s << c->totalObservations << endl;

//	cerr << "total: " << s.str() << endl;

	for(unsigned int i=0; i < c->runs.size(); i++){
//		cerr << c->sim->p_hidden[c->runs[i]->accesionId] << endl;
//		cerr << c->runs[i]->p << endl;
//		cerr << c->runs[i]->score << endl;
	}

	// used this to create reference
//		printToFile("/home/willy/workspace/RNAMv2/UnitTests/Controller/allgorithm/test1.txt",s.str());


//	ASSRT("0:0.480198;1:0.519802;" == );


//	string shouldBe = read_file("/home/willy/workspace/RNAMv2/UnitTests/Controller/allgorithm/test1.txt");

//	ASSRT(shouldBe == s.str());

	// num of uniquely mapped reads
//	ASSRT(854 == sum(c->totalObservations));
}

void test_Controller_runsEmpty(){
	ParameterHandler *par = new ParameterHandler();
		par->read_parameters_from_file(
				"/home/willy/workspace/RNAMv2/UnitTests/Controller/allgorithm/parameters");
		par->verbosityDebug = 0;
		par->ignoreReadNum = 0;
		par->simulation = 1;


		Run *r = new Run();
		r->numOfSpots = 10;
		r->batchSize = 10;
		r->accesionId = "run";
		r->observationSum = 0;
		r->observations[0] = 0;
		r->timesDownloaded = 0;

		Controller *c = new Controller(par);

		c->runs.push_back(r);
		c->downloadableRuns.push_back(r);

		c->updateDownloadableRuns();
		cerr << c->downloadableRuns.size() << endl;

		// simulate download
		r->observations[0] = 10;
		r->observationSum = 10;
		r->timesDownloaded++;


		c->updateDownloadableRuns();
		cerr << c->downloadableRuns.size() << endl;

		c->algorithm();


}

void test_Controller_createDice(){

	DEBUG(0,"Testing dieCreator");
	ParameterHandler *par = new ParameterHandler();
	par->read_parameters_from_file(
			"/home/willy/workspace/RNAMv2/UnitTests/Controller/dieCreator/parameters");


	Controller *c = new Controller(par);
		c->initialize();

//		c->createDiceFromRuns();
}

void test_Controller_exitCondition(){
	ParameterHandler *par = new ParameterHandler();
	par->read_parameters_from_file(
			"/home/willy/workspace/RNAMv2/UnitTests/Controller/ExitCondition/parametersCopy");

	par->numOfBlocks = 2;
	par->pseudoCount = 0;
	par->estimator = 1;
	par->simulation = 1;
	par->verbosityDebug = 0;
	par->ignoreReadNum = 1;

	Controller *con = new Controller(par);

	con->initialize();

	con->algorithm();
	cerr << con->maxProfit << endl;

	// -----------------------------------------
	par->loadAllOnce = 1;
	Controller *con2 = new Controller(par);
	con2->initialize();
	con2->algorithm();


}

void test_Controller(){
	cerr << "testing Controller" << endl;

//	test_Controller_allgorithm();

//	test_Controller_allgorithmSimulation();

//	test_Controller_createDice();

	// test_Controller_runsEmpty();

	test_Controller_exitCondition();
}

// DirichletDistribution
void test_DirichletDistribution_const(){
	ParameterHandler *p = new ParameterHandler();
	p->numOfBlocks = 10;
	DirichletDistribution d(p,(double)1/5);

	ASSRT(d.weight == (double)1/5);
	ASSRT(sum(d.alpha) == p->numOfBlocks*1.0);
	ASSRT(sum(d.alpha) == d.alphaSum);

	ASSRT(d.runs.size() == 0);
}

void test_DirichletDistribution_reset(){
	ParameterHandler *p = new ParameterHandler();
	p->numOfBlocks = 10;
	DirichletDistribution d(p,(double)1/5);

	Run *r = new Run();
	d.runs.push_back(r);

	d.reset();

	ASSRT(sum(d.aggregateCCounts) == 0);
	ASSRT(d.runs.size() == 0);

}

void test_DirichletDistribution_calculateMixtureParameters(){
	ParameterHandler *p = new ParameterHandler();
	p->numOfBlocks = 2;
	DirichletDistribution d(p,(double)1);

	vector<Run*> runs;
	runs.resize(2);
	Run *r1 = new Run();
	r1->observations[0] = 10;
	r1->observations[1] = 1;
	r1->observationSum = 11;

	Run *r2 = new Run();
	r2->observations[0] = 10;
	r2->observations[1] = 1;
	r2->observationSum = 11;

	runs[0] = r1;
	runs[1] = r2;

	// the runs have been sampled into component d
	d.runs.push_back(r1);
	d.runs.push_back(r2);
	d.calculateCsum();
	ASSRT(22 == d.cSum);

	d.calculateCCountsAndQij();

	stringstream s;
	s << d.Qij;
	ASSRT(s.str() == "1:0.0909091;0:0.909091;");

	d.calculateCj();

	UDmap::iterator j;
	stringstream s2;
	for(j = d.alpha.begin(); j != d.alpha.end(); j++){
		s2 << d.cj[j->first] << "\n";
	}

	ASSRT(s2.str() == "1:1;0:1;\n1:10;0:10;\n");
}

string test_DirichletDistribution_Mixtures2runs2dim(int x_1, int x_2, int y_1, int y_2, string path, bool print = false){
	Run *r1 = new Run();
	r1->observations[0] = x_1;
	r1->observations[1] = x_2;
	r1->observationSum = sum(r1->observations);

	Run *r2 = new Run();
	r2->observations[0] = y_1;
	r2->observations[1] = y_2;
	r2->observationSum = sum(r2->observations);

	ParameterHandler *p = new ParameterHandler();
	p->numOfBlocks = 2;
	p->outFileNamePrefix = path;

	DirichletDistribution d(p,(double)1);

	d.runs.push_back(r1);
	d.runs.push_back(r2);

	d.calculateMixtureParameters(2);

	DirichletMixture *mix = new DirichletMixture(p);



	d.alphaSum = mix->newtonsMethod(0,d,fabs(d.alphaSum),1000,0.0001);

	d.alpha = d.Qij;
	scaleMap(d.alpha,d.alphaSum);

	stringstream s2;
	s2 << d.alpha;

	if(print){
		cerr << r1->observations << endl;
		cerr << r2->observations << endl;
		cerr << s2.str() << endl;
	}

	return s2.str();
}

void test_DirichletDistribution_MixturesNruns_n_dim(vector<UUmap> &obs,bool print = false){


}

void test_DirichletDistribution_Mixtures(){

	string path = "/home/willy/workspace/RNAMv2/UnitTests/Estimator/DirichletMixture/mixtures";
	//---------------------------------------------------------

//	test_DirichletDistribution_Mixtures2runs2dim(9,1,1,9,path,true);
//
//	//---------------------------------------------------------
//	test_DirichletDistribution_Mixtures2runs2dim(90,10,10,90,true);
//
//
//	//---------------------------------------------------------
//	test_DirichletDistribution_Mixtures2runs2dim(11,10,8,12,true);


	//---------------------------------------------------------
	test_DirichletDistribution_Mixtures2runs2dim(2,1,1,2,path + "/t4",true);


	//---------------------------------------------------------
	test_DirichletDistribution_Mixtures2runs2dim(5,1,1,5,path + "/t5",true);


	//---------------------------------------------------------
	test_DirichletDistribution_Mixtures2runs2dim(4,1,1,4,path + "/t6",true);

	//---------------------------------------------------------
	test_DirichletDistribution_Mixtures2runs2dim(3,1,1,3,path + "/t7",true);

	//---------------------------------------------------------
	test_DirichletDistribution_Mixtures2runs2dim(30,10,10,30,path + "/t7",true);

	//---------------------------------------------------------
	test_DirichletDistribution_Mixtures2runs2dim(90,10,90,10,path + "/t8",true);

	//---------------------------------------------------------
	test_DirichletDistribution_Mixtures2runs2dim(900,100,900,100,path + "/t9",true);

	//---------------------------------------------------------
	test_DirichletDistribution_Mixtures2runs2dim(10,1,1,10,path + "/t9",true);

	//---------------------------------------------------------
	test_DirichletDistribution_Mixtures2runs2dim(8,4,4,8,path + "/t9",true);

	//---------------------------------------------------------
	test_DirichletDistribution_Mixtures2runs2dim(2,0,0,2,path + "/t9",true);

}

void test_DirichletDistribution_EalphaSum(){
	ParameterHandler *p = new ParameterHandler();
	p->numOfBlocks = 2;
	DirichletDistribution d(p,(double)1);

	vector<Run*> runs;
	runs.resize(2);
	Run *r1 = new Run();
	r1->observations[0] = 10;
	r1->observations[1] = 1;
	r1->observationSum = 11;

	Run *r2 = new Run();
	r2->observations[0] = 10;
	r2->observations[1] = 1;
	r2->observationSum = 11;

	runs[0] = r1;
	runs[1] = r2;

	// the runs have been sampled into component d
	d.runs.push_back(r1);
	d.runs.push_back(r2);


	d.calculateMixtureParameters(runs.size());

//	cerr << d.Qij << endl;

	UUmap cStar;
	UDmap::iterator j;
	for(j =d.alpha.begin(); j != d.alpha.end(); j++) {
		cStar += d.cj[j->first];
	}
	stringstream s1;
	s1 << cStar;

	ASSRT("0:11;1:11;" == s1.str());





	j = d.alpha.begin();
	cerr << d.cj[j->first] << endl;
	double vj = (double)Evariance(d.cj[j->first])/(double)pow(Emean(cStar),2) - ((double)Evariance(cStar)/(double)pow(Emean(cStar),2))*(double)pow(d.Qij[j->first],2);

	cerr << vj << endl;

	UUmap cStar2;
	cStar2 = cStar;
	dotProduct(cStar2);

	cerr << cStar2 << endl;

	double counter = d.Qij[j->first]*(1.0-d.Qij[j->first])
					*(double)Emean(cStar2)/(double)pow(Emean(cStar),2)
			- vj;

	cerr << counter << endl;

	double denominator = vj -(double)d.Qij[j->first]*(1.0-d.Qij[j->first])/(double)Emean(cStar);

	cerr << denominator << endl;
	cerr << counter/denominator << endl;

	//-------------------------------------------
	j++;
	cerr << d.cj[j->first] << endl;
	double vj2 = (double)Evariance(d.cj[j->first])/(double)pow(Emean(cStar),2) - ((double)Evariance(cStar)/(double)pow(Emean(cStar),2))*(double)pow(d.Qij[j->first],2);

	cerr << vj2 << endl;

	UUmap cStar22;
	cStar22 = cStar;
	dotProduct(cStar22);

	cerr << cStar22 << endl;

	double counter2 = d.Qij[j->first]*(1.0-d.Qij[j->first])
					*(double)Emean(cStar22)/(double)pow(Emean(cStar),2)
			- vj2;

	cerr << counter2 << endl;

	double denominator2 = vj2 -(double)d.Qij[j->first]*(1.0-d.Qij[j->first])/(double)Emean(cStar);

	cerr << denominator2 << endl;
	cerr << counter2/denominator2 << endl;




	double alp = d.EalphaSum();

	cerr << alp << endl;

}

void test_DirichletDistribution_EalphaSum2(){
	ParameterHandler *p = new ParameterHandler();
	p->numOfBlocks = 2;
	DirichletDistribution d(p,(double)1);

	vector<Run*> runs;
	runs.resize(2);
	Run *r1 = new Run();
	r1->observations[0] = 10;
	r1->observations[1] = 1;
	r1->observationSum = 11;

	Run *r2 = new Run();
	r2->observations[0] = 10;
	r2->observations[1] = 3;
	r2->observationSum = 13;

	runs[0] = r1;
	runs[1] = r2;

	// the runs have been sampled into component d
	d.runs.push_back(r1);
	d.runs.push_back(r2);


	d.calculateMixtureParameters(runs.size());

//	cerr << d.EalphaSum() << endl;
}

void test_DirichletDistribution(){
//	test_DirichletDistribution_reset();
//	test_DirichletDistribution_const();
//	test_DirichletDistribution_calculateMixtureParameters();
//	test_DirichletDistribution_EalphaSum();
	test_DirichletDistribution_Mixtures();
}

// DirichletMixture
void test_DirichletMixture_train(){
	ParameterHandler *p = new ParameterHandler();
		p->numOfBlocks = 2;
		p->trainingsIterations = 13;
		p->newtonIterations = 10;
		p->newtonPrecision = 0.0001;
		p->components = 2;
		p->outFileNamePrefix = "/home/willy/workspace/RNAMv2/UnitTests/Estimator/DirichletMixture/test0/";
//		p->randomSeed = 130;

		vector<Run*> runs;

		for(unsigned int i = 0; i < 10; i++){
			Run *r = new Run();
			r->accesionId = "1";
			r->observations[0] = 100;
			r->observations[1] = 10;
			r->observationSum = 110;
//			r->p[0] = 0.5;
//			r->p[1] = 0.5;
			runs.push_back(r);

			Run *r2 = new Run();
			r2->accesionId = "2";
			r2->observations[0] = 10;
			r2->observations[1] = 100;
			r2->observationSum = 110;
//			r2->p[0] = 0.5;
//			r2->p[1] = 0.5;
			runs.push_back(r2);
		}

		DirichletMixture *m = new DirichletMixture(p);

		ASSRT(2==m->components.size());
//			cerr << m->components[0].alphaSum << endl;
//			cerr << m->components[0].alpha << endl;
//
//			vector<double> likelihoods;
//			m->calculateLikelihoods(runs[0],likelihoods);
//			cerr << likelihoods << endl;
//			cerr << m->sampler->select_randomly_p(likelihoods) << endl;
//
//
//			// samples the k-th dice into the component corresponding to the number generated by the discrete distribution
//			m->components[m->sampler->select_randomly_p(likelihoods)].runs.push_back(runs[0]);


		p->verbosityDebug = 2;

		USmap tr;
		tr[0] = "0";
		tr[1] = "1";
		m->setTranslate2Str(tr);

		m->train(runs);

		int count = 0;
		for(unsigned int i = 0; i < m->components[0].runs.size(); i++){
			if(m->components[0].runs[i]->accesionId == "1") count++;
		}
		cerr << "count:" << count << endl;
}

void test_DirichletMixture_calculateLikelihood(){

	ParameterHandler * p = new ParameterHandler();
		p->numOfBlocks = 2;
		p->components = 2;


	DirichletDistribution d(p,0.5);
		d.alpha[0] = 1;
		d.alpha[1] = 1;
		d.alphaSum = 2;

	DirichletDistribution d2(p,0.5);
		d2.alpha[0] = 1;
		d2.alpha[1] = 1;
		d2.alphaSum = 2;


		vector<Run*> runs;
	for(unsigned int i = 0; i < 5; i++){
	Run *r1 = new Run();
		r1->observations[0] = 10;
		r1->observations[1] = 90;
		r1->observationSum = 100;
		r1->accesionId = "1";
		runs.push_back(r1);

	Run *r2 = new Run();
		r2->observations[0] = 90;
		r2->observations[1] = 10;
		r2->observationSum = 100;
		r2->accesionId = "2";
		runs.push_back(r2);
	}

	cerr << d.calculateLikelihood(runs[0]) << endl;
	cerr << d2.calculateLikelihood(runs[0]) << endl;

	cerr << d.calculateLikelihood(runs[1]) << endl;
	cerr << d2.calculateLikelihood(runs[1]) << endl;

	//-----------------------------------------------

	// first has 3 of "1" and second 3 of "2"
	// there should be an imbalance now
	d.runs.push_back(runs[0]);
	d.runs.push_back(runs[1]);
	d.runs.push_back(runs[2]);
	d.runs.push_back(runs[6]);
	d.runs.push_back(runs[4]);

	d2.runs.push_back(runs[3]);
	d2.runs.push_back(runs[6]);
	d2.runs.push_back(runs[7]);
	d2.runs.push_back(runs[8]);
	d2.runs.push_back(runs[9]);


	d.calculateMixtureParameters(10);

	DirichletMixture *mix = new DirichletMixture(p);



	d.alphaSum = mix->newtonsMethod(0,d,fabs(d.alphaSum),5,0.001);
	d.alpha = d.Qij;
	scaleMap(d.alpha,d.alphaSum);

	cerr << d.alphaSum << endl;
	cerr << d.alpha << endl;
	cerr << d.Qij << endl;

	cerr << d.calculateLikelihood(runs[0]) << endl;
	cerr << d.calculateLikelihood(runs[1]) << endl;

}

void test_DirichletMixture_calculateP(){
	ParameterHandler *p = new ParameterHandler();
		p->numOfBlocks = 2;
		p->trainingsIterations = 12;
		p->newtonIterations = 10;
		p->newtonPrecision = 0.0001;
		p->components = 2;
		p->outFileNamePrefix = "/home/willy/workspace/RNAMv2/UnitTests/Estimator/DirichletMixture/test0/";
		p->randomSeed = 13809;

		vector<Run*> runs;

		for(unsigned int i = 0; i < 4; i++){
			Run *r = new Run();
			r->accesionId = "1";
			r->observations[0] = 100;
			r->observations[1] = 10;
			r->observationSum = 110;
//			r->p[0] = 0.5;
//			r->p[1] = 0.5;
			runs.push_back(r);

			Run *r2 = new Run();
			r2->accesionId = "2";
			r2->observations[0] = 10;
			r2->observations[1] = 100;
			r2->observationSum = 110;
//			r2->p[0] = 0.5;
//			r2->p[1] = 0.5;
			runs.push_back(r2);

			Run *r3 = new Run();
			r3->accesionId = "3";
			r3->observations[0] = 0;
			r3->observations[1] = 0;
			r3->observationSum = 0;
//			r2->p[0] = 0.5;
//			r2->p[1] = 0.5;
			runs.push_back(r3);
		}

		DirichletMixture *m = new DirichletMixture(p);

		ASSRT(2==m->components.size());

		p->verbosityDebug = 2;

		USmap tr;
		tr[0] = "0";
		tr[1] = "1";
		m->setTranslate2Str(tr);

		m->train(runs);

		int count = 0;
		for(unsigned int i = 0; i < m->components[0].runs.size(); i++){
			if(m->components[0].runs[i]->accesionId == "1") count++;
		}



		//-------------------------------------------------------------
		m->calculateP(runs);
		cerr << count << " " << m->components[0].runs.size() << endl;

		Controller *con = new Controller(p);
		con->runs = runs;
		cerr << con->runs.size() << endl;
		con->sim->translate2str = tr;

		con->batchCount = 0;

		p->exportObservationsToFile = 1;
		p->simulation = 1;

		for(unsigned int i = 0; i < con->runs.size(); i++){
			con->totalObservations += runs[i]->observations;
		}

		con->exportTotalObservationCSV("best");



}

void test_DirichletMixture(){
//	test_DirichletMixture_calculateLikelihood();
//	test_DirichletMixture_train();
	test_DirichletMixture_calculateP();
}


// Estimators

void test_simpleEstimator1Run(UUmap obs, double pseudoCount, string solution, bool shouldPrint=false)
{
	ParameterHandler *par = new ParameterHandler();
		par->read_parameters_from_file(
				"/home/willy/workspace/RNAMv2/UnitTests/Estimator/SimpleEstimator/parameters");

		par->numOfBlocks = obs.size();
		par->pseudoCount = pseudoCount;

		vector<Run*> runs;
		Run* r = new Run();
		runs.push_back(r);

	SimpleEstimator *si = new SimpleEstimator(par);

	UUmap total;
	initMap(total,(U32)0,(U32)obs.size());

	si->initializeRuns(runs,total);

	// simulate observations
	r->observations = obs;

	si->estimateP(runs,0);

	stringstream s;
	s << r->observations;

	stringstream s2;
	s2 << obs;

	ASSRT(s.str() == s2.str());

	stringstream s3;
	s3 << r->p;

	ASSRT(s3.str() == solution);

	if(shouldPrint) { cerr << s3.str() << endl;}
}

void test_simpleEstimator(){

	//----------------------------------------------------------------
	{
		UUmap tot1;
		tot1[0] = 1;
		tot1[1] = 0;

		test_simpleEstimator1Run(tot1,1.0,"0:0.666667;1:0.333333;");
	}

	//----------------------------------------------------------------
	{
		UUmap tot;
		tot[0] = 112;
		tot[1] = 16;

		test_simpleEstimator1Run(tot,1.0,"0:0.869231;1:0.130769;");
	}

	//----------------------------------------------------------------
	{
		UUmap tot;
		tot[0] = 123123123;
		tot[1] = 1;

		test_simpleEstimator1Run(tot,0.0,"0:1;1:8.12195e-09;");
	}

	//----------------------------------------------------------------
	{
		UUmap tot;
		tot[0] = 123123123;
		tot[1] = 0;

		test_simpleEstimator1Run(tot,0.0,"0:1;1:0;");
	}

	//----------------------------------------------------------------
	{
		UUmap tot;
		tot[0] = 14;
		tot[1] = 3;
		tot[2] = 4;
		tot[3] = 18;
		tot[4] = 16;

		test_simpleEstimator1Run(tot,1.0,
				"0:0.25;1:0.0666667;2:0.0833333;3:0.316667;4:0.283333;");

	}
}

void test_AdvancedEstimator2Run(UUmap obs, UUmap obs2, double pseudoCount, double lambda, string solution, bool shouldPrint=false)
{
	ParameterHandler *par = new ParameterHandler();
		par->read_parameters_from_file(
				"/home/willy/workspace/RNAMv2/UnitTests/Estimator/SimpleEstimator/parameters");

		par->numOfBlocks = obs.size();
		par->pseudoCount = pseudoCount;
		par->lambda = lambda;
		par->verbosityDebug = 3;

		vector<Run*> runs;
		Run* r = new Run();
		runs.push_back(r);
		Run* r2 = new Run();
		runs.push_back(r2);

	AdvancedEstimator *si = new AdvancedEstimator(par);

	UUmap total;
	initMap(total,(U32)0,(U32)obs.size());

	si->initializeRuns(runs,total);

	// simulate observations
	r->observations = obs;
	r->observationSum = sum(obs);
	r2->observations = obs2;
	r2->observationSum = sum(obs2);

	si->estimateP(runs,0);

	stringstream s;
	s << r->observations;

	stringstream s2;
	s2 << obs;

	ASSRT(s.str() == s2.str());

	stringstream s3;
	s3 << r->p;

	stringstream s4;
	s4 << r2->p;

	ASSRT(s3.str() == solution);

	if(shouldPrint) { cerr << s3.str() << endl;}
	if(shouldPrint) { cerr << s4.str() << endl;}
}

void test_AdvancedEstimator(){

	//----------------------------------------------------------------
		{
			UUmap tot1;
			tot1[0] = 1;
			tot1[1] = 1;

			UUmap tot2;
			tot2[0] = 1;
			tot2[1] = 0;

			test_AdvancedEstimator2Run(tot1,tot2,1.0,1.0,"0:0.533333;1:0.466667;");
		}
//
//		//----------------------------------------------------------------
//		{
//			UUmap tot;
//			tot[0] = 112;
//			tot[1] = 16;
//
//			test_AdvancedEstimator1Run(tot,1.0,0.0,"0:0.869231;1:0.130769;");
//		}
//
//		//----------------------------------------------------------------
//		{
//			UUmap tot;
//			tot[0] = 123123123;
//			tot[1] = 1;
//
//			test_AdvancedEstimator1Run(tot,0.0,0.0,"0:1;1:8.12195e-09;");
//		}
//
//		//----------------------------------------------------------------
//		{
//			UUmap tot;
//			tot[0] = 123123123;
//			tot[1] = 0;
//
//			test_AdvancedEstimator1Run(tot,0.0,0.0,"0:1;1:0;");
//		}
//
//		//----------------------------------------------------------------
//		{
//			UUmap tot;
//			tot[0] = 14;
//			tot[1] = 3;
//			tot[2] = 4;
//			tot[3] = 18;
//			tot[4] = 16;
//
//			test_AdvancedEstimator1Run(tot,1.0,0.0,
//					"0:0.25;1:0.0666667;2:0.0833333;3:0.316667;4:0.283333;");
//
//		}
//
//		//-----------------------------------------------------------------
//		// tests for lambda
//		//-----------------------------------------------------------------
//
//		{
//			UUmap tot1;
//			tot1[0] = 1;
//			tot1[1] = 0;
//
//			test_AdvancedEstimator1Run(tot1,1.0,100.0,"0:0.666667;1:0.333333;");
//		}
//
//		//----------------------------------------------------------------
}

//void test_DirichletMixture(){
//	ParameterHandler *p = new ParameterHandler();
//
//		p->read_parameters_from_file(
//				"/home/willy/workspace/RNAMv2/UnitTests/Estimator/DirichletMixture/test1/parameters");
//
//
//}

void test_Estimators(){
	test_simpleEstimator();
	test_AdvancedEstimator();
//	test_DirichletMixture();
}

// Simulator

void test_Simulator_readInputRuns(){
	ParameterHandler *p = new ParameterHandler();
	p->pathToRuns = "/home/willy/workspace/RNAMv2/UnitTests/Simulator/readInputRuns";
	//p->dieList = 1;
	Simulator *sim = new Simulator(p);
	sim->readInputRuns();

	ostringstream ss;
	ss << sim->inputRuns;
	string str = ss.str();

//	cerr << str << endl;

	ASSRT("0:die1;1:die2;2:die14;" == str);
}

void test_Simulator_initializeRuns(){
	ParameterHandler *p = new ParameterHandler();
	p->pathToRuns = "/home/willy/workspace/RNAMv2/UnitTests/Simulator/initializeRuns";
	p->pathToDice = "/home/willy/workspace/RNAMv2/UnitTests/Simulator/initializeRuns";

	Simulator *sim = new Simulator(p);

	vector<Run*> runs;

	{
		sim->readInputRuns();
		sim->initializeRuns(runs);

		ostringstream ss;
		stringstream ps;

		for(unsigned int i = 0; i < runs.size(); i++){
			ss << i << ":" << runs[i]->accesionId << ";";
			ps << i << ":" << sim->p_hidden[runs[i]->accesionId] << ";";
		}
		string str = ss.str();

//		cerr << str << endl;
//		cerr << ps.str() << endl;

		ASSRT(runs.size() == 2);
		ASSRT("0:die1;1:die2;" == str);

		ASSRT("0:0:0.5;1:0.5;;1:0:0.793548;1:0.206452;;" == ps.str());

	}
}

void test_Simulator_simulateObservations(){

	ParameterHandler *p = new ParameterHandler();
	p->batchSize = 10000;

	Simulator *sim = new Simulator(p);

	UUmap totalObservations;
	totalObservations[0] = 0;
	totalObservations[1] = 0;

	Run *r = new Run();
	r->observations[0] = 0;
	r->observations[1] = 0;
	r->accesionId = "test";

	vector<double> pH;

	pH.resize(2);
	pH[0] = 9;
	pH[1] = 1;

	sim->p_hidden["test"] = pH;

	sim->ran->seed(0);
	sim->simulateObservations(r,totalObservations);

	stringstream s1;
	stringstream s2;

	s1 << totalObservations;
	s2 << r->observations;

//	cerr << totalObservations << endl;
//	cerr << r->observations << endl;

	ASSRT("1:975;0:9025;" == s1.str());
	ASSRT("1:975;0:9025;" == s2.str());

}

void test_Simulator(){
	cerr << "testing Simulator" << endl;

	test_Simulator_readInputRuns();
	test_Simulator_initializeRuns();
	test_Simulator_simulateObservations();
}

// Functions

void test_Digamma() {
	double s = 0.0000001;

	ASSRT(fabs(digamma(0.5) - (-1.96351003)) < s);
	ASSRT(fabs(digamma(1) - (-0.57721566490153286060651209008240243104215933593992359880)) < s);
	ASSRT(fabs(digamma(10) - (2.251752589066721107647456163885851537211808918028330369448)) < s);

}

void test_Trigamma() {
	double s = 0.0001;

	ASSRT(fabs(trigamma(0.5) - (4.9348022005446793094172454999380755676568497036203953)) < s);
	ASSRT(fabs(trigamma(1) - (1.644934066848226436472415166646025189218949901206798437735)) < s);
	ASSRT(fabs(trigamma(10) - (0.105166335681685746122201006908055927440164312897400604528)) < s);

}

void test_Emean() {
	UUmap m;
	m[0] = 13;
	m[1] = 17;
	m[2] = 5;
	m[3] = 12;
	m[4] = 14;

	ASSRT(Emean(m) == 12.2);
}

void test_Evariance() {
	UUmap m;
	m[0] = 13;
	m[1] = 17;
	m[2] = 5;
	m[3] = 12;
	m[4] = 14;

	ASSRT(19.7 == Evariance(m));
}



void test_Functions(){
	test_Digamma();
	test_Trigamma();
	test_Emean();
	test_Evariance();
}

// ClusterEstimator
void test_ClusterEstimator_const(){
	ParameterHandler *p = new ParameterHandler();
	p->numOfBlocks = 10;
	p->components = 2;
	p->pseudoCount = 1;
	p->verbosityDebug = 2;
	p->outFileNamePrefix = "/home/willy/workspace/RNAMv2/UnitTests/Estimator/ClusterEstimator/test0/";

	ClusterEstimator *cl = new ClusterEstimator(p);


	vector<Run*> runs;

			for(unsigned int i = 0; i < 1; i++){
				Run *r = new Run();
				r->accesionId = "1";
				r->observations[0] = 100;
				r->observations[1] = 10;
				r->observationSum = 110;
				runs.push_back(r);

				Run *r2 = new Run();
				r2->accesionId = "2";
				r2->observations[0] = 0;
				r2->observations[1] = 0;
				r2->observationSum = 0;

				runs.push_back(r2);
			}

			p->verbosityDebug = 2;

			cl->estimateP(runs,0);


	cerr << runs[0]->p << endl;
}


static bool pred(Run* s) {
	if(0 == s->observationSum) return true;
	else return false;
}


void test_ClusterRemoveRuns(){
	std::vector<Run*> v;
	for(unsigned int i = 0; i < 10; i++){
		Run * r = new Run();
		r->observationSum = 0;
		v.push_back(r);
	}
	Run * r = new Run();
	r->observationSum = 10;
	v.push_back(r);
	ASSRT(v.size() == 11);

	v.erase( std::remove_if( v.begin(), v.end(), pred ), v.end() );
	ASSRT(v.size() == 1);
}



void test_ClusterEstimator_calculateP(){

	ParameterHandler *p = new ParameterHandler();
	p->numOfBlocks = 3;
	p->lambda = 2;
	p->verbosityDebug = 3;
	ClusterComponent c(p);

	Run *r = new Run();
	r->observations[0] = 10;
	r->observations[1] = 0;
	r->observations[2] = 12;
	r->observationSum = 22;

	r->q[0] = (double)10/22;
	r->q[1] = 0;
	r->q[2] = (double)12/22;


	c.setSeed(r->q);

	c.runs.push_back(r);

	c.calculateAlphaSum();

	c.calculateP();

	cerr << r->p << endl;
}

void test_ClusterEstimator_estimateP(){
	ParameterHandler *p = new ParameterHandler();


}

void test_ClusterEstimator_calculateVariance(){

}

void test_ClusterEstimator(){
	test_ClusterEstimator_const();
	test_ClusterRemoveRuns();
	test_ClusterEstimator_calculateP();
	test_ClusterEstimator_estimateP();
	test_ClusterEstimator_calculateVariance();
}



// Memory tests
void test_MapSize(){
	UUmap m;

	cerr << sizeof(m) << endl;

	for(unsigned int i = 0; i < 24; i++){
		m[i] = i;
	}


	cerr << sizeof(m) << endl;

	cerr << sizeof(m[0]) << endl;
	cerr << sizeof((U32) 0) << endl;


  size_t count = 0;
  for (unsigned i = 0; i < m.bucket_count(); ++i) {
	size_t bucket_size = m.bucket_size(i);
	if (bucket_size == 0) {
	  count++;
	}
	else {
	  count += bucket_size;
	}
  }

  cerr << count << endl;
}


struct sysinfo memInfo;


int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getValue(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}


void RunCrashTest(){

	sysinfo (&memInfo);
	long long totalVirtualMem = memInfo.totalram;
	//Add other values in next statement to avoid int overflow on right hand side...
	totalVirtualMem += memInfo.totalswap;
	totalVirtualMem *= memInfo.mem_unit;

	cerr << totalVirtualMem/1000000000 << endl;


	long long totalPhysMem = memInfo.totalram;
	//Multiply in next statement to avoid int overflow on right hand side...
	totalPhysMem *= memInfo.mem_unit;


	cerr << totalPhysMem/1000000000 << endl;

	long long physMemUsed = memInfo.totalram - memInfo.freeram;
	//Multiply in next statement to avoid int overflow on right hand side...
	physMemUsed *= memInfo.mem_unit;

	cerr << physMemUsed << endl;


//	{
//		vector<UUmap> maps;
//		unsigned int count = 0;
//		while(count < 2000){
//			UUmap m;
//			initMap(m,(U32)0,(U32)30000);
//			maps.push_back(m);
//			count++;
////			cerr << count << " " << getValue()/1000 << "MB" << endl;
//		}
//		cerr << count << " " << getValue()/1000 << "MB" << endl;
//	}
//	cerr << getValue() << " KB" << endl;


	{
		vector<UUmap > maps;
		unsigned int count = 0;
		while(count < 100){
			unordered_map<U32, U32> m;
			initMap(m,(U32)0,(U32)30000);

			for(U32 i = 0; i < 30000; i++){
				m[i] = i;
			}
			maps.push_back(m);
			count++;
			cerr << count << " " << getValue()/1000 << "MB" << endl;
		}
		cerr << count << " " << getValue()/1000 << "MB" << endl;
	}
	cerr << getValue() << " KB" << endl;
}

int main(int argc, char *argv[]) {
	cerr << "Starting tests..." << endl;


//	test_operators();
//
//	test_ParameterHandler();
//
//	test_Run();
//
//	test_Downloader();
//
//	test_RNAread();
//
	test_Aligner();
//
//	test_myRandomEngine();
//
//	test_ChromosomeInitializer();

//	test_Controller();

//	test_Estimators();
//
//	test_Simulator();
//
//	test_Functions();

//	test_DirichletDistribution();

//	test_DirichletMixture();

//	test_ClusterEstimator();

//	test_MapSize();

//	RunCrashTest();

	cerr << "...finished tests with " << error_count << " errors." << endl;
	return EXIT_SUCCESS;
}
