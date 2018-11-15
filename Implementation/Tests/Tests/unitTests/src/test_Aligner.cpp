/*
 * test_Aligner.cpp
 *
 *  Created on: 31.12.2017
 *      Author: willy
 */


#include "../../../../headers/Alligner.h"
#include "../../../catch/catch.hpp"
#include <iostream>
#include "../../../../headers/Operators.h"

using namespace std;

ParameterHandler *param = new ParameterHandler();

TEST_CASE( "Aligner 3 Files", "For some runs of Pombe runs are split into 3 files" )
{

	param->outFileNamePrefix = "/home/willy/BRAKER/VARUS/Implementation/Tests/Tests/unitTests/Aligner/";
	param->verbosity = -1;
	Alligner all(param);

	Run r;
	r.accesionId = "testRun";
	r.N = 0;
	r.X = 10;
	r.paired = true;

	string s = all.shellCommand(&r);

	string should = "/home/willy/Bachelorarbeit/STAR-2.5.1b/source/STAR --runThreadN 4 --genomeDir /home/willy/Bachelorarbeit/Manager/genome/ --readFilesIn /home/willy/BRAKER/VARUS/Implementation/Tests/Tests/unitTests/Aligner//testRun/N0X10//bla_1.fasta /home/willy/BRAKER/VARUS/Implementation/Tests/Tests/unitTests/Aligner//testRun/N0X10//bla_3.fasta --outFileNamePrefix /home/willy/BRAKER/VARUS/Implementation/Tests/Tests/unitTests/Aligner//testRun/N0X10/";
//	cout << s << endl;

	REQUIRE(should == s);
}

TEST_CASE( "Aligner 2 Files", " the standard case for paired runs" )
{

	param->outFileNamePrefix = "/home/willy/BRAKER/VARUS/Implementation/Tests/Tests/unitTests/Aligner/";
	param->verbosity = -1;
	Alligner all(param);

	Run r;
	r.accesionId = "testSplit2";
	r.N = 0;
	r.X = 10;
	r.paired = true;

	string s = all.shellCommand(&r);

	string should = "/home/willy/Bachelorarbeit/STAR-2.5.1b/source/STAR --runThreadN 4 --genomeDir /home/willy/Bachelorarbeit/Manager/genome/ --readFilesIn /home/willy/BRAKER/VARUS/Implementation/Tests/Tests/unitTests/Aligner//testSplit2/N0X10//bla_1.fasta /home/willy/BRAKER/VARUS/Implementation/Tests/Tests/unitTests/Aligner//testSplit2/N0X10//bla_2.fasta --outFileNamePrefix /home/willy/BRAKER/VARUS/Implementation/Tests/Tests/unitTests/Aligner//testSplit2/N0X10/";
//	cout << s << endl;

	REQUIRE(should == s);
}


