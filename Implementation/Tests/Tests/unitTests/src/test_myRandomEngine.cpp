#include "../../../../src/myRandomEngine.cpp"

#include "../../../catch/catch.hpp"
#include <iostream>
#include "../../../../headers/Operators.h"

using namespace std;

TEST_CASE( "Random", "" )
{
	myRandomEngine r;
	std::vector<int> v;

	v.push_back(190);

	REQUIRE(r.select_randomly(v) == 190);
}

TEST_CASE( "selectRandomly", "" )
{
	myRandomEngine r;

	r.seed(12);
	REQUIRE(r.select_randomly(200) == 30);
}

TEST_CASE( "shuffleExceptLast", "" )
{
	myRandomEngine r;

	std::vector<int> v;

    for(int i = 0; i < 20; i++){
	    v.push_back(i);
    }

	r.seed(167);
    r.shuffleExceptLast(v);

    cout << v << endl;
    cout << "size : " << v.size() << endl;

    REQUIRE(v.size() == 20);
	REQUIRE(v[19] == 19);

    r.shuffleExceptLast(v);
	REQUIRE(v[19] == 19);

    r.shuffleExceptLast(v);
	REQUIRE(v[19] == 19);
}


