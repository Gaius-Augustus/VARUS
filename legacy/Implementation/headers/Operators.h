/*
 * Operators.h
 *
 *  Created on: 25.11.2016
 *      Author: willy
 */

#ifndef OPERATORS_H_
#define OPERATORS_H_

#include <unordered_map>
#include <ostream>
#include <iostream>
#include <sstream>
#include <memory>

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "debug.h"
#include "TypeConventions.h"
#include "Run.h"

#include <algorithm>    // std::random_shuffle
#include <cstdlib>
#include <ctime>        // std::time
#include "math.h"

/* ! \brief: This File holds a few operators and functions for things that are needed a lot.
 * e.g. passing a map to an ostream.
 */


template<typename K, typename V>
std::ostream& operator<<(std::ostream& os, std::unordered_map<K, V>& map) {
    for (const auto& p : map) os << p.first << ":" << p.second << ";";
    return os;
}

template<typename K, typename V>
std::stringstream& operator<<(std::stringstream& os, std::unordered_map<K, V>& map) {
    for (const auto& p : map) os << p.first << ":" << p.second << ";";
    return os;
}

template<template <typename...> class Map, typename K, typename V, typename V2>
void operator += (Map<K, V> &m1, const Map<K, V2> &m2) {
	for (const auto& p : m2) m1[p.first] += p.second;
}

template<template <typename...> class Map, typename K, typename V, typename V2>
void operator *= (Map<K, V> &m, const V2 scalar) {
	for (const auto& p : m) m[p.first] *= scalar;
}

template<template <typename...> class Map, typename K, typename V, typename V2>
void operator /= (Map<K, V> &m, const V2 scalar) {
	for (const auto& p : m) m[p.first] /= scalar;
}

template<template <typename...> class Map, typename K, typename V1, typename V2>
Map<K, V1> operator + (Map<K, V1> &m1, const Map<K, V2> &m2) {
//	Map<K, V1> ret;
	Map<K, V1> ret(m1);
//	for (const auto& p : m1) ret[p.first] = p.second;
	for (const auto& p : m2) ret[p.first] += p.second;
	return ret;
}

template<template <typename...> class Map, typename K, typename V, typename V2>
Map<K,V> operator / (const Map<K,V> &m, V2 n) {
	Map<K, V> ret;
	for (const auto& p : m) ret[p.first] = (V)p.second/n;
	return ret;
}

template<template <typename...> class Map, typename K, typename V>
void initMap(Map<K,V> &m, const V a, const K n) {
	for(K i = 0; i < n; i++) {
		m[i] = a;
	}
}

template<template <typename...> class Map, typename V>
double euklidNorm(Map<U32,V> &m) {
	double ret = 0.0;
	for(const auto& p : m){
		ret += p.second*p.second;
	}
	ret = sqrt(ret);
	return ret;
}

template<template <typename...> class Map, typename K, typename V>
V sum(const Map<K,V> &m) {
	V sum = 0;
	for (const auto& p : m) sum += (V) p.second;
	return sum;
}

//template<template <typename...> class Map, typename V, typename V2>
//void scaleMap(Map<U32,V> &m,  const V2 scaleTo) {
//	double scalar = scaleTo/sum(m);
//	m *= scalar;
//}

void scaleMap(UDmap &m,  const double scaleTo);

template<template <typename...> class Map, typename K, typename V>
void set0(Map<K,V> &m) {
	for (const auto& p : m) m[p.first] = 0;
}

template<template <typename...> class Map, typename K, typename V>
void dotProduct(Map<K,V> &m1) {
//	Map<K,V> ret;
	for (const auto& p : m1) m1[p.first] = p.second*p.second;
//	return m1;
}

template<template <typename...> class Map, typename K, typename V>
std::string toStr(Map<K,V> &m) {
	std::stringstream ss;
	std::unordered_map<unsigned int, double>::iterator i;
	for (const auto& p : m) ss << p.second << ":";
	return ss.str();
}

template <typename T>
unsigned int num_of_digits(T i)
{
    // returns the number of digits of an integer
    unsigned int digits = 1;
    while(i /= 10)
    {
        digits++;
    }
    return digits;
}

//template<typename K, typename V>
//std::vector<V> map2vec(std::unordered_map<K, V> &m) {
//	std::vector<double> v;
//	v.resize(m.size());
//    for (const auto& p : m) {
//		if(p.first >= 0 && p.first < v.size()) {
//			v[p.first] = p.second;
//		}
////		else {
////			DEBUG(1, "WARNING: Wrong index for vector!");
////		}
//	}
//	return v;
//}

template<typename K, typename V>
std::vector<V> map2vec(std::unordered_map<K, V> &m) {
	std::vector<V> v;
    for (const auto& p : m) {
    	v.push_back(p.second);
	}
	return v;
}

//std::ostream& operator<<(std::ostream& os, std::unique_ptr<Run> &r) {
//	os << r->accesionId;
//}


template<typename V>
std::ostream& operator<<(std::ostream& os, std::vector<V>& v) {
    for(unsigned int i = 0; i < v.size(); i++) {
    	os << i << ":" << v[i] << ";";
    }
    return os;
}

//std::ostream& operator<<(std::ostream& os, std::vector<std::unique_ptr<Run>>& v) {
//    for(unsigned int i = 0; i < v.size(); i++) {
//    	os << i << ":" << v[i] << ";";
//    }
//    return os;
//}




//template<template <typename...> class Map, typename V>
//Map<unsigned int,V> operator = (Map<unsigned int,V> &m, std::vector<V>& v) {
//    for(unsigned int i = 0; i < v.size(); i++) {
//    	m[i] = v[i];
//    }
//}

template<typename K, typename V>
void normMap(std::unordered_map<K,V> &m) {
	double s = sum(m);
    for (const auto& p : m) {
		m[p.first] = p.second/s;
	}
}

template<typename V>
V sum(std::vector<V> &v) {
	V sum = 0;
    for(unsigned int i = 0; i < v.size(); i++) {
		sum += v[i];
	}
	return sum;
}

template<typename V>
void normVector(std::vector<V> &v) {
	double s = sum(v);
    for(unsigned int i = 0; i < v.size(); i++) {
		v[i] = v[i]/s;
	}
}

//unsigned int select_randomly(unsigned int max);

//template<typename V>
//V select_randomly(std::vector<V> &vec){
//
//	return vec[select_randomly(vec.size()-1)];
//}



void remove_substr(const std::string &sub, std::string &str);

#endif /* OPERATORS_H_ */
