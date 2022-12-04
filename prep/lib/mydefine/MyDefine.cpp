// ***************************************************************************
// MyDefine.cpp (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <cstdio>
#include <fstream>
#include <vector>
#include <cmath>

#include "MyDefine.h"

using namespace std;

MyDefine::MyDefine() {
}

MyDefine::~MyDefine() {
}

/*** definition of global vars ***/
Config config;
InputParser parser;
GenomeData genomedata;
/*** end of definition ***/

/*** definition of general functions ***/
//****** trim string ******//
string trim(const string &str, const char *charlist) {
	string ret(str);
	size_t indx = ret.find_first_not_of(charlist);
	if(indx != string::npos) {
		ret.erase(0, indx);
		indx = ret.find_last_not_of(charlist);
		ret.erase(indx+1);
	}
	else {
		ret.erase();
	}
	return ret;
}

//****** abbreviation of chromosome name ******//
string abbrOfChr(string chr) {
	string abbr_chr = chr;
	size_t i = abbr_chr.find("chrom");
	if(i == string::npos) {
		i = abbr_chr.find("chr");
		if(i != string::npos) {
			abbr_chr = abbr_chr.substr(i+3,abbr_chr.size()-3);
		}
	}
	else {
		abbr_chr = abbr_chr.substr(i+5,abbr_chr.size()-5);
	}
	return abbr_chr;
}
/*** end of definition ***/

