// ***************************************************************************
// MyDefine.h (c) 2022 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _MYDEFINE_H
#define _MYDEFINE_H

#include <iostream>
#include <cstring>
#include <cstdlib>

#include "Config.h"
#include "InputParser.h"
#include "GenomeData.h"

class MyDefine {
	public:
		MyDefine();
		~MyDefine();
};

/*** declaration of global vars ***/
extern Config config;
extern InputParser parser;
extern GenomeData genomedata;
/*** end of declaration ***/

/*** declaration of general functions ***/
string trim(const string &str, const char *charlist = " \t\r\n");
string abbrOfChr(string chr);
/*** end of declaration ***/


#endif

