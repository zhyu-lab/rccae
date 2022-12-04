// ***************************************************************************
// Config.h (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _CONFIG_H
#define _CONFIG_H

#include <cstdlib>
#include <string>
#include <map>

using namespace std;

class Config {
	private:
		map<string, string> stringParas;
		map<string, long> intParas;
		map<string, double> realParas;
		
	public:
		Config();
		~Config() {}
		
		bool isVerbose() {return !(intParas["verbose"] == 0);}
		
		string getStringPara(string paraName);
		void setStringPara(string paraName, string value);
		long getIntPara(string paraName);
		void setIntPara(string paraName, long value);
		double getRealPara(string paraName);
		void setRealPara(string paraName, double value);
		
};


#endif

