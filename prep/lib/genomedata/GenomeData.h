// ***************************************************************************
// GenomeData.h (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _GENOMEDATA_H
#define _GENOMEDATA_H

#include <vector>
#include <map>
#include <cstring>

#include "Matrix.h"

class GenomeData {
	private:
		void calculateGC(string sequence, double &gc_content, long &length);
		void calculateMapScore(double *chromData, long long spos, long long epos, double &mapScore);
		
		void fetchInputsFromBAM();
		void fetchInputsFromBAMs();
		
	public:
		GenomeData() {}
		~GenomeData() {}
		
		void fetchInputs();
};


#endif

