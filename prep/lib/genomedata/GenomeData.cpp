// ***************************************************************************
// GenomeData.cpp (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <cassert>
#include <sys/stat.h>

#include "GenomeData.h"
#include "Config.h"
#include "MyDefine.h"
#include "Matrix.h"
#include "Fasta.h"
#include "split.h"

#include "api/BamReader.h"
#include "api/BamWriter.h"

extern "C" {
	#include "common.h"
	#include "bbiFile.h"
	#include "bigWig.h"
}

using namespace std;
using namespace BamTools;

void GenomeData::fetchInputs() {
	string bamFile = config.getStringPara("bam");
	struct stat buffer;
	if(stat(bamFile.c_str(), &buffer) == 0 && S_ISREG(buffer.st_mode)) { // merged BAM file
		cerr << "fetch read counts from the merged BAM file " << bamFile << endl;
		fetchInputsFromBAM();
	}
	else if(stat(bamFile.c_str(), &buffer) == 0 && S_ISDIR(buffer.st_mode)) {
		cerr << "fetch read counts from BAM files in directory " << bamFile << endl;
		fetchInputsFromBAMs();
	}
	else {
		cerr << "BAM file or directory containing BAM files is not correctly specified!" << endl;
		return;
	}
}

void GenomeData::fetchInputsFromBAM() {
	int i, j, k;
	long spos, epos;

	//open genome reference file
	string refFile = config.getStringPara("ref");
	FastaReference fr;
    fr.open(refFile, false);
	FastaIndexEntry fie;
	string chr_prefix_fr, ref = fr.index->sequenceNames[0];
	if(ref.find("chr") != string::npos) {
		chr_prefix_fr = "chr";
	}
	else {
		chr_prefix_fr = "";
	}

	//open mappability file
	string mapFile = config.getStringPara("map");
	struct bbiFile *bwf = bigWigFileOpen((char *)mapFile.c_str());
	struct bbiChromInfo *info = bbiChromList(bwf);
	struct bigWigValsOnChrom *chromVals;
	double *chromData;
	string chr_prefix_bw;
	ref = info->name;	
	if(ref.find("chr") != string::npos) {
		chr_prefix_bw = "chr";
	}
	else {
		chr_prefix_bw = "";
	}
	bbiChromInfoFreeList(&info);
	
	string bamFile = config.getStringPara("bam");
	BamReader reader;
	if(!reader.Open(bamFile)) {
	    cerr << "cannot open BAM file " << bamFile << endl;
	    exit(-1);
	}
	if(!reader.LocateIndex()) {
		cerr << "Building BAM index for file " << bamFile << "..." << endl;
		if(!reader.CreateIndex()) {
			cerr << "Index creation failed!" << endl;
			exit(-1);
		}
	}
	vector<RefData> refData_all = reader.GetReferenceData();
	vector<RefData> chromosomes;
	string chrom_list = config.getStringPara("chrlist");
	if(!chrom_list.empty()) {
		bool quit = 0;
		vector<string> chroms = split(chrom_list, ',');
		for(i = 0; i < chroms.size(); i++) {
			int refID = reader.GetReferenceID(chroms[i]);
			if(refID < 0) {
				cerr << "Unrecognized chromosome name: " << chroms[i] << endl;
				cerr << "Use --list for list of valid chromosome names." << endl;
				exit(-1);
			}
			else {
				assert(chroms[i].compare(refData_all[refID].RefName) == 0);
				chromosomes.push_back(refData_all[refID]);
			}
		}
	}
	else {
		chromosomes = refData_all;
	}
	
	//open barcode file
	string barcodeFile = config.getStringPara("barcode");
	ifstream ifs;
	ifs.open(barcodeFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open barcode file " << barcodeFile << endl;
		exit(-1);
	}
	string line, barcode;
	vector<string> barcodes;
	while(getline(ifs, line)) {
		barcode = trim(line);
		barcodes.push_back(barcode);
	}
	ifs.close();
	
	//open output file
	string outputFile = config.getStringPara("output");
    ofstream ofs;
	ofs.open(outputFile.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << outputFile << endl;
		exit(-1);
	}
	//bin size
	int binsize = config.getIntPara("binsize");
	ofs << "#bin size = " << binsize << endl;
	ofs << "chr,gc-content,mappability,";
	for(i = 0; i < barcodes.size()-1; i++) {
		ofs << barcodes[i] << ",";
	}
	ofs << barcodes[barcodes.size()-1] << endl;
	
	//read quality threshold
	int mapQ_th = config.getIntPara("mapQ_th");
	double gc, mapscore;
	long length;
	
	for(i = 0; i < chromosomes.size(); i++) {
		string refName = chromosomes[i].RefName;
		long refLength = chromosomes[i].RefLength;
		int refID = reader.GetReferenceID(refName);
		long bin = binsize;
		
		cerr << "processing chromosome " << refName << "..." << endl;
		
		string chr = abbrOfChr(refName);
		
		if(binsize > refLength) {
			cerr << "Bin size greater than chromosome length of " << refName << ", adjusting to chromosome length: " << refLength << endl;
			bin = refLength;
		}
		
		chromVals = bigWigValsOnChromNew();
		ref = chr_prefix_bw + chr;
		if(bigWigValsOnChromFetchData(chromVals, (char *) ref.c_str(), bwf)){
			chromData = chromVals->valBuf;
		}
		else{
			cerr << "Error: cannot fetch mappability data on chromosome " << refName << endl;
			exit(-1);
		}
		
		ref = chr_prefix_fr + chr;
		
		spos = 0; epos = bin;
		if(reader.SetRegion(refID, 0, refID, refLength)) {
			vector<long> counts(barcodes.size(), 0);
			BamAlignment read;
			while(reader.GetNextAlignmentCore(read)) {
				while(read.Position > epos) {
					if(epos <= refLength) {
						string sequence = fr.getSubSequence(ref, spos, epos-spos);
						calculateGC(sequence, gc, length);
						calculateMapScore(chromData, spos, epos-1, mapscore);
						ofs << refName << "," << gc << "," << mapscore << ",";
						for(j = 0; j < counts.size()-1; j++) {
							ofs << counts[j] << ",";
						}
						ofs << counts[counts.size()-1] << endl;
					}
					
					counts.assign(barcodes.size(), 0);
					spos += bin; 
					epos = spos+bin;
				}
				if(read.MapQuality >= mapQ_th && !read.IsDuplicate()) {
					read.BuildCharData();
					if(read.GetTag("CB", barcode)) {
						for(j = 0; j < barcodes.size(); j++) {
							if(barcode.compare(barcodes[j]) == 0) {
								counts[j] += 1;
								break;
							}
						}
					}
				}
			}
			while(epos <= refLength) {
				string sequence = fr.getSubSequence(ref, spos, epos-spos);
				calculateGC(sequence, gc, length);
				calculateMapScore(chromData, spos, epos-1, mapscore);
				ofs << refName << "," << gc << "," << mapscore << ",";
				for(j = 0; j < counts.size()-1; j++) {
					ofs << counts[j] << ",";
				}
				ofs << counts[counts.size()-1] << endl;
				
				counts.assign(barcodes.size(), 0);
				spos += bin; 
				epos = spos+bin;
			}
			/*
			string sequence = fr.getSubSequence(ref, spos, epos-spos);
			calculateGC(sequence, gc, length);
			calculateMapScore(chromData, spos, epos-1, mapscore);
			ofs << refName << "," << gc << "," << mapscore << ",";
			for(j = 0; j < counts.size()-1; j++) {
				ofs << counts[j] << ",";
			}
			ofs << counts[counts.size()-1] << endl;
			*/
		}
		else {
			cerr << "Error: cannot retrieve chromosome " << refName << " from BAM file." << endl;
			exit(-1);
		}
		bigWigValsOnChromFree(&chromVals);
	}
	
	reader.Close();
	ofs.close();
}

void GenomeData::fetchInputsFromBAMs() {
	int i, j, k, n;
	long spos, epos;

	//open genome reference file
	string refFile = config.getStringPara("ref");
	FastaReference fr;
    fr.open(refFile, false);
	FastaIndexEntry fie;
	string chr_prefix_fr, ref = fr.index->sequenceNames[0];
	if(ref.find("chr") != string::npos) {
		chr_prefix_fr = "chr";
	}
	else {
		chr_prefix_fr = "";
	}

	//open mappability file
	string mapFile = config.getStringPara("map");
	struct bbiFile *bwf = bigWigFileOpen((char *)mapFile.c_str());
	struct bbiChromInfo *info = bbiChromList(bwf);
	struct bigWigValsOnChrom *chromVals;
	double *chromData;
	string chr_prefix_bw;
	ref = info->name;	
	if(ref.find("chr") != string::npos) {
		chr_prefix_bw = "chr";
	}
	else {
		chr_prefix_bw = "";
	}
	bbiChromInfoFreeList(&info);
	
	string bamDir = config.getStringPara("bam");
	
	//open barcode file
	string barcodeFile = config.getStringPara("barcode");
	ifstream ifs;
	ifs.open(barcodeFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open barcode file " << barcodeFile << endl;
		exit(-1);
	}
	string line, barcode;
	vector<string> barcodes;
	while(getline(ifs, line)) {
		barcode = trim(line);
		barcodes.push_back(barcode);
	}
	ifs.close();
	
	//read quality threshold
	int mapQ_th = config.getIntPara("mapQ_th");
	int binsize = config.getIntPara("binsize");
	double gc, mapscore;
	long length;
	
	string chrom_list = config.getStringPara("chrlist");
	
	// calculate GC-content and mappability
	string outputFile = config.getStringPara("output");
    ofstream ofs;
	/*
	string tmpFile = outputFile + ".gcmap";
	ofs.open(tmpFile.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << tmpFile << endl;
		exit(-1);
	}
	*/
	vector<string> chroms = split(chrom_list, ',');
	vector<string> chrs;
	vector<double> gcs, maps;
	for(i = 0; i < chroms.size(); i++) {
		long bin = binsize;
		string chr = abbrOfChr(chroms[i]);
		
		ref = chr_prefix_fr + chr;
		fie = fr.index->entry(ref);
		long refLength = fie.length;
		
		if(bin > refLength) {
			cerr << "Bin size greater than chromosome length of " << chroms[i] << ", adjusting to chromosome length: " << refLength << endl;
			bin = refLength;
		}
		
		chromVals = bigWigValsOnChromNew();
		ref = chr_prefix_bw + chr;
		if(bigWigValsOnChromFetchData(chromVals, (char *) ref.c_str(), bwf)){
			chromData = chromVals->valBuf;
		}
		else{
			cerr << "Error: cannot fetch mappability data on chromosome " << chr << endl;
			exit(-1);
		}
		
		ref = chr_prefix_fr + chr;
		
		spos = 0;
		while(spos+bin <= refLength) {
			string sequence = fr.getSubSequence(ref, spos, bin);
			calculateGC(sequence, gc, length);
			calculateMapScore(chromData, spos, spos+bin-1, mapscore);
			//ofs << chr << "," << gc << "," << mapscore << endl;
			chrs.push_back(chr); gcs.push_back(gc); maps.push_back(mapscore); 
			spos += bin; 
		}
		bigWigValsOnChromFree(&chromVals);
	}
	ofs.close();
	
	long binCount = chrs.size();
	
	// count reads
	/*
	tmpFile = outputFile + ".readcounts";
	ofs.open(tmpFile.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << tmpFile << endl;
		exit(-1);
	}
	*/
	
	int num_cell = barcodes.size();
	Matrix<int> counts(num_cell, binCount, 0);
	for(n = 0; n < num_cell; n++) {
		string bamFile = bamDir + "/" + barcodes[n] + ".bam";
		
		cerr << "processing BAM " << barcodes[n] << "..." << endl;
		
		BamReader reader;
		if(!reader.Open(bamFile)) {
			cerr << "cannot open BAM file " << bamFile << endl;
			exit(-1);
		}
		if(!reader.LocateIndex()) {
			cerr << "Building BAM index for file " << bamFile << "..." << endl;
			if(!reader.CreateIndex()) {
				cerr << "Index creation failed!" << endl;
				exit(-1);
			}
		}
		
		vector<RefData> refData_all = reader.GetReferenceData();
		vector<RefData> chromosomes;
		for(i = 0; i < chroms.size(); i++) {
			int refID = reader.GetReferenceID(chroms[i]);
			if(refID < 0) {
				cerr << "Unrecognized chromosome name: " << chroms[i] << endl;
				exit(-1);
			}
			else {
				assert(chroms[i].compare(refData_all[refID].RefName) == 0);
				chromosomes.push_back(refData_all[refID]);
			}
		}
		
		k = 0;
		for(i = 0; i < chromosomes.size(); i++) {
			string refName = chromosomes[i].RefName;
			long refLength = chromosomes[i].RefLength;
			int refID = reader.GetReferenceID(refName);
			long bin = binsize;
			
			cerr << "\tchromosome " << refName << "..." << endl;
			
			string chr = abbrOfChr(refName);
			
			if(bin > refLength) {
				cerr << "Bin size greater than chromosome length of " << refName << ", adjusting to chromosome length: " << refLength << endl;
				bin = refLength;
			}
			
			spos = 0; epos = bin;
			if(reader.SetRegion(refID, 0, refID, refLength)) {
				long count = 0;
				BamAlignment read;
				while(reader.GetNextAlignmentCore(read)) {
					while(read.Position > epos) {
						/*
						if(epos+bin <= refLength || i < chromosomes.size()-1) {
							ofs << count << ",";
						}
						else {
							ofs << count << endl;
						}
						*/
						counts[n*binCount+k] = count;
						k++;
						count = 0;
						spos += bin; 
						epos = spos+bin;
					}
					if(read.MapQuality >= mapQ_th && !read.IsDuplicate()) {
						count++;
					}
				}
				while(epos <= refLength) {
					/*
					if(epos+bin <= refLength || i < chromosomes.size()-1) {
						ofs << count << ",";
					}
					else {
						ofs << count << endl;
					}
					*/
					counts[n*binCount+k] = count;
					k++;
					count = 0;
					spos += bin; 
					epos = spos+bin;
				}
			}
			else {
				cerr << "Error: cannot retrieve chromosome " << refName << " from BAM file " << bamFile << endl;
				exit(-1);
			}
		}
		reader.Close();
	}
	//ofs.close();
	
	//merge results into a single file
	ofs.open(outputFile.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << outputFile << endl;
		exit(-1);
	}
	// header
	ofs << "#bin size = " << binsize << endl;
	ofs << "chr,gc-content,mappability,";
	for(i = 0; i < barcodes.size()-1; i++) {
		ofs << barcodes[i] << ",";
	}
	ofs << barcodes[barcodes.size()-1] << endl;
	for(i = 0; i < binCount; i++) {
		ofs << chrs[i] << "," << gcs[i] << "," << maps[i] << ",";
		for(n = 0; n < num_cell-1; n++) {
			ofs << counts[n*binCount+i] << ",";
		}
		ofs << counts[(num_cell-1)*binCount+i] << endl;
	}
	
	ofs.close();
}

void GenomeData::calculateGC(string sequence, double &gc_content, long &length) {
	long gc_count = 0;
	long n_count = 0;
	char c;
	length = 0;

	for(int i = 0; i < sequence.length(); i++){
		c = sequence.at(i);
		if(c == '>'){
			break;
		}
		switch(c) {
			case 'G':
			case 'C':
			case 'g':
			case 'c':
				gc_count++;
				length++;
				break;
			case 'A':
			case 'T':
			case 'a':
			case 't':
				length++;
				break;
			default:
				n_count++;
				length++;
				break;
		}
	}

	gc_content = (n_count < length)? (double) gc_count/(length-n_count) : -1;
}

void GenomeData::calculateMapScore(double * chromData, long long spos, long long epos, double &mapScore) {
	mapScore = 0;
	for(int i = spos; i <= epos; i++){
		mapScore += chromData[i];
	}
	mapScore = mapScore/(epos-spos+1);
}


