// ***************************************************************************
// InputParser.cpp (c) 2019 zhenhua yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <limits.h>

#include "InputParser.h"
#include "MyDefine.h"

using namespace std;

void InputParser::parseArgs(int argc, char *argv[]) {
	string refFile = "", bwFile = "", chrlist = "";
	string bamFile = "", barcodeFile = "", outFile = "";
	int binsize = 500000, mapQ_th = 20;

	struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"bam", required_argument, 0, 'b'},
		{"ref", required_argument, 0, 'r'},
		{"map", required_argument, 0, 'm'},
		{"barcode", required_argument, 0, 'B'},
		{"chrlist", required_argument, 0, 'c'},
		{"binsize", required_argument, 0, 's'},
		{"mapQ", required_argument, 0, 'q'},
		{"output", required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};

	int c;
	//Parse command line parameters
	while((c = getopt_long(argc, argv, "hb:r:m:B:c:s:q:o:", long_options, NULL)) != -1){
		switch(c){
			case 'h':
				usage(argv[0]);
				exit(0);
			case 'b':
				bamFile = optarg;
				break;
			case 'r':
				refFile = optarg;
				break;
			case 'm':
				bwFile = optarg;
				break;
			case 'B':
				barcodeFile = optarg;
				break;
			case 'c':
				chrlist = optarg;
				break;
			case 's':
				binsize = atoi(optarg);
				break;
			case 'q':
				mapQ_th = atoi(optarg);
				break;
			case 'o':
				outFile = optarg;
				break;
			default :
				usage(argv[0]);
				exit(1);
		}
	}
	
	if(bamFile.empty()){
        cerr << "Use --bam to specify a merged BAM file (.bam) or a directory containing BAM files of all cells." << endl;
		usage(argv[0]);
        exit(1);
    }

	if(refFile.empty()){
		cerr << "Use --ref to specify reference file (.fasta)." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(bwFile.empty()){
		cerr << "Use --map to specify mappability file(.bigwig)." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(barcodeFile.empty()){
		cerr << "Use --barcode to specify a file that lists the barcodes of all cells or names of all BAM files to be analyzed." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(chrlist.empty()){
		cerr << "Use --chrlist to specify the entries chromosomes to be analyzed (should be separated by commas)." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(outFile.empty()){
		cerr << "Use --output to specify the output file." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(binsize < 100000) {
		cerr << "Error: the value of parameter \"binsize\" should be at least 100000." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(mapQ_th < 0) {
		cerr << "Error: the value of parameter \"mapQ\" should be a nonnegative integer!" << endl;
		usage(argv[0]);
		exit(1);
	}
	
	config.setStringPara("ref", refFile);
	config.setStringPara("map", bwFile);
	config.setStringPara("bam", bamFile);
	config.setStringPara("barcode", barcodeFile);
	config.setStringPara("chrlist", chrlist);
	config.setStringPara("output", outFile);
	config.setIntPara("binsize", binsize);
	config.setIntPara("mapQ_th", mapQ_th);
}

void InputParser::usage(const char* app) {
	cerr << "Usage: " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -b, --bam <string>              a merged BAM file (10x Genomics) or a directory containing BAM files of the cells to be analyzed" << endl
		<< "    -r, --ref <string>              genome reference file(.fasta)" << endl
		<< "    -m, --map <string>              mappability file(.bw)" << endl
		<< "    -B, --barcode <string>          a file listing the barcodes of all cells or names of all BAM files to be analyzed" << endl
		<< "    -c, --chrlist <string>          the entries chromosomes to be analyzed (should be separated by commas)" << endl
		<< "    -s, --binsize <int>             set the size of bin to count reads [default:500000]" << endl
		<< "    -q, --mapQ <int>                threshold value for mapping quality [default:20]" << endl
		<< "    -o, --output <string>           output file" << endl
		<< endl
		<< "Example:" << endl
		<< "   " << app << " -b example.bam -r hg19.fasta -m hg19.bw -B barcodes.txt -c 1,2,3 -o example.out" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}
