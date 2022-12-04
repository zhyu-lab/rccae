#!/bin/bash

echo 

function usage() {
	echo "  $0 <bam> <ref> <map> <barcode> <output>"
	echo "Options:"
	echo "  bam: a merged BAM file (10x Genomics) or a directory containing BAM files of the cells to be analyzed"
	echo "  ref: genome reference file(.fasta)"
	echo "  map: mappability file(.bw)"
	echo "  barcode: a file listing the barcodes of all cells (for 10x Genomics) or names of all BAM files to be analyzed (one BAM per cell)"
	echo "  output: output directory to save results"
	echo "Examples:"
	echo "  $0 merged.bam hg19.fa hg19.bw barcodes.txt results"
	echo "  $0 bam_dir hg19.fa hg19.bw bams.txt results"
}

if [ $# -ne 5 ]; then
	echo "error inputs!"
	echo -n "Usage: " 
	usage
	exit 0
fi

bam=$1
ref=$2
map=$3
barcodes=$4
output_dir=$5

log_file=$output_dir/log.txt

# chromosomes to be analyzed
chroms=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22
# chroms=1,2,3

test -e $log_file && rm $log_file

test ! -e $output_dir && mkdir -p $output_dir

current=`date "+%Y-%m-%d %H:%M:%S"`
seconds_s=`date -d "$current" +%s`

echo "Step 1. get read counts, GC-content and mappability......"
./prep/bin/prepInput -b $bam -r $ref -m $map -B $barcodes -c $chroms -s 500000 -q 20 -o $output_dir/readcounts.txt > $log_file 2>&1

echo "Step 2. detect tumor clones and denoise read counts data......"
python ./cae/train.py --input $output_dir/readcounts.txt --epochs 100 --batch_size 64 --lr 0.0001 --latent_dim 3 --seed 0 --output $output_dir >> $log_file 2>&1

echo "Step 3. detect single-cell copy number alterations......"
# the path to MCR v91 needs to be specified
./hmm/run_SCHMM.sh /usr/local/MATLAB/MATLAB_Runtime/v91 $output_dir/lrc.txt $output_dir 10 >> $log_file 2>&1

current=`date "+%Y-%m-%d %H:%M:%S"`
seconds_e=`date -d "$current" +%s`
let secs=seconds_e-seconds_s
echo "Elapsed time: $secs seconds!"

exit 0
