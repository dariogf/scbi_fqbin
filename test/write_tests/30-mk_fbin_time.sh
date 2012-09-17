#!/bin/bash

if [ "$#" != 3 ]; 
then
        echo "Usage: $0 INPUT_FILE OUTPUT_FOLDER RESULTS_FILE"
        echo "==============================================================="
        exit -1
fi
INPUT_FILE=$1;
OUTPUT_FOLDER=$2;
RESULTS_FILE=$3;


name=`basename $INPUT_FILE`

fasta=`head -c 1 $INPUT_FILE |grep '>'`

if [[ -z $fasta ]]; then
    mk_fbin -o $OUTPUT_FOLDER/${name}.fbin $INPUT_FILE
else
    mk_fbin -o $OUTPUT_FOLDER/${name}.fbin -F $INPUT_FILE
fi

