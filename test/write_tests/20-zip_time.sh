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

zip -r $OUTPUT_FOLDER/${name}.zip $INPUT_FILE

