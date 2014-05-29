#!/bin/bash

cd ${LS_SUBCWD}

echo "LSF job running in: " `pwd`

eval `scram runtime -sh`

echo $CMSSW_VERSION

FILENAME=$1".json"

./buildTemplate.exe $FILENAME
