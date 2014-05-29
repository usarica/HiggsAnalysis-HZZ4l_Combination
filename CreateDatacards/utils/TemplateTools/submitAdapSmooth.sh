#!/bin/bash
LISTDIR=$1
ERGTEV=$2"TeV"
TEMPLATESLISTFILE=$LISTDIR"/List_of_Templates.txt"
TEMPLATESDIR=$LISTDIR"/builder_"$ERGTEV"/"
while read i;
do
	FILE=$TEMPLATESDIR$i
	bsub -q 8nh -o "lsflog_"$i"_"$ERGTEV".txt" -e "lsferr_"$i"_"$ERGTEV".err" produceAdapSmooth3D.lsf.sh $FILE
done < $TEMPLATESLISTFILE
