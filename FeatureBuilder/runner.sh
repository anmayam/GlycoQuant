#!/bin/bash

DATAFolder="/media/sf_VBoxSharedFolder/glycan_quant/data/col_cancer/GlycansList"
sort_unique=1
for file in /$DATAFolder/*.glyde
do
	out="$file.subtrees"
	python TreeBuilder.py -i $file -o $out -u $sort_unique
done
