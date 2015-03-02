#!/usr/bin/env python

''' 
Glycan Feature Builder
Parses a list of glycan structures that are based on LINUCS encoding and returns a set of sub-tree based features for each glycan

Input : A tab-delimited file with two columns :glycan index and LINUCS string
Output : A features 1/0 file indicating presence/absence of a feature
'''

import argparse
import sys
import pandas as pd
from pandas import DataFrame
import pdb
import re


def is_sugar(s):
	sugars = ['Man', 'GlcNAc', 'Gal', 'GalNAc', 'NeuAc', 'NeuGc', 'DeHex']
	if s[1:] in sugars:
		return True;
	#s_index = [i for i, word in enumerate(sugars) if word.endswith(s)]
	#if (s_index > -1):
	#	return True; 
		

# ---- Structure parser ----- #
def parse_structure(g_linucs):
	# get mono structures which will be found in []
	mono = list()
	splits = re.split(r'[\],\[,{}]', g_linucs)
	for s in splits:				
		if(len(s)>1 and is_sugar(s) and s not in mono):		
			mono.append(s)
			print mono
	# get disaccharide structures
	# get rid of the core
	core = '[][bGlcNAc]{[(4+1)][bGlcNAc]{[(4+1)]'
	match = re.sub(core, '',g_linucs)
	gminuscore = matchl
	pdb.set_trace()
	distructures = list()
	splits = re.split('{', gminuscore)
	parent = splits[0]
	for s in splits[1,]:
		if (splits[i].startswith('}')):
			parent = splits[i][1,]
		else:
			di = parent+splits[i]
			print di
			pdb.set_trace()
			


 


	
	
	

# ---- Main ---- #
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--INPUT', required=True, help='Input file(required)')
	parser.add_argument('-o', '--OUTPUT', required=True, help='Output filename (required)')
	
	args = parser.parse_args()


	# read in glycan structure file
#	pdb.set_trace()
	gData=pd.read_table(args.INPUT, sep='\t', header=0, index_col=0)
	print("Read structure data :%s" % args.INPUT)


	# start parsing
	for index, row in gData.iterrows():
		subtrees = parse_structure(row['LINUCS']) 
	
