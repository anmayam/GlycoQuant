#!/usr/bin/env python

''' 
Glycan Feature Builder
Parses a list of glycan structures that are based on LINUCS encoding and returns a set of sub-tree based features for each glycan

Input : A glyde XML file containing a set of glycans 
Output : A features text 1/0 file indicating presence/absence of a feature 
'''

import argparse
import sys
import pdb
import re

try:
	import xml.etree.cElementTree as ET #choose faster method if available
except ImportError:
	import xml.etree.ElementTree as ET

	
def is_sugar(s):
	sugars = ['man', 'glcnac', 'gal', 'galnac', 'neuac', 'neugc', 'dehex']
	if s.lower() in sugars:
		return True;

def add_derivative(d):
	derivatives = ['n-acetyl']
	if d.lower() in derivatives:
		return 'nac'

# ---- Main ---- #
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--INPUT', required=True, help='Input file(required)')
	parser.add_argument('-o', '--OUTPUT', required=True, help='Output filename (required)')
	
	args = parser.parse_args()


	# read in glyde file
	gtree = ET.ElementTree(file=args.INPUT)

	# start parsing
	for nglycan in gtree.iter(tag='molecule'):
		nglycan_id = nglycan.get('id')
		subtrees = list()
		
		# get all mono residues first. 
		base_type =''
		residue_numbers = list()	
		for elem in nglycan.iter('residue'):			
			residue_ref = elem.get('ref')		
			residue_type = elem.get('subtype')			
			if (residue_type == 'base_type'):
				splits = re.split('=|-', residue_ref);			
				base_type = splits[2]+splits[3][1:]
				residue_numbers.append(elem.get('partid'))
			elif (residue_type =='substituent'):
				splits = re.split('=', residue_ref)
				sub_type = splits[1]
				base_type = base_type+add_derivative(sub_type)
			else:
				raise Exception('Unrecognizable text in residue subtype in Glyde format')
			if(is_sugar(base_type[1:])):
				subtrees.append(base_type)

			
		# getting all the other trees		
		for start_node in residue_numbers:
			pdb.set_trace()
			r = start_node
			tree_string = r
			while (not r):
				for elem in nglycan.iterfind('residue_link[@to="'+r+'"]'):
					from_residue = elem.get('from')
					if from_residue in residue_numbers:
						for child in elem:
							link = child.get('from')+'-'+child.get('to')
							tree_string = from_residue+link+tree_string
							subtree.append(tree_string)
				r = from_residue
								
	subtrees.append(splits[2]+splits[3][1:]) # should be form bMan, aglcnac etc. 
				
	print "Done!\n"

			


			

	
