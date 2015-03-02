#!/usr/bin/env python

''' 
Glycan Subtree Builder
Parses a list of glycan structures in GLYDE format and build subtrees for each glycan

Input : A glyde XML file containing a set of glycans 
Output : A text file indicating list of all sub-trees 
'''

import argparse
import sys
import pdb
import re
import operator

from collections import deque, OrderedDict

try:
	import xml.etree.cElementTree as ET #choose faster method if available
except ImportError:
	import xml.etree.ElementTree as ET


def is_sugar(s):
	''' Build terms for monosaccharide recognition here
	gro : indicates Neu5Ac, see a-dgro-dgal-NON-2:6|1:a|2:keto|3:d||(5d:1)n-acetyl in MonosaccharidesDB.org'''
	
	sugars = ['man', 'glcnac', 'gal', 'galnac', 'neunac', 'neugc', 'dehex'] 
	if s.lower() in sugars:
		return True;

def add_derivative(d):
	''' adds derivative to base type '''
	derivatives = ['n-acetyl']
	if d.lower() in derivatives:
		return 'nac'


def get_neighbors(elements):
 	''' get neighbors of an iterable'''
	''' to do : incorporate this for GalNAc containing sugars, since GalNac will currently split into two nodes'''
	iterator = iter(elements)
	previous = None
	item = iterator.next()
	for next in iterator:
		yield (prev, item, next)
		prev = item
		item - next
	yield (prev, item, None)

def get_carb_id(s):
	''' returns simple names for sugars in monosaccharideDB form '''
	carb_names = {'x-dglc-HEX-1:5': 'xglc',				
		       	'b-dglc-HEX-1:5' : 'bglc',
                       	'b-dman-HEX-1:5' : 'bman',
			'a-dman-HEX-1:5' : 'aman',
			'b-dgal-HEX-1:5' : 'bgal',				
			'a-dgro-dgal-NON-2:6|1:a|2:keto|3:d' :'aneu',
			'a-lgal-HEX-1:5|6:d' : 'adehex'}
	if (s in carb_names.keys()):
		return(carb_names[s])
	else:
		print s
		raise Exception('Unrecognized base type monosaccharide')

def make_residue_link(Gtree, res1, res2, link_type):
	'''
	Generates a link between two sugar residues 
	link_type indicates  direction from node 1 to node 2, 
	forward	i.e from reducing to non-reducing end of the glycan.
	reverse from nonreducing to reducing
	The reverse edge is also created and added to Gtree
	'''
	if res1 not in Gtree:
		Gtree[res1] ={ }
	(Gtree[res1])[res2] = 'forward-'+link_type 
	if res2 not in Gtree:
		Gtree[res2] = { }
	(Gtree[res2])[res1] = 'reverse-'+link_type[::-1]
	return Gtree

def glycan_bfs(Gtree, start_res, end_res, names_of_res):
	''' Returns shortest path between two nodes.
	Uses both forward and reverse edges and original saccharide names
	'''
	visited_res = set()	
	first = list()
	first.append("_"+start_res) 
	queue = deque()
	queue.append(first) 	
	while queue:
		path = queue.pop()
		res = re.split("_", path[-1])[-1]# path[-1][-1]
		if res in visited_res:			
			continue
		visited_res.add(res)
		if res == end_res:
			final_path = path[-1]				
			visited_res = sorted(visited_res, key=int, reverse=True)
			for res in visited_res: #names_of_res.keys():
				final_path = re.sub(r'_'+res, names_of_res[res], final_path)
			return final_path
		for neighbor in Gtree[res]:		
			new_path = list()			
			new_path.append(path[-1] +"_"+Gtree[res][neighbor]+"_"+neighbor)
			queue.append(new_path)	
				

# ---- Main ---- #
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--INPUT', required=True, help='Input file(required)')
	parser.add_argument('-o', '--OUTPUT', required=True, help='Output filename (required)')
	parser.add_argument('-u', '--UNIQUE', required =False, help ='returns unique subtrees (0/1)')
	
	args = parser.parse_args()


	# params
	print_only_unique=False
	if (args.UNIQUE):
		print_only_unique = int(args.UNIQUE)
	
	# read in glyde file
	g_glyde = ET.ElementTree(file=args.INPUT)

	# set output
	fo = open(args.OUTPUT, "w")

	# start parsing'14': 'reverse-6
	residue_dict = {}
	for nglycan in g_glyde.iter(tag='molecule'):
		subtrees = list()
		fo.write(">"+nglycan.get('id')+"\n")
		
		# get all mono residues first. 
		base_type =''
		residue_numbers = list()	
		all_elements = nglycan.iter('residue')
		for elem in all_elements:			
			residue_ref = elem.get('ref')		
			residue_type = elem.get('subtype')			
			if (residue_type == 'base_type'):
				splits = re.split('=', residue_ref);			
				base_type = splits[1]
				base_type = get_carb_id(base_type)			
				residue_numbers.append(elem.get('partid'))			
			elif (residue_type =='substituent'):				
				splits = re.split('=', residue_ref)
				sub_type = splits[1]
				base_type = base_type+add_derivative(sub_type)		
			else:
				raise Exception('Unrecognizable type in Glyde format')
		
			if(is_sugar(base_type[1:])):			
				subtrees.append(base_type)
				residue_dict[residue_numbers[-1]] = base_type

		# build a tree for all nodes using a dictionary
		pdb.set_trace()
		residue_numbers.sort()
		g_tree ={}
		for r in residue_numbers:
			for elem in nglycan.iterfind('residue_link[@to="'+r+'"]'):
				from_residue = elem.get('from')				
				if from_residue in residue_numbers:
					for child in elem:
						link = child.get('from')+'-'+child.get('to')
				 		make_residue_link(g_tree, from_residue, r, link)
						
		# Return all sub-trees using a BFS		
		pdb.set_trace()
		residues = sorted(g_tree.keys(), reverse=False)		
		residues.sort(key=int, reverse=True)		
		for i in range(0, len(residues)-1):
			res1 = residues[i]			
			for j in range(i+1, len(residues)):
				res2 = residues[j]
				res_path  = glycan_bfs(g_tree,res1, res2, residue_dict)
				print "res  = " + res1 + " to " + res2  + ":" + res_path
				subtrees.append(res_path)
		
 
		#get uniques, sort and print
		if (print_only_unique):
			subtrees = list(set(subtrees))
		subtrees.sort(key = lambda s:len(s))
		for stree in subtrees:			
			fo.write(stree +"\n")

	
	fo.close()			
	print "Done!\n"

			


			

	
