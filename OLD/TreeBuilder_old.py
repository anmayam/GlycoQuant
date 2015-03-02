#!/usr/bin/env python

''' 
Glycan Feature Builder
Parses a list of glycan sub-trees  for each glycan

Input : A glyde XML file containing a set of glycans 
Output : A text file indicating list of subtrees for each glycan 
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
	sugars = ['man', 'glcnac', 'gal', 'galnac', 'neuac', 'neugc', 'dehex']
	if s.lower() in sugars:
		return True;

def add_derivative(d):
	derivatives = ['n-acetyl']
	if d.lower() in derivatives:
		return 'nac'

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
	(Gtree[res1])[res2] = 'forward_'+link_type 
	if res2 not in Gtree:
		Gtree[res2] = { }
	(Gtree[res2])[res1] = 'reverse_'+link_type
	return Gtree

def glycan_dfs(Gtree, res, has_traversed):
	has_traversed[res] = True
	chain = "res =" + res ; 
	for neighbour in Gtree[res]:
		#print neighbour
		#pdb.set_trace()	
		if neighbour not in has_traversed:
			glycan_dfs(Gtree, neighbour, has_traversed)
		chain = chain + "res= " + neighbour + ", edge= " + Gtree[res][neighbour]
		return chain
		#pdb.set_trace()

def glycan_bfs(Gtree, start_res, end_res):	
	visited_res = set()
	queue = deque([start_res]) 
	chain = start_res
	while queue:
		res = queue.pop()
		pdb.set_trace()
		if res in visited_res:
			continue
		visited_res.add(res)
		if res == end_res:
			return chain
		for neighbor in Gtree[res]:
			if neighbor not in visited_res:
				chain= chain + "_" +Gtree[res][neighbor]+"_"+neighbor
				if neighbor == end_res:
					return chain
				else:
					queue.appendleft(neighbor)	
				

		


# ---- Main ---- #
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--INPUT', required=True, help='Input file(required)')
	parser.add_argument('-o', '--OUTPUT', required=True, help='Output filename (required)')
	
	args = parser.parse_args()


	# read in glyde file
	g_glyde = ET.ElementTree(file=args.INPUT)

	# start parsing
	for nglycan in g_glyde.iter(tag='molecule'):
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

		
		# build a tree for all nodes using a dictionary
		residue_numbers.sort()
		g_tree ={}
		for r in residue_numbers:
			for elem in nglycan.iterfind('residue_link[@to="'+r+'"]'):
				from_residue = elem.get('from')				
				if from_residue in residue_numbers:
					for child in elem:
						link = child.get('from')+'-'+child.get('to')
				 		make_residue_link(g_tree, from_residue, r, link)
						
		# Return all trees using a DFS
		#chain = ''
		#has_traversed={}
		#for key, value in sorted(g_tree.iteritems(), key=lambda (k,v): (v,k)):
		#	print "%s %s" % (key,value)

		pdb.set_trace()

		for k in g_tree:
			all_neighbors = g_tree[k]
			all_neighbors_sorted = sorted(all_neighbors.items(), key=operator.itemgetter(0), reverse=False)
			g_tree[k] = OrderedDict(all_neighbors_sorted)		
		residues = sorted(g_tree.keys(), reverse = False)


		
		for i in range(0, len(residues)-1):
			res1 = residues[i]
			pdb.set_trace()
			for j in range(i+1, len(residues)):
				j=4				
				res2 = residues[j]
				print "res  = " + res1 + " to " + res2			 
				print glycan_bfs(g_tree,res1, res2)+"\n"
			pdb.set_trace()




#if res not in has_traversed:
#				fchain = glycan_dfs(g_tree, res, has_traversed)
#				print fchain
								
	 
				
	print "Done!\n"

			


			

	
