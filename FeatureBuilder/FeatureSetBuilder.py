''' 
Glycan feature set builder
Parses a set of glycan sub-tree files (i.e., output of TreeBuilder.py) and returns a 0/1 feature matrix

Input : A directory path containg a list of .subtrees files 
Output : A feature matrix printed as a csv. Each row is a glycan and columns are feature vector with 1 or 0 indicating observance of that feature vector
'''

import argparse
import sys
import pdb
import re
import pandas as pd
from os import listdir
from os.path import isfile, join



# ---- Main ---- #
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--DIRECTORY', required=True, help='Input filedirectory (required)')
	args = parser.parse_args()

	# read in list of files
	path = args.DIRECTORY
	subtree_files = [f for f in listdir(path) if f.endswith('.subtrees')] #(join(path,

	# parse each file into a dictionary['subtree']['glycan'] = number
	trees = list()
	
	subtree_dictionary = {}
	glycans = list()
	for  filename in subtree_files:
		filename = join(path, filename)
#		pdb.set_trace()
		with open(filename) as f:
			# grab glycan id
			glycan_id = f.readline().rstrip()[1:]
			glycans.append(glycan_id)
			for line in f:
				line = line.rstrip()
				if line not in subtree_dictionary:
					subtree_dictionary[line]={}
					(subtree_dictionary[line])[glycan_id] = 1
				else:
					if glycan_id in subtree_dictionary[line]:
						subtree_dictionary[line][glycan_id] += 1
					else:
						(subtree_dictionary[line])[glycan_id] = 1
	
		print "Read file %s " % filename

#	subtrees =subtree_dictionary.keys()
#	subtrees.sort(key=lambda s:len(s))
#	glycans.sort(key=int)	
	print "Done reading files\n"

	# Change to a dictionary and print
	df = pd.DataFrame(subtree_dictionary)
	cols = df.columns.tolist()
	cols.sort(key=lambda s:len(s))
#	pdb.set_trace()
	df = df[cols]
	df = df.fillna(0)
#	df.reindex_axis(sorted(df.columns, key=lambda s:len(s)), axis=1)
	df.to_csv(join(path, 'Feature_matrix.csv'))

	
#	pdb.set_trace()
#	# Building the feature matrix
#	fo = open(join(path, 'Feature_matrix.csv'), "wb")
#	writer = csv.writer(fo)
#	pdb.set_trace()
#	column_headers = subtrees
#	column_headers.insert(0, 'glycan')
#	writer.writerows(column_headers)
#	for g in glycans:
#		list_feature = list(g)
#		for s in subtrees:
#			list_feature.append(subtree_dictionary[s][g])		
#		pdb.set_trace()
#		writer.writerows(list_feature)
#	fo.close()
#
#
	print "Done!"
#
