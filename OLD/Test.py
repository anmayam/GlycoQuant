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

# graph is in adjacent list representation
graph = {
        '1': ['2', '3', '4'],
        '2': ['5', '6'],
        '5': ['9', '10'],
        '4': ['7', '8'],
        '7': ['11', '12']
        }

def bfs(graph, start, end):
    # maintain a queue of paths
    queue = []
    # push the first path into the queue
    queue.append([start])
    while queue:
	#pdb.set_trace()
        # get the first path from the queue
        path = queue.pop(0)
        # get the last node from the path
        node = path[-1]
        # path found
        if node == end:
            return path
        # enumerate all adjacent nodes, construct a new path and push it into the queue
        for adjacent in graph.get(node, []):
            new_path = list(path)
            new_path.append(adjacent)
            queue.append(new_path)

print bfs(graph, '5', '6')
