#!/usr/bin/env python
#-*- coding:utf-8 -*-

# TING:
# The module contains a Node class with methods for adding nodes and querying
# intervals, and functions for building and traversing the tree of nodes.
# 
# Generate an interval tree from an input file containing intervals:
# Eeach line in the input file must contain just one interval.
# Each line/interval must follow this format: Chromosome, Start, End,
# Information [and any optional extra Information columns at the end], separated
# by whitespace (/[\s\t]+/).
# Example:
# Chromosome    Start   End     Info                    [Extra Info]
#                               (e.g. chromatin color)  (e.g. mappability value)
# chr1          456675  456994  red                     172 
# Information in the Extra Info column/s is treated as a single column and must
# be explicitly unpacked if needed.

import sys
from collections import defaultdict

# Really big trees could cause RuntimeError: maximum recursion depth exceeded.
# Change this value to an integer > 1000 if you get this error.
MAX_RECURSION_DEPTH = sys.getrecursionlimit()

class Node():
    """
    Create Node instances that store the center position of the node, a list of
    centered intervals overlapping this center, and up to two children nodes.
    """
    def __init__(self, center):
        # If center has not a proper value, implicitly raise a ValueError.
        self.center = int(center)
        self.centered = []
        self.children = [None, None]

    def __repr__(self):
        return 'Node at %d. Overlaps with %d intervals. %d children nodes' \
                % (self.center,
                   len(self.centered),
                   len([1 for child in self.children if child is not None])
                   )

    def add_child(self, new_child, index):
        """
        Hang a new_child Node instance from the current (parent) node.
        Index is 0 or 1 if new_child is left or right child, respectively.
        """
        assert self.children[index] is None, 'Child overwriting!'
        self.children[index] = new_child

    def query_point(self, query):
        """
        Generator function. Yield a list of intervals (iv) in the tree that
        contain the query point.
        """
        yield [iv for iv in self.centered if iv[0] <= query <= iv[1]]

        if query < self.center and self.children[0]:
            for found in self.children[0].query_point(query):
                yield found

        elif self.center < query and self.children[1]:
            for found in self.children[1].query_point(query):
                yield found

    def query_interval(self, query):
        """
        Generator funtion. Yield a list of intervals (iv) in the tree that
        overlap, at least partially, with the query interval.
        """
        yield [iv for iv in self.centered if not (iv[1] < query[0] \
                                               or query[1] < iv[0])]

        if query[0] < self.center and self.children[0]:
            for found in self.children[0].query_interval(query):
                yield found

        if self.center < query[1] and self.children[1]:
            for found in self.children[1].query_interval(query):
                yield found

def set_recursion_limit():
    """
    Leave maximum recursion depth at its default value (1000) if it is not
    specified or is set at a value < 1000, change it to the specified value
    otherwise.
    """
    if MAX_RECURSION_DEPTH and MAX_RECURSION_DEPTH > 1000:
        sys.setrecursionlimit(MAX_RECURSION_DEPTH)

def build_tree(intervals):
    """
    Create a Node instance and initialize its centered list, which contains all
    the intervals that overlap with the center position of the node. Recursively
    create children nodes containing the remainig intervals. 
    """
    set_recursion_limit()

    # Centers have been previously sorted for calculating the median.
    # s: start, e: end, c: color, i: info.
    try:
        centers = [(s+e)/2.0 for (s,e,c,i) in intervals]
    except ValueError:
        centers = [(s+e)/2.0 for (s,e,c) in intervals]
    center = centers[len(centers)//2]
                
    centered = [i for i in intervals if i[0] <= center <= i[1]]
    left = [i for i in intervals if i[1] < center]
    right = [i for i in intervals if center < i[0]]

    node = Node(center)
    node.centered = centered
    if left:
        node.add_child(build_tree(left), 0)
    if right:
        node.add_child(build_tree(right), 1)
    
    return node

def create_trees_dict(intervals_file):
    """
    Create a dictionary in which each key is a chromosome and its value is an
    interval tree containing all the intervals in that chromosome.
    """
    # Temp dictionary that holds a list of all the intervals in each chromosme.
    chromosomes = defaultdict(list)
    with open(intervals_file) as fin:
        sline = (line.rstrip().split(None,4) for line in fin)
        for inter in sline:
            if len(inter) == 5:
                interval = (int(inter[1]),int(inter[2]),inter[3],inter[4])
            elif len(inter) == 4:
                interval = (int(inter[1]),int(inter[2]),inter[3])
            else:
                raise ValueError('Tree file must have the following format:\n\
                                 chromosome\tstart\tend\tinfo\t[other_info]')
            chromosomes[inter[0]].append(interval)

    # Final dictionary with one tree per chromosome.
    trees = defaultdict()
    for chromosome in chromosomes:
        chromosomes[chromosome].sort()
        trees[chromosome] = build_tree(chromosomes[chromosome])

    del chromosomes
 
    return trees

def DFS_point_query(trees, query):
    """
    Perform a depth-first search of the query point on the interval tree and
    yield the tree intervals that overlap with the query.
    query must be a 2-tuple like (chromosome, position).
    """
    pass

def DFS_interval_query(trees, query):
    """
    Generator funtion. Perform a depth-first search of the query interval on the
    interval tree and yield the tree intervals that overlap with the query. 
    query must be a 3-tuple like (chromosome, start, end).
    """
    set_recursion_limit()

    chrom, start, end = query # Implicitly raises a ValueError if not 3-tuple.
    position = (int(start), int(end))
    total_found = []
    for chromosome in trees:
        if chrom != chromosome:
            continue
        for found in trees[chromosome].query_interval(position):
            total_found.extend(found)
    
    if not total_found:
        return
    for hit in total_found:
        if len(hit) == 5:
            yield (chrom, hit[0], hit[1], hit[2], hit[3])
        elif len(hit) == 4:
            yield (chrom, hit[0], hit[1], hit[2])
    
if __name__ == '__main__':
    try:
        trees = create_trees_dict(sys.argv[1])
        query_file = sys.argv[2]
    
        with open(query_file) as quf:
            slines = (line.rstrip().split() for line in quf)
            queries = ((sl[0], int(sl[1]), int(sl[2])) for sl in slines)

            for query in queries:
                for overlap in  DFS_interval_query(trees, query):
                    print overlap
    
    except IndexError:
        print ('Please use the following format:\n' +
               'ting.py <tree_file> <query_file>')
