#!/usr/bin/env python
#-*- coding:utf-8 -*-

# TING: TING Is Not GenomicRanges
# ((Tree|Thousands) of INtervals based on GenomicRanges)
# Generate an interval tree from a BED file containing intervals.
# The module contains a Node class with methods for adding nodes and querying
# intervals, and functions for building and traversing the tree of nodes.
# 
# Input files have to be in BED format:
# http://genome.ucsc.edu/FAQ/FAQformat#format1
# http://www.ensembl.org/info/website/upload/bed.html
# Information in the optional fields is treated as a single column and must
# be explicitly unpacked if needed.

import sys
from collections import defaultdict

# gzopen decompresses the file before calling open() if the file is compressed.
try:
    from gzopen import gzopen as open
except ImportError:
    pass

class Node():
    """
    Create Node instances that store the center position of the node, a list of
    centered intervals overlapping this center, and up to two children nodes.
    """

    def __init__(self, center):
        # If center is not convertible to int, implicitly raise a ValueError.
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

        if query < self.center \
        and self.children[0]: 
 
            for found in self.children[0].query_point(query):
                yield found

        elif self.center < query \
        and self.children[1]:

            for found in self.children[1].query_point(query):
                yield found

    def query_interval(self, query):
        """
        Generator function. Yield a list of intervals (iv) in the tree that
        overlap, at least partially, with the query interval.
        """

        yield [iv for iv in self.centered if iv[0] <= query[1] \
                                         and iv[1] >= query[0]]

        if query[1] < self.center \
        and self.children[0]:

            for found in self.children[0].query_interval(query):
                yield found

        if self.center < query[0] \
        and self.children[1]:

            for found in self.children[1].query_interval(query):
                yield found


class InputFormatError(Exception):
    """
    Raise this exception if input format is not BED.
    """

    pass

def build_tree(intervals):
    """
    Create a Node instance and initialize its centered list, which contains all
    the intervals that overlap with the center position of the node. Recursively
    create children nodes with the remainig intervals. 
    """

    # s: start, e: end, i: info.
    try:
        centers = [(s+e)/2.0 for (s,e,i) in intervals]
    except ValueError:
        centers = [(s+e)/2.0 for (s,e) in intervals]
    center = sorted(centers)[len(centers) // 2]
                
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
        slines = (line.rstrip().split(None, 3) for line in fin)
        for inter in slines:
            try:
                interval = [int(inter[1]), int(inter[2])]
            
            # Check that the input file is in BED format.
            except IndexError:
                raise InputFormatError('Tree file must be in BED format.')
            except ValueError:
                raise InputFormatError('Tree file must be in BED format.')

            if len(inter) > 3:
                interval.append(inter[3])
                interval = tuple(interval)
            
            chromosomes[inter[0]].append(interval)

    # Final dictionary with one tree per chromosome.
    trees = defaultdict()
    for chromosome in chromosomes:
        trees[chromosome] = build_tree(chromosomes[chromosome])

    del chromosomes
 
    return trees

def find_overlaps(trees, query):
    """
    Generator function. Perform a depth-first search of the query interval on
    the interval tree and yield the tree intervals that overlap with the query. 
    query must be a 3-tuple like (chromosome, start, end) for intervals and a
    2-tuple like (chromosome, position) for points.
    """

    if len(query) == 3:
        chrom, start, end = query
        position = (int(start), int(end))
        is_interval = True
    elif len(query) == 2:
        chrom, position = query
        is_interval = False
    else:
        raise InputFormatError('Query file must be in BED format.')

    total_found = []
    for chromosome in trees:
        if chrom != chromosome:
            continue
        if is_interval:
            for found in trees[chromosome].query_interval(position):
                total_found.extend(found)
        else:
            for found in trees[chromosome].query_point(position):
                total_found.extend(found)
    
    if not total_found:
        return
    for hit in total_found:
        yield tuple([chrom] + [field for field in hit])

def _query_from_main(trees, interval):
    """
    This function is defined for readability, it is called from main().
    """
    try:
        query = [interval[0], int(interval[1])]
    except IndexError:
        raise InputFormatError('Query file must be in BED format.')
    except ValueError:
        raise InputFormatError('Query file must be in BED format.')
    
    if len(interval) > 2:
        try:
            query.append(int(interval[2]))
        except ValueError:
            raise InputFormatError('Query file must be in BED format.')
        
    query = tuple(query)
    for overlap in find_overlaps(trees, query):
        print '\t'.join(str(i) for i in overlap)


if __name__ == '__main__':
    try:
        tree_file = sys.argv[1]
        query_file = sys.argv[2]

        trees = create_trees_dict(tree_file)

        with open(query_file) as quf:
            slines = (line.rstrip().split(None, 3) for line in quf)
            for sl in slines:
                _query_from_main(trees, sl)

    except IndexError:
        print ('Please use the following format:\n' +
               'ting.py <tree_file> <query_file>')
