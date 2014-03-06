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
        self.is_checked = False

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

    def query_point(self, query, is_sorted):
        """
        Generator function. Yield a list of intervals (iv) in the tree that
        contain the query point.
        """

        yield [iv for iv in self.centered if iv[0] <= query <= iv[1]]

        if query < self.center \
        and self.children[0] \
        and not self.children[0].is_checked:
 
            is_empty = True
            for found in self.children[0].query_point(query):
                if found:
                    is_empty = False
                yield found

            if is_empty and is_sorted:
                self.children[0].is_checked = True

        elif self.center < query \
        and self.children[1] \
        and not self.children[1].is_checked:

            is_empty = True
            for found in self.children[1].query_point(query):
                if found:
                    is_empty = False
                yield found

            if is_empty and is_sorted:
                self.children[1].is_checked = True

    def query_interval(self, query, is_sorted):
        """
        Generator funtion. Yield a list of intervals (iv) in the tree that
        overlap, at least partially, with the query interval.
        """

        yield [iv for iv in self.centered if not (iv[1] < query[0] \
                                               or query[1] < iv[0])]

        if query[0] < self.center \
        and self.children[0] \
        and not self.children[0].is_checked:

            is_empty = True
            for found in self.children[0].query_interval(query):
                if found:
                    is_empty = False
                yield found

            if is_empty and is_sorted:
                self.children[0].is_checked = True

        if self.center < query[1] \
        and self.children[1] \
        and not self.children[1].is_checked:

            is_empty = True
            for found in self.children[1].query_interval(query):
                if found:
                    is_empty = False
                yield found

            if is_empty and is_sorted:
                self.children[1].is_checked = True


class InputFormatError(Exception):
    """
    Raise this exception if input format is not BED.
    """

    pass

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
    create children nodes with the remainig intervals. 
    """

    set_recursion_limit()

    # Centers are sorted to find the median (center).
    intervals.sort()
    # s: start, e: end, i: info.
    try:
        centers = [(s+e)/2.0 for (s,e,i) in intervals]
    except ValueError:
        centers = [(s+e)/2.0 for (s,e) in intervals]
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
        slines = (line.rstrip().split(None,3) for line in fin)
        for inter in slines:
            try:
                interval = [int(inter[1]),int(inter[2])]
            
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
        chromosomes[chromosome].sort()
        trees[chromosome] = build_tree(chromosomes[chromosome])

    del chromosomes
 
    return trees

def DFS_interval_query(trees, query, is_sorted):
    """
    Generator funtion. Perform a depth-first search of the query interval on the
    interval tree and yield the tree intervals that overlap with the query. 
    query must be a 3-tuple like (chromosome, start, end) for intervals and a
    2-tuple like (chromosome, position) for points.
    """

    set_recursion_limit()

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
            for found in trees[chromosome].query_interval(position, is_sorted):
                total_found.extend(found)
        else:
            for found in trees[chromosome].query_point(position, is_sorted):
                total_found.extend(found)
    
    if not total_found:
        return
    for hit in total_found:
        yield tuple([chrom] + [field for field in hit])

def query_from_main(trees, interval, is_sorted):
    try:
        query = [sl[0], int(sl[1])]
    except IndexError:
        raise InputFormatError('Query file must be in BED format.')
    except ValueError:
        raise InputFormatError('Query file must be in BED format.')
    
    if len(query) > 2:
        try:
            query.append(int(sl[2]))
        except ValueError:
            raise InputFormatError(
                'Query file must be in BED format.'
                )
        
    query = tuple(query)
    for overlap in DFS_interval_query(trees, query, is_sorted):
        print overlap


if __name__ == '__main__':
    try:
        trees = create_trees_dict(sys.argv[1])
        query_file = sys.argv[2]
        if sys.argv[3] == '--sorted':
            query_is_sorted = True
        elif sys.argv[3] == '--unsorted':
            query_is_sorted = False
        else:
            raise IndexError # Jump to the bottom line.
    
        with open(query_file) as quf:
            slines = (line.rstrip().split(None,3) for line in quf)
            for sl in slines:
                query_from_main(trees, sl, query_is_sorted)

    except IndexError:
        print ('Please use the following format:\n' +
               'ting.py <tree_file> <query_file> --(sorted|unsorted)')
