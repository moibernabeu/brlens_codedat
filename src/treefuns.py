#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
treefuns.py -- Functions useful to manipulate phylome trees

Requirements:
 - operator
 - scipy
 - numpy

Written by Moisès Bernabeu <moigil.bernabeu.sci@gmail.com>
April 2022
'''

# Import libraries ----
from operator import itemgetter
from scipy import stats
import numpy as np


# Define functions ----
def get_species(node):
    '''
    Get species name

    Args:
        node (TreeNode): tree node with sequence name

    Returns:
        string: the tree node species label
    '''

    if '_' in node:
        return node.split("_")[1]
    else:
        return node


def root(tree, root_dict):
    '''
    Root the tree according to a rooting dictionary

    The tree is rooted with a dictionary containing species-to-age information,
    the farthest sequence from an species in the tree which has maximum age is
    selected to be the outgroup of the tree.

    Args:
        tree (PhyloTree): ete3 PhyloTree object with a get_species_tag function

        root_dict (dictionary): a dictionary containing the species age,
        indexes refers to the phylome.
        Eg.: 3: {'SP1': 1, 'SP2': 2}
        Where 3 is the phylome number and the SP2 is older than SP1.

    Returns:
        string: the outgroup sequence and the tree object is rooted without
        being returned.

    Raises:
        Exception: description
    '''

    # Checking whether any species is in the tree
    if any(sp in root_dict for sp in tree.get_species()):
        # Getting the rooting subdictionary with the species in the tree
        ogdval = max([root_dict.get(sp, 0) for sp in tree.get_species()])

        # Getting the outgroup species
        ogsps = [k for k, val in root_dict.items()
                 if val == ogdval and k in tree.get_species()][0]

        # Getting the outgroup sequences dictionary
        ogseqdict = {seq: tree.get_distance(seq) for seq in
                     tree.get_leaf_names() if ogsps in seq}

        # Getting the farthest oldest leaf
        ogseq = max(ogseqdict)
    else:
        # Setting the root in the farthest leaf
        ogseq = tree.get_farthest_leaf()[0].get_leaf_names()[0]

    tree.set_outgroup(ogseq)

    return ogseq


def annotate_tree(tree, df, spcol, cols):
    '''
    Tree leaves annotation

    Tree leaves are annotated with the columns

    Args:
        tree (PhyloTree): ete3 PhyloTree object with a get_species_tag function

        df (DataFrame): pandas dataframe containing the information to
        annotate the tree

        spcol (string): the column containing the species names

        cols (string or list of strings): column or columns containing the
        annotations

    Returns:
        string: the input tree is annotated, the funtion returns 0

    Raises:
        Exception: description
    '''

    # Converting string to single element list
    if type(cols) == 'str':
        cols = [cols]

    # Iterating the leaves
    for leaf in tree.get_leaves():
        # Iterating the dataframe columns
        for col in cols:
            # Annotating the leaf with its column feature
            sp = list(leaf.get_species())[0]
            feat_val = list(df[col][df[spcol] == sp])[0]
            leaf.add_feature(col, feat_val)

    return 0


def tree_stats(tree):
    '''
    Get tree branch stats

    From a tree the function retrieves basic numerical information about the
    branches lengths. First, it creates a list where stores all the root to
    tip distances, then calculates the median, the mean, the width and the sum

    Args:
        tree (PhyloTree): ete3 PhyloTree object with a get_species_tag function

    Returns:
        dictionary: dictionary with the tree statistics

    Raises:
        Exception: description
    '''

    # Getting all the leaves
    ndlf = tree.get_leaves()

    # Retrieving the root to tip distances
    distl = list()
    for leaf in ndlf:
        distl.append(tree.get_distance(leaf))

    # Generating the output dictionary
    nodedict = dict()
    nodedict['leafno'] = len(distl)
    nodedict['median'] = np.median(distl)
    nodedict['mean'] = np.mean(distl)
    nodedict['width'] = tree.get_farthest_leaf()[1]
    nodedict['sum'] = sum(distl)
    nodedict['kurt'] = stats.kurtosis(distl)
    nodedict['skew'] = stats.skew(distl)

    return nodedict


def get_group_mrca(tree, tree_id, feature, value, sp_in=None):
    '''
    Get the greatest monophyletic subtree of a labeled group

    The function retrieves a subtree which all leaves have the value of value
    of the feature indicated. Check all groups with these conditions and
    returns one that maximizes the number of leaves.

    Args:
        tree (PhyloTree): ete3 PhyloTree object with a get_species_tag function

        tree_id (string): the phylome tree idea

        feature (string): the name of the feature where the MRCA label
        is stored

        value (string): the name of the value of the feature that defines
        the clade

        sp_in (string): the name of the species that has to be inside the
        MRCA group

    Returns:
        dictionary: dictionary with some information about the node and
        the node

    Raises:
        Exception: description
    '''

    # Getting the number of leaves in the tree
    tlno = len(tree.get_leaf_names())

    # Getting the monophyletic sets containing all species with col attribute
    mphysets = list()
    for st in tree.traverse():
        # Get the tree leaves
        leaves = st.get_leaves()
        stlno = len(leaves)

        # Width of the monophyletic group
        stwdth = st.get_farthest_leaf()[1]

        # Generating list of col attributes in subtree
        feat_list = list()
        for stleag in st.get_leaves():
            feat_list.append(getattr(stleag, feature))

        # Checking the subtree only contains 1 'col' attributes, it has more
        # than 1 leaves, it is not the entire tree and it is not a polytomy
        lnames = st.get_leaf_names()
        if sp_in is None:
            sptoincl = lnames[0]
        else:
            sptoincl = sp_in

        if (sptoincl in lnames and len(set(feat_list)) == 1 and stlno > 1 and
                stlno != tlno and stwdth != 0):
            # Appending to a list a dictionary with the basic information of
            # the group monophyletic group
            mphy = dict()
            mphy['tree'] = tree_id
            mphy['node'] = st
            mphy['seq_no'] = stlno
            mphy[feature] = feat_list[0]
            mphysets.append((stlno, mphy))

    # Getting the maximum leaves monophyletic group for each normalising group
    mphylist = sorted(mphysets, key=itemgetter(0), reverse=True)

    return mphylist[0][1]


def count_dupl_specs(tree):
    '''
    The function iterates over the nodes and counts the number of speciation
    and duplication events

    Args:
        tree (PhyloTree): ete3 phylogenetic tree annotated with the
                          evolutionary events

    Returns:
        dict: duplication and speciation counts dictionary
    '''

    # Creating an empty dictionary
    cdict = {'D': 0, 'S': 0}

    # Iterating through nodes and increasing the counter in the dict according
    # to the event
    for node in tree.iter_search_nodes():
        if len(node.get_leaf_names()) > 1:
            cdict[node.evoltype] += 1

    return cdict
