#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
event_dist.py -- Tip to event distance calculation

This script retrieves a distances from a phylogenetic tree tip to an event
on the same tree.

Requirements:

Written by Moisès Bernabeu <moises.bernabeu.sci@gmail.com>
Barcelona, 2022
'''


# Importing libraries ----
from optparse import OptionParser
from rooted_phylomes import ROOTED_PHYLOMES as root_dict
import ete3
import pandas as pd
from multiprocessing import Process, Manager
from treefuns import get_species, root, annotate_tree, \
    tree_stats, get_group_mrca, count_dupl_specs
from utils import file_exists, create_folder


# Definitions ----
def get_ndists(tree, phylome_id, rootdict, gnmdf, spcol,
               normcol, normtag, evcol, evtag):
    treel = tree.split('\t')
    print('Calculating: ', treel[0])
    t = ete3.PhyloTree(treel[3], sp_naming_function=get_species)

    root(t, rootdict)
    t.get_descendant_evol_events()

    annotate_tree(t, gnmdf, spcol, [normcol, evcol])

    norm_group = get_group_mrca(t, treel[0], normcol, normtag)

    norm_stats = tree_stats(norm_group['node'])
    nsd = count_dupl_specs(norm_group['node'])
    nfactor = norm_stats['median']

    evdict = get_group_mrca(t, treel[0], evcol, evtag, treel[0])

    odict = dict()
    odict['seed'] = treel[0]
    odict['species'] = get_species(treel[0])
    odict['event_dist'] = evdict['node'].get_distance(treel[0])
    odict['event_ndist'] = evdict['node'].get_distance(treel[0]) / nfactor
    odict['seed_dist'] = t.get_distance(treel[0])
    odict['seed_ndist'] = t.get_distance(treel[0]) / nfactor
    odict['nfactor'] = nfactor

    odict = {**odict,
             **{'norm_' + k: v for k, v in norm_stats.items()},
             **{'norm_' + k: v for k, v in nsd.items()},
             **{'whole_' + k: v for k, v in tree_stats(t).items()},
             **{'whole_' + k: v for k, v in count_dupl_specs(t).items()}}

    return odict


class dist_process(Process):
    def __init__(self, tree_row, phylome_id, rootdict, gnmdf,
                 spcol, normcol, normtag, evcol, evtag, olist):
        Process.__init__(self)
        self.tree_row = tree_row
        self.phylome_id = phylome_id
        self.rootdict = rootdict
        self.gnmdf = gnmdf
        self.spcol = spcol
        self.normcol = normcol
        self.normtag = normtag
        self.evcol = evcol
        self.evtag = evtag
        self.olist = olist

    def run(self):
        odict = get_ndists(self.tree_row, self.phylome_id, self.rootdict,
                           self.gnmdf, self.spcol, self.normcol, self.normtag,
                           self.evcol, self.evtag)
        self.olist.append(odict)


def main():
    # Script options definition ----
    parser = OptionParser()
    parser.add_option('-d', '--def', dest='default',
                      help='Default configurations to execute with our data.',
                      action='store_true')
    parser.add_option('-f', '--file', dest='ifile',
                      help='In file',
                      metavar='<path/to/file.txt>')
    parser.add_option('-g', '--groups', dest='groups',
                      help='Groups file (.csv)',
                      metavar='<path/to/file.csv>')
    parser.add_option('-o', '--out', dest='output',
                      help='output directory',
                      metavar='<path/to/folder>')
    parser.add_option('-c', '--cpu', dest='cpus',
                      help='Number of CPUs', type='int',
                      metavar='<N>')
    (options, args) = parser.parse_args()

    if options.default:
        infile = '../test/0076_108.txt'
        gnmdffile = '../test/0076_norm_groups.csv'
        outdir = '../outputs/'
        cpus = 4
    else:
        infile = options.ifile
        gnmdffile = options.groups
        outdir = options.output
        cpus = options.cpus

    ofilenm = infile.rsplit('/', 1)[1].split('.', 1)[0]
    ofile = '%s/%s_dist.csv' % (outdir, ofilenm)

    if not file_exists(ofile):
        create_folder(outdir)

        gnmdf = pd.read_csv(gnmdffile)
        phylome_id = infile.rsplit('/', 1)[1].split('_', 1)[0]

        with Manager() as manager:
            olist = manager.list()

            processes = list()
            for tree_row in open(infile, 'r'):
                if tree_row != '':
                    if len(processes) >= cpus:
                        done = False
                        while not done:
                            for process in processes:
                                if not process.is_alive():
                                    processes.remove(process)
                                    done = True

                    process = dist_process(tree_row, phylome_id,
                                           root_dict[int(phylome_id)], gnmdf,
                                           'Proteome', 'Normalising group',
                                           'A', 'Metazoan', 'metazoan', olist)
                    processes.append(process)
                    process.start()

            for process in processes:
                process.join()

            # Writing output files
            odf = pd.DataFrame(list(olist))
            odf.to_csv(ofile, index=False)

    return 0


if __name__ == '__main__':
    main()
