#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
get_trees.py -- Get trees from phylome db

This script runs in parallel threads the download of the phylomes data and
trees. It requires a list in plain text format separated by line breaks.

Requirements: download.py and thread.py

Written by Moisès Bernabeu <mail@mail.com>
January 2022
'''

# Import libraries ----
from optparse import OptionParser
from os.path import isfile
from os import stat
from thread import thread
import time
# import ftplib
from rooted_phylomes import ROOTED_PHYLOMES


# Definitions ----
def yet_downloaded(dir, name):
    '''
    Check the existence and size of the file, returns True or False
    '''
    fp = '/'.join([dir, name])

    # Check file existence
    if isfile(fp):
        # Check the file is not empty
        if stat(fp).st_size != 0:
            return True
    else:
        return False


def main():
    # Script options definition ----
    parser = OptionParser()
    parser.add_option('-d', '--def', dest='default',
                      help='Default configurations to execute with our data.',
                      action='store_true')
    parser.add_option('-a', '--all', dest='all',
                      help='Download all phylomes.',
                      action='store_true')
    parser.add_option('-f', '--file', dest='ifile',
                      help='In file',
                      metavar='<path/to/file.txt>')
    parser.add_option('-w', '--workdirir', dest='workdir',
                      help='Working directory in which data will be stored',
                      metavar='<path/to/workdir>')
    parser.add_option('-t', '--th', dest='threads',
                      help='Number of threads to be used.',
                      metavar='<N>', type='int')
    (options, args) = parser.parse_args()

    # FTP initial direction
    ftp = 'ftp://phylomedb.org/phylomedb/phylomes/phylome_'
    trees = 'best_trees.txt.gz'
    data = 'phylome_info.txt.gz'
    gene = 'all_gene_names.txt.gz'
    prot = 'all_protein_names.txt.gz'

    # Default settings testing
    if options.default:
        workdir = '../outputs'
        threads = 4
        ifile = '../data/phylome_list.txt'
        pdbids = open(ifile)
    elif options.all:
        # pdbftp = ftplib.FTP('phylomedb.org')
        # pdbftp.login()
        # pdbftp.cwd('phylomedb/phylomes/')
        # files = pdbftp.nlst()
        # pdbids = [item.replace('phylome_', '') for item in files]

        pdbids = [str(item).zfill(4) for item in list(ROOTED_PHYLOMES.keys())]
        threads = 4
        workdir = '../outputs'
    else:
        workdir = options.workdir
        threads = options.threads
        ifile = options.ifile
        pdbids = open(ifile)

    # Execution of all threads in the threadpool
    sleeptime = 1
    tasks = list()

    for line in pdbids:
        line = line.replace('\n', '')
        if line != '':
            if (len(tasks) >= threads):
                # Wait for a process to finish
                done = False
                while not done:
                    time.sleep(sleeptime)
                    for task in tasks:
                        if not task.is_alive():
                            # With these conditions the thread is ended, print
                            # and free a slot
                            task.show()
                            tasks.remove(task)
                            done = True

            # Downloading filenames
            tree_url = ftp + line + '/' + trees
            data_url = ftp + line + '/' + data
            gene_url = ftp + line + '/' + gene
            prot_url = ftp + line + '/' + prot

            tree_ofile = line + '_' + trees
            data_ofile = line + '_' + data
            gene_ofile = line + '_' + gene
            prot_ofile = line + '_' + prot

            # Download trees
            if not yet_downloaded(workdir, tree_ofile):
                tree_thread = thread(line, tree_url, workdir, tree_ofile)
                tree_thread.show()
                tree_thread.start()
                tasks.append(tree_thread)

            # Download phylome data
            if not yet_downloaded(workdir, data_ofile):
                data_thread = thread(line, data_url, workdir, data_ofile)
                data_thread.show()
                data_thread.start()
                tasks.append(data_thread)

            # Download protein data
            if not yet_downloaded(workdir, gene_ofile):
                gene_thread = thread(line, gene_url, workdir, gene_ofile)
                gene_thread.show()
                gene_thread.start()
                tasks.append(gene_thread)

            # Download protein data
            if not yet_downloaded(workdir, prot_ofile):
                prot_thread = thread(line, prot_url, workdir, prot_ofile)
                prot_thread.show()
                prot_thread.start()
                tasks.append(prot_thread)

    for task in tasks:
        # Finish remaining tasks
        task.join()
        task.show()


if __name__ == '__main__':
    main()
