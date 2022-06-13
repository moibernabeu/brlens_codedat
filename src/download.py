#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
download.py -- Download trees and data from phylome db functions

Written by Moisès Bernabeu <mail@mail.com>
January 2022
'''

# Import libraries ----
from os.path import isfile
from urllib.request import urlretrieve
import socket
from utils import create_folder


# Definitions ----
class download_error(Exception):
    '''
    Download error
    '''
    pass


class file_download(object):
    '''
    Tree download class
    '''

    # Class constructor
    def __init__(self, url, workdir, ofile, timeout=180):
        '''
        Construct the download_tree class with its own objects
        '''

        self.url = url
        self.timeout = timeout
        self.workdir = workdir
        self.ofile = ofile

    # Make containing folder
    def make_workdir(self, workdir):
        create_folder(self.workdir)

    # Download the item
    def get_item(self):
        '''
        Download the tree from PhylomeDB if it is necessary
        '''

        # Create the directory
        self.make_workdir(self.workdir)

        # Check the existence of the file and download if it is not created
        dest = '/'.join([self.workdir, self.ofile])

        if not isfile(dest):
            try:
                socket.setdefaulttimeout(self.timeout)
                urlretrieve(self.url, dest)
            except socket.timeout:
                raise download_error(('Timeout (%d): %s'
                                      % (self.timeout, self.url)))
            except Exception:
                raise download_error('Cannot download %s' % self.url)

    # Run the download
    def run(self):
        '''
        Download
        '''

        self.get_item()
