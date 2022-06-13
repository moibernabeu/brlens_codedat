#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
thread.py -- Construct the threads to download the trees and its data

Written by Name <mail@mail.com>
Month 2022
'''

# Import libraries ----
import threading
from download import file_download


# Definitions ----
class thread(threading.Thread):
    '''
    Run download in a thread
    '''

    # Class constructor
    def __init__(self, dbid, url, workdir, ofile):
        threading.Thread.__init__(self)
        self.dbid = dbid
        self.url = url
        self.workdir = workdir
        self.ofile = ofile
        self.message = 'thread for the item %s is created' % self.dbid

    # Show method
    def show(self):
        '''
        Prints thread output
        '''

        print('Download: %s' % self.message)

    # Run method
    def run(self):
        '''
        Run the thread
        '''

        self.message = 'Thread for the item %s is running' % self.dbid

        try:
            # Download the item
            item_download = file_download(self.url, self.workdir, self.ofile)
            item_download.run()
        except Exception as e:
            self.message = 'Thread %s exception: %s' % (self.dbid, e.args)
        else:
            self.message = 'Thread %s ended successfully' % self.dbid
