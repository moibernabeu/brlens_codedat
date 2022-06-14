#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
join_normalise.py -- Phylome partitions join and calculate normalised distances



Requirements: pandas

Written by Moisès Bernabeu <moigil.bernabeu.sci@gmail.com>
March 2022
'''

# Import libraries ----
import pandas as pd
from glob import glob

# Script
def main():
    files = glob('outputs/*.csv')
    print(files)

    ids = list()
    for file in files:
        ids.append(file.rsplit('/', 1)[1].split('_', 1)[0])
    ids = set(ids)

    for phyid in ids:
        print('Parsing: ', phyid)
        distfiles = glob('outputs/%s*_dist.csv' % phyid)

        for i, file in enumerate(distfiles):
            if i == 0:
                distdf = pd.read_csv(file)
            else:
                distdf = pd.concat([distdf, pd.read_csv(file)])

        distdf.to_csv('outputs/%s_dist.csv' % phyid, index=False)


if __name__ == '__main__':
    main()
