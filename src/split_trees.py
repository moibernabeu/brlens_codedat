#!/usr/bin/env python3

import sys
from glob import glob
from utils import create_folder
import gzip

create_folder('splitted')

files = glob('outputs/*_best_trees.txt.gz')

print(files)

for file in files:
    print(file)
    ph_id = file.rsplit('/', 1)[1].split('_', 1)[0]
    print(ph_id)
    fs = sys.getsizeof(gzip.open(file, 'r').read())
    opts = 500000
    if fs > opts:
        efno = fs / opts
        fmaxsize = efno / int(efno) * opts + 100

    olines = ''
    ofileno = 0
    for line in gzip.open(file, 'r'):
        if sys.getsizeof(olines) < fmaxsize:
            olines += line.decode('utf-8')
        else:
            print(ofileno, 'written')
            ofile = open('splitted/%s_%s.txt' % (ph_id, ofileno), 'w')
            ofile.write(olines)
            ofile.close()
            olines = ''
            ofileno += 1

    ofile = open('splitted/%s_%s.txt' % (ph_id, ofileno), 'w')
    ofile.write(olines)
    ofile.close()
