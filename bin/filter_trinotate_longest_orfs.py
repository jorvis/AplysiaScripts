#!/usr/bin/env python3

import sys
import re

infile = sys.argv[1]

genes = dict()

for line in open(infile):
    line = line.rstrip()
    cols = line.split("\t")

    asmbl_id = cols[0]
    coords = cols[4]

    if asmbl_id not in genes:
        genes[asmbl_id] = {'length': 0, 'line': None}

    m = re.match("(\d+)\-(\d+)\[", coords)
    if m:
        orf_length = abs(int(m.group(2)) - int(m.group(1)))
        if orf_length > genes[asmbl_id]['length']:
            genes[asmbl_id]['length'] = orf_length
            genes[asmbl_id]['line'] = line
    else:
        raise Exception("Unhandled coords: ({0})".format(coords))

for asmbl_id in genes:
    print(genes[asmbl_id]['line'])
