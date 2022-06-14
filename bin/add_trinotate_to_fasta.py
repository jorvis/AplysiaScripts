#!/usr/bin/env python3

"""

Put some general, high-level documentation here

"""

import argparse
import os
import re

from biocode import utils

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    parser.add_argument('-if', '--input_fasta', type=str, required=True, help='Path to an input fasta file to be read' )
    parser.add_argument('-it', '--input_trinotate', type=str, required=True, help='Path to an input trinotate file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    TRINOTATE_PRODUCT_COLUMN = 5

    trin = dict()
    
    for line in open(args.input_trinotate):
        cols = line.split("\t")
        id = cols[0]
        product = 'hypothetical protein'
        
        m = re.search("Full=(.+?)\;", cols[TRINOTATE_PRODUCT_COLUMN])
        if m:
            product = m.group(1)

        trin[id] = product

    seqs = utils.fasta_dict_from_file(args.input_fasta)

    ofh = open(args.output_file, 'wt')

    for seq_id in seqs:
        if seq_id not in trin:
            raise Exception("Seq id {0} not found in trinotate input.".format(seq_id))

    ofh.close()
    


if __name__ == '__main__':
    main()







