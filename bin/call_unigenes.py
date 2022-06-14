#!/usr/bin/env python3

"""

From Dr. Abrams:

1)  Run CD-HIT EST with current parameter set on PASA & DN Trinity contigs that passed 
    TransRate with our scaled Cseg score (0.4 cutoff)

2)  For each cluster call a unigene, based on ORF length.  When there is a tie among 
    longest ORFs in a single cluster, choose longest contig among those that are tied

3)  For these new unigenes for PASA + DN Trinity:

    A) Post full contigs on new Blast site "PASA + DN Trinity Unigenes"
       Wayne and Tom will assess these new unigenes & the clusters
    B) Analyze TOI coverage - both nucleotides and AAs
    C)  Run BUSCO on this unigene set  (lower priority)

4)  If no new issues are found, launch annotation

=================================

This script is written to address #2 using the output from #1.

Creates a list of IDs for transcripts which will be the unigenes.

"""

import argparse
import os
import re
from biocode import utils

def main():
    parser = argparse.ArgumentParser( description='Aplysia unigene caller')

    parser.add_argument('-if', '--input_fasta', type=str, required=True, help='Path to one or more input files to be read, comma-separated with no spaces' )
    parser.add_argument('-c', '--cluster_file', type=str, required=True, help='The .clstr file output from CD-HIT' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    seqs = dict()

    for fasta_file in args.input_fasta.split(','):
        print("Reading FASTA records from: {0}".format(fasta_file))
        seqs.update(utils.fasta_dict_from_file(fasta_file))

    uni_ids_fh = open(args.output_file, 'wt')

    uni_ids = list()

    # each member is keyed by transcript/orf id with value as its length
    cluster = dict()

    for line in open(args.cluster_file):
        line = line.rstrip()

        if line.startswith('>'):
            if len(cluster):
                uni_ids.append(process_cluster(cluster, seqs))

            # Init new cluster
            cluster = dict()
        else:
            m = re.match("\d+\s+(\d+)nt, >(.+).p1\.\.\.", line)
            if m:
                seq_len, seq_id = m.groups()
                seq_len = int(seq_len)
                cluster[seq_id] = seq_len
            else:
                raise Exception("This line didn't match the expected format: {0}".format(line))
        

    # Don't forget the last one
    uni_ids.append(process_cluster(cluster, seqs))

    for id in uni_ids:
        uni_ids_fh.write("{0}\n".format(id))

    uni_ids_fh.close()


def process_cluster(cluster, seqs):
    # this could be written much better, but I'm writing fast.
    longest_length = 0
    longest_ids = list()

    # get longest length
    for seq_id in cluster:
        if cluster[seq_id] > longest_length:
            longest_length = cluster[seq_id]

    for seq_id in cluster:
        if cluster[seq_id] == longest_length:
            longest_ids.append(seq_id)

    if len(longest_ids) == 1:
        return longest_ids[0]
    else:
        longest_transcript_len = 0
        longest_transcript_id = None
        for seq_id in longest_ids:
            if len(seqs[seq_id]['s']) > longest_transcript_len:
                longest_transcript_id = seq_id

        return longest_transcript_id
                
    
if __name__ == '__main__':
    main()







