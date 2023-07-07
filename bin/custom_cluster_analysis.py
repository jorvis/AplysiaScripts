#!/usr/bin/env python3

import re, sys
from biocode import utils

SOURCE_FASTA = '/usr/local/projects/aplysia/cd-hit/all_merged/A1-B6.unigenes.cdhitest.fasta'
CLUSTER_FILE = '/usr/local/projects/aplysia/cd-hit/all_merged/A1-B6.unigenes.cdhitest.fasta.clstr'

ORF_FILES = [ '/usr/local/projects/aplysia/cd-hit/all_merged/A2-A10.unigenes.faa',
              '/usr/local/projects/aplysia/cd-hit/all_merged/A1.unigenes.prefixed.faa',
              '/usr/local/projects/aplysia/cd-hit/all_merged/B4-B6.unigenes.faa'
            ]

seqs = utils.fasta_dict_from_file(SOURCE_FASTA)

current_cluster_ids = list()
current_cluster_rep = None
cluster_count = 0
longest_orf = dict()
cluster_categories = {
    'cns_only': 0,
    'both_cns_wins': 0,
    'both_peripheral_wins': 0,
    'peripheral_only': 0,
    'ovo_only': 0,
    'peripheral_with_ovo': 0,
    'peripheral_sans_ovo': 0
}

def process_cluster(ids, rep):
    cluster_type = 'competitive'
    
    if len(ids) == 1:
        cluster_type = 'singular'

    has_CNS = False
    has_nonCNS = False
    longest_orf_len = 0
    longest_orf_id = None
    ovo_count = 0

    for id in ids:
        print("DEBUG: {0} - {1}".format(longest_orf[id], id))
        if id.startswith('CNS'):
            has_CNS = True
        else:
            if id.startswith('ovotestis'):
                ovo_count += 1
                
            has_nonCNS = True

        this_orf_len = longest_orf[id]

        if this_orf_len > longest_orf_len:
            longest_orf_len = this_orf_len
            longest_orf_id = id

    ## CNS only
    if has_CNS and not has_nonCNS:
        cluster_categories['cns_only'] += 1
        print("Winner: CNS only\n")

    ## Both
    elif has_CNS and has_nonCNS:
        # but CNS wins
        if longest_orf_id.startswith('CNS'):
            cluster_categories['both_cns_wins'] += 1
            print("Winner: Both, CNS wins\n")
            
        ## but peripheral wins
        else:
            cluster_categories['both_peripheral_wins'] += 1
            print("Winner: Both, periphery wins\n")

    ## Peripheral only
    else:
        cluster_categories['peripheral_only'] += 1
        print("Winner: Peripheral only\n")

        # We want a sub-filter to discern mixed peripheral vs Ovotestis only
        if ovo_count:
            if ovo_count == len(ids):
                cluster_categories['ovo_only'] += 1
            else:
                cluster_categories['peripheral_with_ovo'] += 1
        else:
            cluster_categories['peripheral_sans_ovo'] += 1
            
    
for orf_file in ORF_FILES:
    print("Processing ORF file: {0}".format(orf_file))
    for line in open(orf_file):
        if line[0] == '>':
            m = re.match(">(.*?).p\d+ .* len:(\d+) .*", line)
            if m:
                (seq_id, seq_len) = m.groups()
                seq_len = int(seq_len)
                
                if seq_id in longest_orf:
                    if seq_len > longest_orf[seq_id]:
                        longest_orf[seq_id] = seq_len
                else:
                    longest_orf[seq_id] = seq_len
            else:
                raise Exception("Unexpected format in peptide file with line:\n{0}".format(line))

for line in open(CLUSTER_FILE):
    line = line.rstrip()

    if line.startswith(">"):
        cluster_count += 1
        # purge cluster and reset
        process_cluster(current_cluster_ids, current_cluster_rep)
        current_cluster_ids = list()
        current_cluster_rep = None

    else:
        m = re.match(".*>(.*?)\.\.\. .*", line)
        if m:
            seq_id = m.group(1)
            current_cluster_ids.append(seq_id)

            # initialize this one if we haven't seen it
            if seq_id not in longest_orf:
                longest_orf[seq_id] = 0

            if line.endswith('*'):
                current_cluster_rep = seq_id
        else:
            raise Exception("Error: Unexpected format in cluster member def line:\n{0}".format(line))

# Don't forget the last one:
process_cluster(current_cluster_ids, current_cluster_rep)

print("total seq: {0}".format(len(longest_orf)))
print("total clusters: {0}".format(cluster_count))

for category in cluster_categories:
    print("{0}\t{1}".format(category, cluster_categories[category]))
