#!/usr/bin/env python3

"""
- Trinoviewer - put on all of 1A-2C
- BLAST
  - Add, all of above in one file
  - Add, short ORFs in separate file
- Trinotate spreadsheet filtered by the 2A-2C groups

TODO:
- Export full PFAM list to annotation site
- Make sure counts are accounting for all PFAMs
- TRINITY_DN377_c0_g1_i11 has 2 PFAMs displayed currently.  How?

"""

import re
from collections import Counter
from biocode import utils
import numpy as np
import json
import sys

source_fasta = '/usr/local/projects/aplysia/trinotate/final_unigenes/20210222.merged.unigenes.fasta'
trinotate_report = '/usr/local/projects/aplysia/trinotate/final_unigenes/trinotate_annotation_report.xls'
BLASTPcol = 6
BLASTXcol = 2
PFAMcol = 7
GOcol = 12

CONSERVED_NAME = 'Conserved hypothetical protein'

transdecoder_pep_fasta = '/usr/local/projects/aplysia/trinotate/final_unigenes/transdecoder.pep'
small_unigenes_file = '/usr/local/projects/aplysia/trinotate/final_unigenes/20210222.unigenes.small.fasta'
ORF_AA_LENGTH_CUTOFF = 133

## Manually created by those in Tom's lab and colleagues
manual_annot_files = [
    '/usr/local/projects/aplysia/trinotate/final_unigenes/manual_annotation/aplysia_annotation_vs_danio.tab',
    '/usr/local/projects/aplysia/trinotate/final_unigenes/manual_annotation/aplysia_annotation_vs_drosophila_celegans.tab',
    '/usr/local/projects/aplysia/trinotate/final_unigenes/manual_annotation/aplysia_annotation_vs_pomacea.tab'
]

def main():
    seqs = utils.fasta_dict_from_file(source_fasta)
    print("INFO: seq count from source fasta: {0}".format(len(seqs)))
    
    small_count = apply_transdecoder_size_filter(seqs)
    print("INFO: seq count after applying transdecoder size filter: {0}".format(len(seqs) - small_count))

    # Wayne wants lists of the most common BLASTP and BLASTX accessions
    blastp_c = Counter()
    blastx_c = Counter()
    pfam_c   = Counter()
    
    trinotate = parse_trinotate(trinotate_report, blastp_c, blastx_c, pfam_c)
    print("Trinotate gene count: {0}".format(len(trinotate) - 1))

    res1a = open("./annotation_1A.fasta", "wt")
    tab1a = open("./annotation_1A.tab", "wt")
    res1b = open("./annotation_1B.fasta", "wt")
    tab1b = open("./annotation_1B.tab", "wt")
    res1c = open("./annotation_1C.fasta", "wt")
    tab1c = open("./annotation_1C.tab", "wt")
    res2a = open("./annotation_2A.fasta", "wt")
    res2b = open("./annotation_2B.fasta", "wt")
    res2c = open("./annotation_2C.fasta", "wt")
    annot = open("./annotation.1A-2C.trinotate.xls", "wt")
    small = open(small_unigenes_file, "wt")

    mannot = load_manual_annotations(manual_annot_files)
    print("Loaded {0} manual annotation entries".format(len(mannot)))

    annot_export = list()
    
    for seq_id in trinotate:
        match_found = False
        
        # it may have been removed in a previous filter
        if seq_id not in seqs:
            continue
        
        rec = trinotate[seq_id]

        #print("DEBUG: {0} has {1} blastx matches and {2} pfam matches".format(seq_id, len(rec['blastx_acc']), len(rec['pfam_acc'])))
        transcript = {
            'id': seq_id, 'product': 'Hypothetical protein',
            'pfams': list(set(filter(None, trinotate[seq_id]['pfam_acc']))), 'go': list(set(rec['go_terms']))
        }

        # is there a BLASTx hit?
        if any(rec['blastx_acc']):
            # are there any BLASTp hits?
            if any(rec['blastp_acc']):
                # use first blastp row which matches blastx hit, Category 1A match
                for (blastx_acc, blastp_acc, blastp_prod, source_line) in zip(rec['blastx_acc'], rec['blastp_acc'], rec['blastp_prod'], rec['line']):
                    if blastx_acc == blastp_acc:
                        if not match_found:
                            if seq_id in mannot:
                                blastp_prod = mannot[seq_id]

                            if seqs[seq_id]['small']:
                                print_record(seq_id, blastp_prod, seqs[seq_id]['s'], small)
                            else:
                                print_record(seq_id, blastp_prod, seqs[seq_id]['s'], res1a)
                                
                            transcript['product'] = blastp_prod
                            tab1a.write(source_line)
                            annot.write(source_line)
                            match_found = True

                if not match_found:
                    # Category 1C match
                    # get first BLASTx accession
                    for (blastx_acc, blastx_prod, source_line) in zip(rec['blastx_acc'], rec['blastx_prod'], rec['line']):
                        if blastx_acc:
                            if not match_found:
                                if seq_id in mannot:
                                    blastx_acc = mannot[seq_id]

                                if seqs[seq_id]['small']:
                                    print_record(seq_id, blastx_prod, seqs[seq_id]['s'], small)
                                else:
                                    print_record(seq_id, blastx_prod, seqs[seq_id]['s'], res1c)
                                    
                                transcript['product'] = blastx_prod
                                tab1c.write(source_line)
                                annot.write(source_line)
                                match_found = True

                ## sanity check
                if not match_found:
                    raise Exception("Unhandled: BLASTp hits found but no passing ones or back-up BLASTx")
            else:
                # No BLASTp hits, use BlastX call correlating to a Pfam call if available.  Category 1B match
                for (blastx_prod, pfam_prod) in zip(rec['blastx_prod'], rec['pfam_prod']):
                    if pfam_prod:
                        if not match_found:
                            if seq_id in mannot:
                                blastx_prod = mannot[seq_id]

                            if seqs[seq_id]['small']:
                                print_record(seq_id, blastx_prod, seqs[seq_id]['s'], small)
                            else:
                                print_record(seq_id, blastx_prod, seqs[seq_id]['s'], res1b)
                                
                            transcript['product'] = blastx_prod
                            tab1b.write(rec['line'][0])
                            annot.write(rec['line'][0])
                            match_found = True
                else:
                    if not match_found:
                        if seq_id in mannot:
                            if seqs[seq_id]['small']:
                                print_record(seq_id, annot[seq_id], seqs[seq_id]['s'], small)
                            else:
                                print_record(seq_id, annot[seq_id], seqs[seq_id]['s'], res1b)
                                
                            transcript['product'] = annot[seq_id]
                        else:
                            if seqs[seq_id]['small']:
                                print_record(seq_id, "Hypothetical protein", seqs[seq_id]['s'], small)
                            else:
                                print_record(seq_id, "Hypothetical protein", seqs[seq_id]['s'], res1b)
                            
                        annot.write(rec['line'][0])
                        match_found = True

        # No BLASTx hit
        else:
            # are there any BLASTp hits?
            if any(rec['blastp_acc']):
                for (blastp_acc, line) in zip(rec['blastp_acc'], rec['line']):
                    if blastp_acc:
                        if not match_found:
                            if seq_id in mannot:
                                if seqs[seq_id]['small']:
                                    print_record(seq_id, mannot[seq_id], seqs[seq_id]['s'], small)
                                else:
                                    print_record(seq_id, mannot[seq_id], seqs[seq_id]['s'], res2a)
                                    
                                transcript['product'] = mannot[seq_id]
                            else:
                                if seqs[seq_id]['small']:
                                    print_record(seq_id, 'Conserved hypothetical protein', seqs[seq_id]['s'], small)
                                else:
                                    print_record(seq_id, 'Conserved hypothetical protein', seqs[seq_id]['s'], res2a)

                                print(f"Transcript {seq_id} was marked for manual annotation but one wasn't found", file=sys.stderr)
                                #raise Exception(f"Transcript {seq_id} was marked for manual annotation but one wasn't found")

                            annot.write(line)
                            match_found = True
            else:
                if any(rec['pfam_acc']):
                    for (pfam_acc, line) in zip(rec['pfam_acc'], rec['line']):
                        if pfam_acc:
                            if not match_found:
                                if seq_id in mannot:
                                    if seqs[seq_id]['small']:
                                        print_record(seq_id, mannot[seq_id], seqs[seq_id]['s'], small)
                                    else:
                                        print_record(seq_id, mannot[seq_id], seqs[seq_id]['s'], res2b)
                                        
                                    transcript['product'] = mannot[seq_id]
                                else:
                                    if seqs[seq_id]['small']:
                                        print_record(seq_id, 'Conserved hypothetical protein', seqs[seq_id]['s'], small)
                                    else:
                                        print_record(seq_id, 'Conserved hypothetical protein', seqs[seq_id]['s'], res2b)
                                        
                                    print(f"Transcript {seq_id} was marked for manual annotation but one wasn't found", file=sys.stderr)
                                    #raise Exception(f"Transcript {seq_id} was marked for manual annotation but one wasn't found")
                                    
                                annot.write(line)
                                match_found = True
                else:
                    if not match_found:
                        if seqs[seq_id]['small']:
                            print_record(seq_id, 'Hypothetical protein', seqs[seq_id]['s'], small)
                        else:
                            print_record(seq_id, 'Hypothetical protein', seqs[seq_id]['s'], res2c)
                            
                        annot.write(rec['line'][0])
                        match_found = True

        if not match_found:
            print("Transcript {0} fell through the ladder".format(seq_id))

        annot_export.append(transcript)
                    
    res1a.close()
    tab1a.close()
    res1b.close()
    tab1b.close()
    res1c.close()
    tab1c.close()
    res2a.close()
    res2b.close()
    res2c.close()
    annot.close()
    small.close()

    blastx_c_fh = open("./annotation.blastx.top.accessions.txt", "wt")
    blastx_c_fh.write("Most common BLASTx accessions:\n")
    for (acc, n) in blastx_c.most_common(100):
        blastx_c_fh.write("\t{0}\t{1}\n".format(acc, n))
    blastx_c_fh.close()

    pfam_counts = Counter()
    for entry in annot_export:
        for id in entry['pfams']:
            pfam_counts[id] += 1

    pfam_c_fh = open("./annotation.pfam.top.accessions.txt", "wt")
    pfam_c_fh.write("Most common PFAM accessions:\n")
    for (acc, n) in pfam_counts.most_common(1000000):
        pfam_c_fh.write("\t{0}\t{1}\n".format(acc, n))
    pfam_c_fh.close()

    json_annot_fh = open("./current.annotation.json", "wt")
    json_annot_fh.write(json.dumps(annot_export, indent=3))
    json_annot_fh.close()



def apply_transdecoder_size_filter(seqs):
    longest_orf_per_seq = dict()
    small_count = 0

    for line in open(transdecoder_pep_fasta):
        m = re.match(">(\S+)\.p.* len:(\d+).*", line)
        if m:
            seq_id = m.group(1)
            orf_len = int(m.group(2))

            if seq_id not in longest_orf_per_seq:
                longest_orf_per_seq[seq_id] = orf_len
            elif orf_len > longest_orf_per_seq[seq_id]:
                longest_orf_per_seq[seq_id] = orf_len

    ids_to_delete = list()
                
    for id in seqs:
        if id not in longest_orf_per_seq or longest_orf_per_seq[id] < ORF_AA_LENGTH_CUTOFF:
            small_count += 1
            seqs[id]['small'] = 1
        else:
            seqs[id]['small'] = 0

    return small_count

def load_manual_annotations(files):
    mannot = dict()

    for file in files:
        for line in open(file):
            line = line.rstrip()
            cols = line.split("\t")
            acc = cols[0]

            if acc in mannot:
                continue

            prod = process_product_name(cols[1])
            mannot[acc] = prod

    return mannot

def get_go_terms(go_str):
    go_terms = list()

    entries = go_str.split('`')
    for entry in entries:
        m = re.match("^(.*?)\^", entry)
        if m:
            go_terms.append(m.group(1))

    return go_terms

def print_record(id, desc, seq, fh):
    fh.write(">{0} {1}\n".format(id, desc))
    for i in range(0, len(seq), 60):
        fh.write(seq[i : i + 60] + "\n")

def process_product_name(prod_in):
    prod_out = prod_in

    if re.search('[Uu]ncharacterized', prod_out):
        prod_out = CONSERVED_NAME

    if re.search('-like', prod_out):
        prod_out = prod_out.replace('-like', '')

    if re.search('_like', prod_out):
        prod_out = prod_out.replace('_like', '')

    if re.search(' LOW_QUALITY_PROTEIN:_', prod_out):
        prod_out = prod_out.replace(' LOW_QUALITY_PROTEIN:_', '')

    # Si:zfos-1056e6.1

    prod_out = prod_out.replace(' dom domain', ' domain')
    prod_out = prod_out.replace('ethyltranfer ', 'ethyltransferase ')
    prod_out = prod_out.replace('Acetyltransf ', 'Acetyltransferase ')

    m = re.match(".._.._\d+\.\d(.+)", prod_out)
    if m:
        prod_out = m.group(1)
        prod_out = prod_out.replace('_', ' ')

    m = re.match("XP_\d+\.\d+_(.+)", prod_out)
    if m:
        prod_out = m.group(1)
        prod_out = prod_out.replace('_', ' ')

    prod_out = prod_out.replace('  ', ' ')
    prod_out = prod_out.replace('  ', ' ')
    prod_out = prod_out.strip().capitalize()

    #print(prod_out)

    return prod_out
    
def parse_trinotate(infile, blastp_c, blastx_c, pfam_c):
    ds = dict()

    for line in open(infile):
        cols = line.rstrip().split("\t")
        seq_id = cols[0]

        #print("Processing seq_id: {0}".format(seq_id))

        if seq_id not in ds:
            ds[seq_id] = { 'blastp_acc': [],  'blastx_acc': [], 'pfam_acc': [],
                           'blastp_prod': [], 'blastx_prod': [], 'pfam_prod': [],
                           'line': []
            }

        # save the entire line for export reasons
        ds[seq_id]['line'].append(line)

        if '^' in cols[BLASTPcol]:
            ds[seq_id]['blastp_acc'].append(cols[BLASTPcol].split('^')[0])
        else:
            ds[seq_id]['blastp_acc'].append(None)

        if '^' in cols[BLASTXcol]:
            ds[seq_id]['blastx_acc'].append(cols[BLASTXcol].split('^')[0])
        else:
            ds[seq_id]['blastx_acc'].append(None)

        if '^' in cols[PFAMcol]:
            #ds[seq_id]['pfam_acc'].append(cols[PFAMcol].split('^')[0].split('.')[0])
            #ds[seq_id]['pfam_prod'].append(cols[PFAMcol].split('^')[2])

            entries = cols[PFAMcol].split('`')
            for entry in entries:
                pfam_acc = entry.split('^')[0].split('.')[0]
                pfam_prod = entry.split('^')[2]
                
                if pfam_acc not in ds[seq_id]['pfam_acc']:
                    ds[seq_id]['pfam_acc'].append(pfam_acc)
                    ds[seq_id]['pfam_prod'].append(pfam_prod)
            
        else:
            ds[seq_id]['pfam_acc'].append(None)
            ds[seq_id]['pfam_prod'].append(None)

        m = re.search("RecName: Full=(.+?)\;", cols[BLASTPcol])
        if m:
            product = process_product(m.group(1))
            ds[seq_id]['blastp_prod'].append(product)
        else:
            ds[seq_id]['blastp_prod'].append(None)

        m = re.search("RecName: Full=(.+?)\;", cols[BLASTXcol])
        if m:
            product = process_product(m.group(1))
            ds[seq_id]['blastx_prod'].append(product)
        else:
            ds[seq_id]['blastx_prod'].append(None)

        ## GO terms
        go_terms = get_go_terms(cols[GOcol])
        go_terms.extend(get_go_terms(cols[GOcol + 1]))
        go_terms.extend(get_go_terms(cols[GOcol + 2]))
        ds[seq_id]['go_terms'] = np.unique(go_terms)

        ## now do the counts
        for pfam_entry in cols[PFAMcol].split('`'):
            m = re.match("^(.+?)\.", pfam_entry)
            if m:
                #if m.group(1) == "PF00759":
                #    print(cols[PFAMcol])

                pfam_c[m.group(1)] += 1
                #print("\t{0}".format(m.group(1)))
        
        for blastx_entry in cols[BLASTXcol].split('`'):
            m = re.match("^(.+?)\.", blastx_entry)
            if m:
                blastx_c[m.group(1)] += 1

        #print("Seq id:{0} has {1} blastp entries and {2} HMM entries".format(seq_id,
        #                                                                     len(ds[seq_id]['blastp_acc']),
        #                                                                     len(ds[seq_id]['pfam_acc'])))
                
    return ds

def process_product(prod):
    # if the product are like this, strip off the brackets and everything within:
    #  Dynein heavy chain 5, axonemal {ECO:0000305}
    m = re.match("(.+) \{ECO.*", prod)
    if m:
        return m.group(1)
    else:
        return prod

if __name__ == '__main__':
    main()


















    
