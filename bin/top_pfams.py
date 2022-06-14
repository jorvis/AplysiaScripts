#!/usr/bin/env python3

"""
Count top PFAMs from all evidence
"""

import re
from collections import Counter
from biocode import utils

source_fasta = '/usr/local/projects/aplysia/trinotate/final_unigenes/20210222.merged.unigenes.fasta'
trinotate_report = '/usr/local/projects/aplysia/trinotate/final_unigenes/trinotate_annotation_report.xls'
BLASTPcol = 6
BLASTXcol = 2
PFAMcol = 7

transdecoder_pep_fasta = '/usr/local/projects/aplysia/trinotate/final_unigenes/transdecoder.pep'
transdecoder_dest_file = '/usr/local/projects/aplysia/trinotate/final_unigenes/20210222.unigenes.small.fasta'
ORF_AA_LENGTH_CUTOFF = 133

def main():
    seqs = utils.fasta_dict_from_file(source_fasta)
    apply_transdecoder_size_filter(seqs)

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
    
    for seq_id in trinotate:
        match_found = False
        
        # it may have been removed in a previous filter
        if seq_id not in seqs:
            continue
        
        rec = trinotate[seq_id]

        #print("DEBUG: {0} has {1} blastx matches and {2} pfam matches".format(seq_id, len(rec['blastx_acc']), len(rec['pfam_acc'])))

        # is there a BLASTx hit?
        if any(rec['blastx_acc']):
            # are there any BLASTp hits?
            if any(rec['blastp_acc']):
                # use first blastp row which matches blastx hit, Category 1A match
                for (blastx_acc, blastp_acc, blastp_prod, source_line) in zip(rec['blastx_acc'], rec['blastp_acc'], rec['blastp_prod'], rec['line']):
                    if blastx_acc == blastp_acc:
                        if not match_found:
                            print_record(seq_id, blastp_prod, seqs[seq_id]['s'], res1a)
                            tab1a.write(source_line)
                            annot.write(source_line)
                            match_found = True

                if not match_found:
                    # Category 1C match
                    # get first BLASTx accession
                    for (blastx_acc, source_line) in zip(rec['blastx_acc'], rec['line']):
                        if blastx_acc:
                            if not match_found:
                                print_record(seq_id, blastx_acc, seqs[seq_id]['s'], res1c)
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
                            print_record(seq_id, blastx_prod, seqs[seq_id]['s'], res1b)
                            tab1b.write(rec['line'][0])
                            annot.write(rec['line'][0])
                            match_found = True
                else:
                    if not match_found:
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
                            print_record(seq_id, 'Pending manual annotation', seqs[seq_id]['s'], res2a)
                            annot.write(line)
                            match_found = True
            else:
                if any(rec['pfam_acc']):
                    for (pfam_acc, line) in zip(rec['pfam_acc'], rec['line']):
                        if pfam_acc:
                            if not match_found:
                                print_record(seq_id, 'Pending manual annotation', seqs[seq_id]['s'], res2b)
                                annot.write(line)
                                match_found = True
                else:
                    if not match_found:
                        print_record(seq_id, 'Hypothetical protein', seqs[seq_id]['s'], res2c)
                        annot.write(rec['line'][0])
                        match_found = True
                    
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

    blastx_c_fh = open("./annotation.blastx.top.accessions.txt", "wt")
    blastx_c_fh.write("Most common BLASTx accessions:\n")
    for (acc, n) in blastx_c.most_common(100):
        blastx_c_fh.write("\t{0}\t{1}\n".format(acc, n))
    blastx_c_fh.close()

    pfam_c_fh = open("./annotation.pfam.top.accessions.txt", "wt")
    pfam_c_fh.write("Most common PFAM accessions:\n")
    for (acc, n) in pfam_c.most_common(1000000):
        pfam_c_fh.write("\t{0}\t{1}\n".format(acc, n))
    pfam_c_fh.close()
        


def apply_transdecoder_size_filter(seqs):
    longest_orf_per_seq = dict()

    for line in open(transdecoder_pep_fasta):
        m = re.match(">(\S+)\.p.* len:(\d+).*", line)
        if m:
            seq_id = m.group(1)
            orf_len = int(m.group(2))

            if seq_id not in longest_orf_per_seq:
                longest_orf_per_seq[seq_id] = orf_len
            elif orf_len > longest_orf_per_seq[seq_id]:
                longest_orf_per_seq[seq_id] = orf_len

    fout = open(transdecoder_dest_file, 'wt')

    ids_to_delete = list()
                
    for id in seqs:
        if id not in longest_orf_per_seq or longest_orf_per_seq[id] < ORF_AA_LENGTH_CUTOFF:
            print_record(id, seqs[id]['h'], seqs[id]['s'], fout)
            ids_to_delete.append(id)

    for id in ids_to_delete:
        del seqs[id]

    fout.close()

def print_record(id, desc, seq, fh):
    fh.write(">{0} {1}\n".format(id, desc))
    for i in range(0, len(seq), 60):
        fh.write(seq[i : i + 60] + "\n")
    
def parse_trinotate(infile, blastp_c, blastx_c, pfam_c):
    ds = dict()

    for line in open(infile):
        cols = line.rstrip().split("\t")
        seq_id = cols[0]

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
            ds[seq_id]['pfam_acc'].append(cols[PFAMcol].split('^')[0])
            ds[seq_id]['pfam_prod'].append(cols[PFAMcol].split('^')[2])
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

        ## now do the counts
        for pfam_entry in cols[PFAMcol].split('`'):
            m = re.match("^(.+?)\.", pfam_entry)
            if m:
                if m.group(1) == "PF00759":
                    print(cols[PFAMcol])

                pfam_c[m.group(1)] += 1
        
        for blastx_entry in cols[BLASTXcol].split('`'):
            m = re.match("^(.+?)\.", blastx_entry)
            if m:
                blastx_c[m.group(1)] += 1
                
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


















    
