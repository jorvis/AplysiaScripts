#!/usr/bin/env python3

import re
from biocode import utils, things

def log(level, msg):
    if LOGGING:
        print("{0}: {1}".format(level, msg))

SOURCE_FASTA = '/usr/local/projects/aplysia/ncbi/20210222.merged.unigenes.vecscreened.filtered.notempty.ncbifiltered.fasta'
#TRANSDECODER_FASTA = '/usr/local/projects/aplysia/transdecoder/final_unigenes/20210222.merged.unigenes.fasta.transdecoder_dir/longest_orfs.cds'
TRANSDECODER_FASTA = '/usr/local/projects/aplysia/ncbi/20210222.merged.unigenes.vecscreened.filtered.notempty.ncbifiltered.fasta.transdecoder_dir/longest_orfs.cds'
OUTPUT_FASTA = '/usr/local/projects/aplysia/ncbi/20210222.merged.unigenes.vecscreened.filtered.notempty.ncbifiltered.orf400.fasta'
MIN_NT_ORF_LEN = 400
LOGGING = False

source = utils.fasta_dict_from_file(SOURCE_FASTA)
transd = utils.fasta_dict_from_file(TRANSDECODER_FASTA)

longest_orf = dict()

for seq_id in transd:
    log('DEBUG', "Processing transdecoder id: {0}".format(seq_id))
    # need to strip off the .p1, .p2, etc from the IDs to make them comparable
    m = re.match("(.+)\.p\d+", seq_id)
    if m:
        transcript_id = m.group(1)
        log('DEBUG', "\tchopped name to: {0}".format(transcript_id))
    else:
        raise Exception("Unexpected transcript ID convention")

    # The len() value here is amino acid length
    m = re.match(".* len\:(\d+) .*", transd[seq_id]['h'])
    if m:
        orf_len = int(m.group(1))
        
        if transcript_id in longest_orf:
            if orf_len > longest_orf[transcript_id]:
                longest_orf[transcript_id] = orf_len
        else:
            longest_orf[transcript_id] = orf_len
    else:
        raise Exception("Unexpected length header line pattern")

ofh = open(OUTPUT_FASTA, 'wt')

export_count = 0
skipped_count = 0

for seq_id in source:
    log('DEBUG', "Processing source sequence id: {0}".format(seq_id))
    log('DEBUG', "\tLongest predicted ORF: {0}".format(longest_orf[seq_id]))
    
    if longest_orf[seq_id] >= (MIN_NT_ORF_LEN / 3):
        export_count += 1
        log('DEBUG', "\t\tIaM exporting {0} because longest ORF was {1}".format(seq_id, longest_orf[seq_id]))
        ofh.write(">{0}\n".format(seq_id))
        ofh.write(utils.wrapped_fasta(source[seq_id]['s']))
        ofh.write("\n")
    else:
        skipped_count += 1
        log('DEBUG', "\t\tNot exporting {0} because longest ORF was {1}".format(seq_id, longest_orf[seq_id]))

    log('DEBUG', "\t\tExported/skipped count: {0}/{1}".format(export_count, skipped_count))
