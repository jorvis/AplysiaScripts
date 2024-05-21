#!/usr/bin/env python3

import argparse
import os
import re
from biocode import utils, things

def log(level, msg):
    if LOGGING:
        print("{0}: {1}".format(level, msg))


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')
    parser.add_argument('-f', '--source_fasta', type=str, required=True, help='Probably transcript sequences' )
    parser.add_argument('-t', '--transdecoder_fasta', type=str, required=True, help='Output of transdecoder using --source_fasta as input' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-min', '--min_nt_orf_length', type=int, required=False, default=0, help='ORFs less than this (in nt) will not be exported' )
    parser.add_argument('-max', '--max_nt_orf_length', type=int, required=False, default=9999999999999, help='ORFs greater than this (in nt) will not be exported' )
    args = parser.parse_args()

    #MIN_NT_ORF_LEN = 400
    LOGGING = False

    source = utils.fasta_dict_from_file(args.source_fasta)
    transd = utils.fasta_dict_from_file(args.transdecoder_fasta)

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

    ofh = open(args.output_file, 'wt')

    export_count = 0
    skipped_count = 0

    for seq_id in source:
        log('DEBUG', "Processing source sequence id: {0}".format(seq_id))
        log('DEBUG', "\tLongest predicted ORF: {0}".format(longest_orf[seq_id]))

        if longest_orf[seq_id] >= (args.min_nt_orf_length / 3) and longest_orf[seq_id] <= (args.max_nt_orf_length / 3):
            export_count += 1
            log('DEBUG', "\t\tIaM exporting {0} because longest ORF was {1}".format(seq_id, longest_orf[seq_id]))
            ofh.write(">{0}\n".format(seq_id))
            ofh.write(utils.wrapped_fasta(source[seq_id]['s']))
            ofh.write("\n")
        else:
            skipped_count += 1
            log('DEBUG', "\t\tNot exporting {0} because longest ORF was {1}".format(seq_id, longest_orf[seq_id]))

        log('DEBUG', "\t\tExported/skipped count: {0}/{1}".format(export_count, skipped_count))

if __name__ == '__main__':
    main()
