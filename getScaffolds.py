#!/usr/bin/env python

########################################################################################################################
#
#	getScaffolds.py -i sequnce_id_file -s fasta_file -o output_file
#
#	Written by: Matt Gitzendanner
#				University of Florida
#				Department of Biology
#				magitz@ufl.edu
#
#
# Script to pull sequences from a fasta file based on a list and write to a new file.
# Adapted from: https://gist.github.com/brentp/477969
########################################################################################################################

from Bio import SeqIO
import sys
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file of sequence IDs to get")
parser.add_argument("-s", help="input fasta sequence file")
parser.add_argument("-o", help="Output file of sequences")
args = parser.parse_args()

in_file= args.i
seq_file=args.s
out_file=args.o

wanted = [line.strip() for line in open(in_file)]
seqiter = SeqIO.parse(open(seq_file), 'fasta')
SeqIO.write((seq for seq in seqiter if seq.id in wanted), out_file, "fasta") 

