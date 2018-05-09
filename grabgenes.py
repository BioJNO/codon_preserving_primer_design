#!/usr/bin/env python
# grabgenes.py

from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import os

parser = argparse.ArgumentParser(description = "given a list of predicted genes pull out correspondimng fasta records")
parser.add_argument('names', help = 'path to text file with desired fasta record ids as lines')
parser.add_argument('fasta', help = 'path to fasta file containing sequences to be extracted from')
parser.add_argument('-o', default = "grabgenesout", help = 'path to desired output directory')
parser.add_argument('--prefix', default = "extracted_records", help = "output file prefix")
parser.add_argument('-s', action="count", default=0, help = "include to write each gene as a single fasta file")
parser.add_argument('-r', default='', help = "string to be used to rename gene records pulled")
parser.add_argument('-c', action='count', default=0, help = "include to create a config file for find_differential_primers.py, requires -s, default = off")
args = parser.parse_args()

os.mkdir(args.o)

if args.s >=1:
	record_path = os.path.join(args.o, "single_fasta")
	os.mkdir(record_path)
	
if args.c >=1:
	configout = open(os.path.join(args.o, "diff_primer_config.txt"), "w")
	
fasta_handle = open(args.fasta)
f_out = open(os.path.join(args.o, args.prefix + ".fasta"), "w")

gene_list_fh = open(args.names, "r")
gene_list_data = gene_list_fh.readlines()
gene_list = [line.rstrip() for line in gene_list_data
              if line.strip() != ""]
			  
n = 0
lenr = len(args.r)

for title, seq in SimpleFastaParser(fasta_handle):
	if title in gene_list and lenr >= 1:
		n = n + 1
		newname = args.r + str(n)
		f_out.write(">%s\n%s\n" % (newname, seq))
		if args.s >= 1:
			outpath = os.path.join(record_path, newname + ".fasta")
			genout = open(outpath, "w")
			genout.write(">%s\n%s\n" % (newname, seq))
			if args.c >=1:
				configpath = "single_fasta/" + newname + ".fasta"
				configout.write("%s\t%s\t%s\t-\t-\t-\n" % (genname, genname, configpath))
	elif title in gene_list and lenr == 0:
		f_out.write(">%s\n%s\n" % (title, seq))
		if args.s >= 1:
			outpath = os.path.join(record_path, title + ".fasta")
			genout = open(outpath, "w")
			genout.write(">%s\n%s\n" % (title, seq))
			if args.c >=1:
				configpath = "single_fasta/" + title + ".fasta"
				configout.write("%s\t%s\t%s\t-\t-\t-\n" % (title, title, configpath))

f_out.close()

