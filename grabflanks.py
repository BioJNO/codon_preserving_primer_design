#!/usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser(description = "pull out desired genes from genomic sequence with flanking regions")
parser.add_argument('genes', help = 'path to fasta file containing desired gene sequences')
parser.add_argument('contigs', help = 'path to fasta file containing genomic contigs')
parser.add_argument('-f', type=int, default = 200, help = 'desired length of genomic flanking region to return either side of desired gene sequences, default=200')
parser.add_argument('-o', default="grabflanksout", help = 'path to desired output directory, default=grabflanksout')
parser.add_argument('-v', action="count", default=0, help = 'set verbose mode, default=off')
parser.add_argument('--prefix', default = "genes", help = 'add output file prefix, default=genes')
parser.add_argument('-s', action="count", default=0, help = "include to write each recovered gene with falnking region as a single fasta file, default=off")
parser.add_argument('-c', action="count", default=0, help = "include to write config files to pass to find_differential_primers.py for each gene flanking region, requires -s, default=off")
args = parser.parse_args()

genehandle = open(args.genes)
genomefilename = args.contigs
os.mkdir(args.o)

if args.c >=1:
	shellout = open("diff_primers.sh", "w")

if args.s >=1:
	record_path = os.path.join(args.o, "single_fasta")
	os.mkdir(record_path)


outfile = open(os.path.join(args.o, args.prefix + "_w_flanking_regions.fasta"), "w")
tabout = open(os.path.join(args.o, "grabflanks_tabular_out.csv"), "w")
tabout.write("gene name,contigname,gene start,gene end,strand,gene length\n")
logfile = open(os.path.join(args.o, "grabflanks_log.txt"), "w")

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna 
from Bio.SeqIO.FastaIO import SimpleFastaParser

def grabflanks (genefasta, genomefasta, flank):
	for title, seq in SimpleFastaParser(genefasta):
		genname = title
		genseq = Seq(seq, generic_dna)
		genstring = str(genseq)
		genseqrev = str(genseq.reverse_complement())
		genlen = len(genstring)
		if args.v >= 1:
			print("length of " + genname + " is " + str(genlen))
		logfile.write("length of " + genname + " is " + str(genlen) + "\n")
		genome = open(genomefasta)
		
		for title,seq in SimpleFastaParser(genome):
			if args.v >= 1:
				print("searching for " + genname + " in contig " + title)
			findgen = seq.find(genstring.upper())
			findgenrev = seq.find(genseqrev.upper())
			contigname = title.replace("|", "_")
			contigseq = seq
			
			if findgen > -1:
				end = findgen + genlen
				if args.v >= 1:
					print("found " + genname + " in contig " + contigname + " from " + str(findgen) + " to " + str(end) + " on the forward strand")
				logfile.write("found " + genname + " in contig " + contigname + " from " + str(findgen) + " to " + str(end)  + " on the forward strand\n")
				genflank = contigseq[findgen - flank:findgen+genlen + flank]
				flankname = genname + "_" + contigname + "_flanking"
				outfile.write(">%s\n%s\n" % (flankname, genflank))
				tabout.write("%s,%s,%s,%s,forward,%s\n" % (genname, contigname, findgen, end, genlen))
				
				if args.s >= 1:
					singlepath = os.path.join(record_path, flankname + ".fasta")
					singleout = open(singlepath, "w")
					singleout.write(">%s\n%s\n" % (flankname, genflank))
					
					if args.c >=1:
						configname = os.path.join(args.o, genname + "_config.txt")
						configout = open(configname, "w")
						configout.write("%s\t%s\t%s\t-\t-\t-\n" % (genname, genname, singlepath))
						shellout.write("\nfind_differential_primers.py -i " + configname + " -o " + genname + "_amplicon_primers -v --prodigal ${PRODIGAL} --eprimer3 ${EPRIMER3} --blast_exe ${BLAST} --blastdb nt --cpus 8 --noprodigal --nocds --psizemin " + str(genlen) + " --psizemax " + str(genlen + args.f) + " --target " + str(args.f) + "," + str(genlen + args.f) + "\nwait\n")
						
			elif findgenrev > -1:
				end = findgenrev + genlen
				if args.v >=1:
					print("found " + genname + " in contig " + contigname + " from " + str(findgenrev) + " to " + str(end) + " on the reverse strand")
				logfile.write("found " + genname + " in contig " + contigname + " from " + str(findgenrev) + " to " + str(end) + " on the reverse strand\n")
				genflank = contigseq[findgenrev - flank:findgenrev+genlen + flank]
				flankname = genname + "_" + contigname + "_flanking"
				outfile.write(">%s\n%s\n" % (flankname, genflank))
				tabout.write("%s,%s,%s,%s,reverse,%s\n" % (genname, contigname, findgenrev, end, genlen))
				
				if args.s >= 1:
					singlepath = os.path.join(record_path, flankname + ".fasta")
					singleout = open(singlepath, "w")
					singleout.write(">%s\n%s\n" % (flankname, genflank))
					
					if args.c >= 1:
						configname = os.path.join(args.o, genname + "_config.txt")
						configout = open(configname, "w")
						configout.write("%s\t%s\t%s\t-\t-\t-\n" % (genname, genname, singlepath))
						shellout.write("\nfind_differential_primers.py -i " + configname + " -o " + genname + "_amplicon_primers -v --prodigal ${PRODIGAL} --eprimer3 ${EPRIMER3} --blast_exe ${BLAST} --blastdb nt --cpus 8 --noprodigal --nocds --psizemin " + str(genlen) + " --psizemax " + str(genlen + args.f) + " --target " + str(args.f) + "," + str(genlen + args.f) + "\nwait\n")

	
	
grabflanks(genehandle, genomefilename, args.f)



outfile.close()
logfile.close()
tabout.close()

