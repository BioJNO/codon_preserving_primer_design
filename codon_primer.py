#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna 
from Bio.SeqIO.FastaIO import SimpleFastaParser
import primer3
import argparse

parser = argparse.ArgumentParser(description = "given a fasta file of coding sequences presents all possible PCR primers which preserve the coding sense of the sequence")
parser.add_argument('fasta', help = "path to fasta file containing desired coding sequences")
parser.add_argument('-pl', type=int, default = 22, help = "desired length of primers, default=22")
parser.add_argument('-f', type=int, default=0, help = "number of base pairs flanking coding region, default = 0")
parser.add_argument('-o', type=str, default="CDS_primers_out.txt", help = "specify output filename, default = CDS_primers_out.txt")
parser.add_argument('-tmax', type=int, default=58, help = "maximum melting temp for primers to report, default = 58")
parser.add_argument('-tmin', type=int, default=52, help = "minimum melting temp for primers to report, default = 52")
parser.add_argument('-os', default = "single_rec_output.txt", help = "output table for unique primers, default = single_rec_output.txt")
parser.add_argument('-om', default = "multi_rec_output.txt", help = "output table for non-unique primers, default = multi_rec_output.txt")

args = parser.parse_args()

seqhandle = open(args.fasta)
outfile = open(args.o, "w")
outfile.write("Target CDS\tPrimer type\tprimer sequence\tTm\tprimer length\tamp start pos\tframelen\torientation\n")

poslist = ['ATG', 'GTG', 'TTG']
neglist = ['TTA', 'CTA', 'TCA']

for title, seq in SimpleFastaParser(seqhandle):
	lenseq = len(seq)
	frame = seq[args.f:lenseq -  args.f].upper()
	if frame[:3] in poslist:
		orientation = 1
	elif frame[:3] in neglist:
		orientation = -1
	framelen = len(frame)
	postframe = seq[lenseq - args.f:].upper()
	preframe = seq[: args.f].upper()
	preframe_obj = Seq(preframe, generic_dna)
	preframerev = str(preframe_obj.reverse_complement())
	frame_obj = Seq(frame, generic_dna)
	framerev = str(frame_obj.reverse_complement())
	postframe_obj = Seq(postframe, generic_dna)
	postframerev = str(postframe_obj.reverse_complement())
	
	for index in range(0, args.f - args.pl, 1):
		fcandidate = preframe[index:index+args.pl]
		fcand_len = str(len(fcandidate))
		Tm = int(primer3.calcTm(fcandidate))
		if Tm > args.tmin and Tm < args.tmax:
			Hairpin = primer3.calcHairpin(fcandidate)
			if Hairpin.structure_found is False:
				outfile.write(title + "\tPreframe forward\t" + fcandidate + "\t" + str(Tm) + "\t" + fcand_len + "\t" + str(index) + "\t" + str(framelen) + "\t" + str(orientation) + "\n")
		
	for index in range(0, framelen - args.pl, 3):
		fcandidate = frame[index:index+args.pl]
		fcand_len = str(len(fcandidate))
		Tm = int(primer3.calcTm(fcandidate))
		if Tm > args.tmin and Tm < args.tmax:
			Hairpin = primer3.calcHairpin(fcandidate)
			if Hairpin.structure_found is False:
				outfile.write(title + "\tIn frame forward\t" + fcandidate + "\t " + str(Tm) + "\t" + fcand_len + "\t" + str(index + args.f) + "\t" + str(framelen) + "\t" + str(orientation) + "\n")
			
	for index in range(0, framelen - args.pl, 3):
		rcandidate = framerev[index:index+args.pl]
		rcand_len = str(len(rcandidate))
		Tm = int(primer3.calcTm(rcandidate))
		if Tm > args.tmin and Tm < args.tmax:
			Hairpin = primer3.calcHairpin(rcandidate)
			if Hairpin.structure_found is False:
				outfile.write(title + "\tIn frame reverse\t" + rcandidate + "\t" + str(Tm) + "\t" + rcand_len + "\t" + str(args.f + (framelen - index)) + "\t" + str(framelen) + "\t" + str(orientation) + "\n")
			
	for index in range(0, args.f - args.pl, 1):
		rcandidate = postframerev[index:index+args.pl]
		rcand_len = str(len(rcandidate))
		Tm = int(primer3.calcTm(rcandidate))
		if Tm > args.tmin and Tm < args.tmax:
			Hairpin = primer3.calcHairpin(rcandidate)
			if Hairpin.structure_found is False:
				outfile.write(title + "\tPost frame reverse\t" + rcandidate + "\t" + str(Tm) + "\t" + rcand_len + "\t" + str(args.f + framelen + (args.f - index)) + "\t" + str(framelen) + "\t" + str(orientation) + "\n")



outfile.close()

primertab = open(args.o, "r")

singleoutfile = open(args.os, "w")
multioutfile = open(args.om, "w")

primer_dict = defaultdict(int)

for line in primertab:
	columns = line.split("\t")
	primer = columns[2]
	primer_dict[primer] += 1

multiples = {}

for key, vals in primer_dict.items():
	if vals > 1:
		print(key,vals)
		multiples[key] = vals


primertab = open(args.input, "r")

for line in primertab:
	columns = line.split("\t")
	primer = columns[2]
	if primer in multiples.keys():
		multioutfile.write(line)
	else:
		singleoutfile.write(line)


multioutfile.close()
singleoutfile.close()




