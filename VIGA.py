#!/usr/bin/env python

# -*- coding: utf-8 -*-

# VIGA - De-novo VIral Genome Annotator
#
# Copyright (C) 2018 -  Enrique Gonzalez-Tortuero
#			Vimalkumar Velayudhan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Importing python libraries
from __future__ import print_function
import argparse
import BCBio.GFF
import csv
import fractions
import glob
import multiprocessing
import numpy
import os
import os.path
import re
import sys
import subprocess
import time
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import OrderedDict, defaultdict
from itertools import product
from scipy import signal
from time import strftime

## Defining the program version
version = "0.11.0"

## Preparing functions
# A batch iterator
def batch_iterator(iterator, batch_size):
	entry = True
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = iterator.next()
			except StopIteration:
				entry = None
			if entry is None:
				break
			batch.append(entry)
		if batch:
			yield batch

# Function equivalent to 'cat' command in UNIX systems (used to concatenate all files in a single one)
def cat_all(listinputfiles, outputfile):
	with open(outputfile, 'w') as finalfile:
		for fname in listinputfiles:
			with open(fname) as infile:
				for line in infile:
					finalfile.write(line)

# Function to verify if origin/terminus peaks are too close or too far apart (suggesting that these peaks are probably wrong)
def check_peaks(peaks, length): 
	closest, farthest = int(length * float(0.45)), int(length * float(0.55))
	pairs = []
	for pair in list(product(*peaks)):
		tr, pk = sorted(list(pair), key = lambda x: x[1], reverse = False) # trough and peak
		a = (tr[0] - pk[0]) % length
		b = (pk[0] - tr[0]) % length
		pt = abs(tr[1] - pk[1]) # distance between values
		if (a <= farthest and a >= closest) or (b <=farthest and b >= closest):
			pairs.append([pt, tr, pk]) # Added this to make sure gets origin and ter right
	if len(pairs) == 0:
		return [False, False]
	pt, tr, pk = sorted(pairs, reverse = True)[0]
	return [tr[0], pk[0]]

# Function equivalent to 'whereis'/'which' commands in UNIX systems (used to find a command in your computer)
def cmd_exists(cmd):
	return subprocess.call("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

# Function to calculate the GC skew
def GCskew(name, length, seq, window, slide):
	replacements = {'G':1, 'C':-1, 'A':0, 'T':0, 'N':0}
	gmc = [] # G - C
	for base in seq:
		try:
			gmc.append(replacements[base])
		except:
			gmc.append(0)
	gpc = [abs(i) for i in gmc] # convert to G + C
	weights = numpy.ones(window)/window # calculate sliding windows for (G - C) and (G + C)
	gmc = [[i, c] for i, c in enumerate(signal.fftconvolve(gmc, weights, 'same').tolist())] # calculate sliding windows for (G - C) and (G + C)
	gpc = [[i, c] for i, c in enumerate(signal.fftconvolve(gpc, weights, 'same').tolist())] # calculate sliding windows for (G - C) and (G + C)
	skew = [[], []] # x and y for gc skew 	# calculate gc skew and cummulative gc skew sum
	c_skew = [[], []] # x and y for gc skew cummulative sums 	# calculate gc skew and cummulative gc skew sum
	cs = 0 # cummulative sum 	# calculate gc skew and cummulative gc skew sum
	for i, m in gmc[0::slide]: # select windows to use based on slide
		p = gpc[i][1]
		if p == 0:
			gcs = 0
		else:
			gcs = m/p
		cs += gcs
		skew[0].append(i)
		c_skew[0].append(i)
		skew[1].append(gcs)
		c_skew[1].append(cs)
	ori, ter = find_ori_ter(c_skew, length)
	return ori, ter, skew, c_skew

# Function to print all messages in standard error instead of standard output
def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Function to find origin and terminus of replication based on cumulative GC skew min and max peaks
def find_ori_ter(c_skew, length): 
	c_skew_min = signal.argrelextrema(numpy.asarray(c_skew[1]), numpy.less, order = 1)[0].tolist()
	c_skew_max = signal.argrelextrema(numpy.asarray(c_skew[1]), numpy.greater, order = 1)[0].tolist()
	if len(c_skew_min) == 0 or len(c_skew_min) == 0:
		return [False, False] # return False if no peaks were detected
	else:
		c_skew_min = [[c_skew[0][i], c_skew[1][i]] for i in c_skew_min]
		c_skew_max = [[c_skew[0][i], c_skew[1][i]] for i in c_skew_max]
		ori, ter = check_peaks([c_skew_min, c_skew_max], length)
	return ori, ter

# Function to split a string numerically and not by alphabetic order
def stringSplitByNumbers(x):
	r = re.compile('(\d+)')
	l = r.split(x)
	return [int(y) if y.isdigit() else y for y in l]

## Processing the arguments
parser = argparse.ArgumentParser(description='VIGA is an automatic de novo VIral Genome Annotator.')
basic_group = parser.add_argument_group('Basic options for VIGA [REQUIRED]')
basic_group.add_argument("--input", dest="inputfile", type=str, required=True, help='Input file as a FASTA file', metavar="FASTAFILE")
basic_group.add_argument("--diamonddb", dest="diamonddatabase", type=str, required=True, help='DIAMOND Database that will be used for the protein function prediction. The database must be created from a amino acid FASTA file as indicated in https://github.com/bbuchfink/diamond.', metavar="DIAMONDDB")
basic_group.add_argument("--blastdb", dest="blastdatabase", type=str, required=True, help='BLAST Database that will be used to refine the protein function prediction in hypothetical proteins. The database must be an amino acid one, not nucleotidic.', metavar="BLASTDB")
basic_group.add_argument("--rfamdb", dest="rfamdatabase", type=str, required=True, help='RFAM Database that will be used for the ribosomal RNA prediction. RFAMDB should be in the format "/full/path/to/rfamdb/Rfam.cm" and must be compressed accordingly (see INFERNAL manual) before running the script.', metavar="RFAMDB")
basic_group.add_argument("--modifiers", dest="modifiers", type=str, required=True, help='Input file as a plain text file with the modifiers per every FASTA header according to SeqIn (https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html). All modifiers must be written in a single line and are separated by a single space character. No space should be placed besides the = sign. For example: [organism=Serratia marcescens subsp. marcescens] [sub-species=marcescens] [strain=AH0650_Sm1] [topology=linear] [moltype=DNA] [tech=wgs] [gcode=11] [country=Australia] [isolation-source=sputum]. This line will be copied and printed along with the record name as the definition line of every contig sequence.', metavar="TEXTFILE")

advanced_general_group = parser.add_argument_group('Advanced general options for VIGA [OPTIONAL]')
advanced_general_group.add_argument("--out", dest="rootoutput", type=str, help='Name of the outputs files (without extension)', metavar="OUTPUTNAME")
advanced_general_group.add_argument("--locus", dest="locus", type=str, default='LOC', help='Name of the sequences. If the input is a multiFASTA file, please put a general name as the program will add the number of the contig at the end of the name (Default: %(default)s)', metavar="STRING")
advanced_general_group.add_argument("--gff", dest="gffprint", action='store_true', default=False, help='Printing the output as GFF3 file (Default: False)')
advanced_general_group.add_argument("--threads", dest="ncpus", default=multiprocessing.cpu_count(), help='Number of threads/cpus (Default: %(default)s cpu [Max. number of cpus])', metavar="INT")
advanced_general_group.add_argument('--mincontigsize', dest="mincontigsize", type=int, default = 200, help = 'Minimum contig length to be considered in the final files (Default: 200 bp)', metavar="INT")

advanced_circularity_group = parser.add_argument_group('Advanced options for contig shape prediction in VIGA [OPTIONAL]')
advanced_circularity_group.add_argument("--readlength", dest="read_length", type=int, default=101, help='Read length for the circularity prediction (default: 101 bp)', metavar="INT")
advanced_circularity_group.add_argument("--windowsize", dest="window", type=int, default=100, help='Window size used to determine the origin of replication in circular contigs according to the cumulative GC skew (default: 100 bp)', metavar="INT")
advanced_circularity_group.add_argument("--slidingsize", dest="slide", type=int, default=10, help='Window size used to determine the origin of replication in circular contigs according to the cumulative GC skew (default: 10 bp)', metavar="INT")

advanced_repeats_group = parser.add_argument_group('Advanced options for detection of repeats in VIGA [OPTIONAL]')
advanced_repeats_group.add_argument("--minrepeat", dest="minrepeat", type=int, default=16, help="Minimum repeat length for CRISPR detection (Default: 16)", metavar="INT")
advanced_repeats_group.add_argument("--maxrepeat", dest="maxrepeat", type=int, default=64, help="Maximum repeat length for CRISPR detection (Default: 64)")
advanced_repeats_group.add_argument("--minspacer", dest="minspacer", type=int, default=8, help="Minimum spacer length for CRISPR detection (Default: 8)")
advanced_repeats_group.add_argument("--maxspacer", dest="maxspacer", type=int, default=64, help="Maximum spacer length for CRISPR detection (Default: 64)")

advanced_gcode_group = parser.add_argument_group('Advanced options for genetic code in VIGA [OPTIONAL]')
gcode_choices = {'1': 'Standard genetic code [Eukaryotic]', '2': 'Vertebrate mitochondrial code', '3': 'Yeast mitochondrial code', '4': 'Mycoplasma/Spiroplasma and Protozoan/mold/coelenterate mitochondrial code', '5': 'Invertebrate mitochondrial code', '6': 'Ciliate, dasycladacean and hexamita nuclear code', '9': 'Echinoderm/flatworm mitochondrial code', '10': 'Euplotid nuclear code', '11': 'Bacteria/Archaea/Phages/Plant plastid', '12': 'Alternative yeast nuclear code', '13': 'Ascidian mitochondrial code', '14': 'Alternative flatworm mitochondrial code', '16': 'Chlorophycean mitochondrial code', '21': 'Trematode mitochondrial code', '22': 'Scedenesmus obliquus mitochondrial code', '23': 'Thraustochytrium mitochondrial code', '24': 'Pterobranquia mitochondrial code', '25': 'Gracilibacteria & Candidate division SR1', '26': 'Pachysolen tannophilus nuclear code', '27': 'Karyorelict nuclear code', '28': 'Condylostoma nuclear code', '29': 'Mesodinium nuclear code', '30': 'Peritrich nuclear code', '31': 'Blastocrithidia nuclear code'}
gcode_help = ('Number of GenBank translation table. At this moment, the available options are {0}. (Default: %(default)s)'.format('{}'.format(', '.join('{0} ({1})'.format(k, v) for k, v in sorted(gcode_choices.items())))))
advanced_gcode_group.add_argument("--gcode", dest="gcode", type=str, default='11', help=gcode_help, metavar="NUMBER")

advanced_type_group = parser.add_argument_group('Advanced options for GenBank division in VIGA [OPTIONAL]')
type_choices = {'BCT': 'Prokaryotic chromosome', 'CON': 'Contig', 'PHG': 'Phages', 'VRL': 'Eukaryotic/Archaea virus'}
type_help = ('GenBank Division: One of the following codes - {0}. (Default: %(default)s)'.format(', '.join('{0} ({1})'.format(k, v) for k, v in type_choices.items())))
advanced_type_group.add_argument("--typedata", dest="typedata", type=str, default='CON', help=type_help, metavar="BCT|CON|VRL|PHG")

advanced_diamond_group = parser.add_argument_group('Advanced options for the first protein function prediction based on homology in VIGA [OPTIONAL]')
advanced_diamond_group.add_argument("--diamondevalue", dest="diamondevalue", default=0.00001, help='DIAMOND e-value threshold (Default: 0.00001)', metavar="FLOAT")
advanced_diamond_group.add_argument("--diamondidthr", dest="diamondidthreshold", default=50.00, help='DIAMOND ID threshold (Default: 50.0)', metavar="FLOAT")
advanced_diamond_group.add_argument("--diamondcoverthr", dest="diamondcovthreshold", default=50.00, help='DIAMOND Coverage threshold (Default: 50.0)', metavar="FLOAT")

advanced_blast_group = parser.add_argument_group('Advanced options for the second protein function prediction based on homology in VIGA [OPTIONAL]')
advanced_blast_group.add_argument("--blastexh", dest="blastexh", action='store_true', default=False, help='Use of exhaustive BLAST to predict the proteins by homology according to Fozo et al. (2010) Nucleic Acids Res (Default=FALSE)')
advanced_blast_group.add_argument("--blastevalue", dest="blastevalue", default=0.00001, help='BLAST e-value threshold (Default: 0.00001)', metavar="FLOAT")
advanced_blast_group.add_argument("--blastidthr", dest="blastidthreshold", default=50.00, help='BLAST ID threshold (Default: 50.0)', metavar="FLOAT")
advanced_blast_group.add_argument("--blastcoverthr", dest="blastcovthreshold", default=50.00, help='BLAST Coverage threshold (Default: 50.0)', metavar="FLOAT")

advanced_hmm_group = parser.add_argument_group('Advanced options for the third protein function prediction based on HMM in VIGA [OPTIONAL]')
advanced_hmm_group.add_argument("--nohmmer", dest="nohmmer", action='store_true', default=False, help='Running only BLAST to predict protein function. (Default: False)')
advanced_hmm_group.add_argument("--hmmerdb", dest="hmmerdatabase", type=str, help='HMMER Database that will be used to add additional information for all proteins according to Hidden Markov Models. This database must be in HMM format and it is mandatory if --nohmmer is disabled', metavar="HMMDB")
advanced_hmm_group.add_argument("--hmmerevalue", dest="hmmerevalue", default=0.001, help='HMMER e-value threshold (Default: 0.001)', metavar="FLOAT")
advanced_hmm_group.add_argument("--hmmeridthr", dest="hmmeridthreshold", default=50.00, help='HMMER ID threshold (Default: 50.0)', metavar="FLOAT")
advanced_hmm_group.add_argument("--hmmercoverthr", dest="hmmercovthreshold", default=50.00, help='HMMER Coverage threshold (Default: 50.0)', metavar="FLOAT")

args = parser.parse_args()

root_output = args.rootoutput
if not root_output:
	root_output = '{}_annotated'.format(os.path.splitext(args.inputfile)[0])

if args.nohmmer == False and args.hmmerdatabase == None:
	sys.exit('You MUST specify HMMER database using the parameter --hmmerdb if you are not using --nohmmer option')

## Printing the header of the program 
eprint("This is VIGA %s" % str(version))
eprint("Written by Enrique Gonzalez Tortuero & Vimalkumar Velayudhan")
eprint("Homepage is https://github.com/EGTortuero/viga")
eprint("Local time: ", strftime("%a, %d %b %Y %H:%M:%S"))

## checking the presence of the programs in the system
if not cmd_exists("lastz")==True:
	sys.exit("You must install LASTZ before using this script")
elif not cmd_exists("cmscan")==True:
	sys.exit("You must install INFERNAL before using this script")
elif not cmd_exists("aragorn")==True:
	sys.exit("You must install ARAGORN before using this script")
elif not cmd_exists("pilercr")==True:
	sys.exit("You must install PILERCR before using this script")
elif not cmd_exists("trf")==True:
	sys.exit("You must install Tandem Repeats Finder before using this script")
elif not cmd_exists("irf")==True:
	sys.exit("You must install Inverted Repeats Finder before using this script")
elif not cmd_exists("prodigal")==True:
	sys.exit("You must install prodigal before using this script")
elif not cmd_exists("diamond")==True:
	sys.exit("You must install DIAMOND before using this script")
elif not cmd_exists("blastp")==True:
	sys.exit("You must install BLAST before using this script")
elif not cmd_exists("hmmsearch")==True or not cmd_exists("hmmbuild"):
	sys.exit("You must install HMMER 3 before using this script")

eprint("Data type is {0} and GenBank translation table no is {1}\n".format(args.typedata, args.gcode))

## Correcting the original file (long headers problem + multiple FASTA files)
record_iter = SeqIO.parse(open(args.inputfile, "rU"),"fasta")
counter = 1
newnamessequences = {}
for i, batch in enumerate(batch_iterator(record_iter, 1)):
	seq_index = "CONTIG_%i" % (i + 1)
	with open("%s.temp.fna" % seq_index, "w") as handle:
		count = SeqIO.write(batch, handle, "fasta")
	with open("%s.temp.fna" % seq_index, "rU") as original, open("%s.fna" % seq_index, "w") as corrected:
		sequences = SeqIO.parse(original, "fasta", IUPAC.ambiguous_dna)
		for record in sequences:
			original_name = record.id
			record.id = "%s_%i" % (args.locus, counter)
			record.description = record.description
			counter += 1
			newnamessequences[record.id] = original_name
			eprint("WARNING: %s was renamed as %s" % (original_name, record.id))
		SeqIO.write(record, corrected, "fasta")
	with open("logfile.txt", "w") as logfile:
		logfile.write("#Original\tNew\n")
		for oldname in sorted(newnamessequences, key = stringSplitByNumbers):
			logfile.write("%s\t%s\n" % (newnamessequences[oldname], oldname))
	os.remove("%s.temp.fna" % seq_index)

## Predicting the shape of the contig (code based on Alex Crits-Christoph's find_circular.py script [https://github.com/alexcritschristoph/VICA/blob/master/find_circular.py])
starttime1 = time.time()
eprint("\nPredicting the shape for all contigs using LASTZ")
genomeshape = {}
with open("logfile.txt", "a") as logfile:
	logfile.write("\n#Contig\tShape\n")
	for newfile in sorted(glob.glob("CONTIG_*.fna")):
		with open(newfile, "rU") as targetfile:
			Sequence = SeqIO.parse(targetfile, "fasta")
			for record in Sequence:
				genomeshape[record.id] = {}
				seq_beginning = str(record.seq[0:args.read_length])
				seq_ending = str(record.seq[len(record.seq)-args.read_length:len(record.seq)])
				combined_seqs = SeqRecord(Seq(seq_beginning + seq_ending, IUPAC.ambiguous_dna), id = record.description)
				SeqIO.write(combined_seqs, "temporal_circular.fasta", "fasta")
				outputlastz = subprocess.check_output(["lastz", "temporal_circular.fasta", "--self", "--notrivial", "--nomirror", "--ambiguous=iupac", "--format=general-:start1,end1,start2,end2,score,strand1,strand2,identity,length1"])
				resultslastz = outputlastz.split("\n")
				for resultlastz in resultslastz:
					if resultlastz != '':
						start1 = resultlastz.split()[0]
						end1 = resultlastz.split()[1]
						start2 = resultlastz.split()[2]
						end2 = resultlastz.split()[3]
						strand1 = resultlastz.split()[5]
						strand2 = resultlastz.split()[6]
						identity = resultlastz.split()[7]
						length = int(resultlastz.split()[9])
						if strand1 == strand2 and length > 0.4 * args.read_length and float(fractions.Fraction(identity)) > 0.95 and int(start1) < 5 and int(start2) > args.read_length and int(end1) < args.read_length and int(end2) > args.read_length * 2 * 0.9:
							genomeshape[record.id]['genomeshape'] = "circular"
							try:
								if genomeshape[record.id]['identity'] >= float(fractions.Fraction(identity)):
									genomeshape[record.id]['identity'] = float(fractions.Fraction(identity))
									genomeshape[record.id]['length'] = length
							except KeyError:
								genomeshape[record.id]['identity'] = float(fractions.Fraction(identity))
								genomeshape[record.id]['length'] = length
						else:
							continue
						if strand1 == strand2 and length > 0.4 * args.read_length and float(fractions.Fraction(identity)) > 0.95 and int(start1) < 5 and int(start2) > args.read_length and int(end1) < args.read_length and int(end2) > args.read_length * 2 * 0.9:
							genomeshape[record.id]['genomeshape'] = "circular"
							try:
								if genomeshape[record.id]['identity'] >= float(fractions.Fraction(identity)):
									genomeshape[record.id]['identity'] = float(fractions.Fraction(identity))
									genomeshape[record.id]['length'] = length
							except KeyError:
								genomeshape[record.id]['identity'] = float(fractions.Fraction(identity))
								genomeshape[record.id]['length'] = length
				try:
					if genomeshape[record.id]['genomeshape'] == "":
							genomeshape[record.id]['genomeshape'] = "linear"
				except KeyError:
					genomeshape[record.id]['genomeshape'] = "linear"
				else:
					genomeshape[record.id]['genomeshape'] = "circular"
					with open("temp.fasta", "w") as correctedcircular:
						Corrseq = str(record.seq[int(genomeshape[record.id]['length'])//2:-int(genomeshape[record.id]['length'])//2])
						Newseq = SeqRecord(Seq(Corrseq, IUPAC.ambiguous_dna), id = record.description)
						SeqIO.write(Newseq, correctedcircular, "fasta")
					os.rename("temp.fasta", "temp2.fasta")
				eprint("%s seems to be a %s contig according to LASTZ" % (record.id, genomeshape[record.id]['genomeshape']))
				logfile.write("%s\t%s\n" % (record.id, genomeshape[record.id]['genomeshape']))
				
## Calculate the cumulative GCskew in circular contigs to determine the origin of replication (Based on iRep -- Brown CT, Olm MR, Thomas BC, Banfield JF (2016) Measurement of bacterial replication rates in microbial communities. Nature Biotechnology 34: 1256-63.)
	if genomeshape[record.id]['genomeshape'] == "circular":
		eprint("Determining the origin of replication of %s according to the GC Skew" % newfile)
		for gotocircularize in SeqIO.parse("temp2.fasta", "fasta"):
			length_contig = len(gotocircularize.seq)
			oric, term, gcskew, cgcskew = GCskew(gotocircularize.id, length_contig, gotocircularize.seq, args.window, args.slide)
			valoric = oric
			if oric == False:
				oric, term = 'n/a', 'n/a'
			else:
				oric, term = '{:,}'.format(oric), '{:,}'.format(term)
			eprint('%s -> Origin: %s Terminus: %s' % (gotocircularize.id, oric, term))
			if valoric != False:
				firstpartseq = str(gotocircularize.seq[valoric:-1])
				secondpartseq = str(gotocircularize.seq[0:(valoric-1)])
				combinedcorrectedseq = SeqRecord(Seq(firstpartseq + secondpartseq, IUPAC.ambiguous_dna), id = gotocircularize.description)
				SeqIO.write(combinedcorrectedseq, newfile, "fasta")
			else:
				eprint("VIGA was unable to predict the origin of replication: %s was not modified!" % record.id)
				os.rename("temp2.fasta", newfile)
		if os.path.isfile("temp2.fasta"):
			os.remove("temp2.fasta")
cat_all(sorted(glob.glob('CONTIG_*.fna')), 'CONTIGS_ALL.fasta')
endtime1 = time.time()
durationtime1 = endtime1 - starttime1
eprint("Done: shape prediction took %s seconds" % str(durationtime1))

## Predicting the rRNA sequences
starttime2 = time.time()
subunits = dict()
countLSU = 0
countSSU = 0
count5S = 0
eprint("\nAnalyzing the presence of host contamination according to the presence of rRNA")
with open("/dev/null", "w") as stderr:
	subprocess.call(["cmscan", "--rfam", "--cut_ga", "--nohmmonly", "--tblout", "rrnafile.csv", "--cpu", str(args.ncpus), args.rfamdatabase, "CONTIGS_ALL.fasta"], stdout=stderr)
with open("rrnafile.csv", "rU") as rrnafile:
	namedict = {"SSU_rRNA_archaea": "16s_rRNA", "SSU_rRNA_bacteria": "16s_rRNA", "SSU_rRNA_eukarya": "18s_rRNA", "SSU_rRNA_microsporidia": "16s_rRNA", "LSU_rRNA_archaea": "23s_rRNA", "LSU_rRNA_bacteria": "23s_rRNA", "LSU_rRNA_eukarya": "28s_rRNA", "5S_rRNA": "5s_rRNA"}
	elementsrRNA = {}
	for line in rrnafile:
		if not line.startswith("#"):
			InfoLINE = re.sub("\s+", ",", line)
			line_splitted = InfoLINE.split(",")
			contig_id = line_splitted[2]
			elementsrRNA[contig_id] = {'LSU': {}, 'SSU': {}, '5S': {}}
			subunits[contig_id] = {}
with open("rrnafile.csv", "rU") as rrnafile:
	for line in rrnafile:
		if not line.startswith("#"):
			InfoLINE = re.sub("\s+", ",", line)
			line_splitted = InfoLINE.split(",")
			item_type = line_splitted[0]
			contig_id = line_splitted[2]
			for saved_contig in elementsrRNA:
				if saved_contig == contig_id:
					if item_type.startswith(('LSU', 'SSU', '5S')):
						if item_type.startswith("LSU"):
							try:
								if elementsrRNA[contig_id]['LSU']['bigscore'] <= float(line_splitted[14]):
									elementsrRNA[contig_id]['LSU'] = {'type': item_type, 'bigscore': line_splitted[14]}
							except KeyError:
								elementsrRNA[contig_id]['LSU'] = {'type': item_type, 'bigscore': line_splitted[14]}
								subunits[contig_id]['LSU'] = {}
						elif item_type.startswith("SSU"):
							try:
								if elementsrRNA[contig_id]['SSU']['bigscore'] <= float(line_splitted[14]):
									elementsrRNA[contig_id]['SSU'] = {'type': item_type, 'bigscore': line_splitted[14]}
							except KeyError:
								elementsrRNA[contig_id]['SSU'] = {'type': item_type, 'bigscore': line_splitted[14]}
								subunits[contig_id]['SSU'] = {}
						else:
							try:
								if elementsrRNA[contig_id]['5S']['bigscore'] <= float(line_splitted[14]):
									elementsrRNA[contig_id]['5S'] = {'type': item_type, 'bigscore': line_splitted[14]}
							except KeyError:
								elementsrRNA[contig_id]['5S'] = {'type': item_type, 'bigscore': line_splitted[14]}
								subunits[contig_id]['5S'] = {}
with open("rrnafile.csv", "rU") as rrnafile:
	for line in rrnafile:
		if not line.startswith("#"):
			InfoLINE = re.sub("\s+", ",", line)
			line_splitted = InfoLINE.split(",")
			item_type = line_splitted[0]
			contig_id = line_splitted[2]
			for saved_contig in elementsrRNA:
				if saved_contig == contig_id:
					if item_type.startswith("LSU"):
						if item_type == elementsrRNA[contig_id]['LSU']['type']:
							countLSU += 1
					elif item_type.startswith("SSU"):
						if item_type == elementsrRNA[contig_id]['SSU']['type']:
							countSSU += 1
					elif item_type.startswith("5S"):
						if item_type == elementsrRNA[contig_id]['5S']['type']:
							count5S += 1
for contig in subunits:
	for element in subunits[contig]:
		i = 0
		j = 0
		k = 0
		if element == 'LSU':
			while i < countLSU:
				subunits[contig][element][i] = {}
				i += 1
		elif element == 'SSU':
			while j < countSSU:
				subunits[contig][element][j] = {}
				j += 1
		elif element == '5S':
			while k < count5S:
				subunits[contig][element][k] = {}
				k += 1
with open("rrnafile.csv", "rU") as rrnafile:
	a = 0
	b = 0
	c = 0
	for line in rrnafile:
		if not line.startswith("#"):
			InfoLINE = re.sub("\s+", ",", line)
			line_splitted = InfoLINE.split(",")
			item_type = line_splitted[0]
			contig_id = line_splitted[2]
			for saved_contig in elementsrRNA:
				if saved_contig == contig_id:
					if item_type.startswith("LSU"):
						if item_type == elementsrRNA[contig_id]['LSU']['type']:
							if line_splitted[9] == "+":
								subunits[saved_contig]['LSU'][a] = {'type': item_type, 'product': namedict[item_type], 'begin': int(line_splitted[7]), 'end': int(line_splitted[8]), 'score': line_splitted[14], 'strand': 1}
							else:
								subunits[saved_contig]['LSU'][a] = {'type': item_type, 'product': namedict[item_type], 'begin': int(line_splitted[8]), 'end': int(line_splitted[7]), 'score': line_splitted[14], 'strand': -1}	
						a += 1
					elif item_type.startswith("SSU"):
						if item_type == elementsrRNA[contig_id]['SSU']['type']:
							if line_splitted[9] == "+":
								subunits[saved_contig]['SSU'][b] = {'type': item_type, 'product': namedict[item_type], 'begin': int(line_splitted[7]), 'end': int(line_splitted[8]), 'score': line_splitted[14], 'strand': 1}
							else:
								subunits[saved_contig]['SSU'][b] = {'type': item_type, 'product': namedict[item_type], 'begin': int(line_splitted[8]), 'end': int(line_splitted[7]), 'score': line_splitted[14], 'strand': -1}	
						b += 1
					elif item_type.startswith("5S"):
						if item_type == elementsrRNA[contig_id]['5S']['type']:
							if line_splitted[9] == "+":
								subunits[saved_contig]['5S'][c] = {'type': item_type, 'product': namedict[item_type], 'begin': int(line_splitted[7]), 'end': int(line_splitted[8]), 'score': line_splitted[14], 'strand': 1}
							else:
								subunits[saved_contig]['5S'][c] = {'type': item_type, 'product': namedict[item_type], 'begin': int(line_splitted[8]), 'end': int(line_splitted[7]), 'score': line_splitted[14], 'strand': -1}	
						c += 1
with open("logfile.txt", "a") as logfile:
	logfile.write("\n#Contig\trRNA product\tStart\tEnd\n")
	for suspectcontig in subunits:
		for rRNA in subunits[suspectcontig]:
			count = 0
			try:
				lengthlist = len(subunits[suspectcontig][rRNA])
			except KeyError:
				continue
			else:
				while count < lengthlist:
					eprint("WARNING: %s: %s (%i ... %i) -> Potential host contamination!" % (suspectcontig, rRNA, int(subunits[suspectcontig][rRNA][count]['begin']), int(subunits[suspectcontig][rRNA][count]['end'])))
					logfile.write("%s\t%s\t%i\t%i\n" % (suspectcontig, subunits[suspectcontig][rRNA][count]['product'], int(subunits[suspectcontig][rRNA][count]['begin']), int(subunits[suspectcontig][rRNA][count]['end'])))
					count +=1
endtime2 = time.time()
durationtime2 = endtime2 - starttime2
eprint("Done: rRNA detection took %s seconds" % str(durationtime2))

## Predicting the tRNA sequences using ARAGORN
starttime3 = time.time()
tRNAdict = dict()
tmRNAdict = dict()
genetictable = "-gc%s" % str(args.gcode)
eprint("\nRunning ARAGORN to predict tRNA-like sequences for all contigs")
for contigfile in sorted(glob.glob("CONTIG_*.fna")):
	with open(contigfile, "rU") as targetfasta:
		for newrecord in SeqIO.parse(targetfasta, "fasta"):
			putativetrnafile = "trnafile_%s.fasta" % newrecord.id
			with open(putativetrnafile, "w") as trnafile:
				if genomeshape[newrecord.id]['genomeshape'] == "circular":
					subprocess.call(["aragorn", "-c", "-fon", genetictable, contigfile], stdout=trnafile)
				else:
					subprocess.call(["aragorn", "-l", "-fon", genetictable, contigfile], stdout=trnafile)
			num_tRNA = len(list(SeqIO.parse(putativetrnafile, "fasta")))
			eprint("Detected %i tRNAs in %s" % (num_tRNA, newrecord.id))
			with open(putativetrnafile, "rU") as trnafile:
				tRNAdict[newrecord.id] = {}
				tmRNAdict[newrecord.id] = {}
				for tRNAseq in SeqIO.parse(trnafile, "fasta"):
					indtRNA = {}
					indtmRNA = {}
					tRNA_information = tRNAseq.description.split(" ")
					tRNApat = re.compile("^tRNA-")
					if tRNA_information[1] == "tmRNA":
						if str(tRNA_information[2]) == "(Permuted)":
							indtmRNA['product'] = "tmRNA"
							tmRNA_coords = str(tRNA_information[3])
							Beginningrevcomppat = re.compile("^c")
							if re.match(Beginningrevcomppat, tmRNA_coords):
								indtmRNA['strand'] = -1
								tmRNA_coords = tmRNA_coords.replace("c[","").replace("]","").split(",")
							else:
								indtmRNA['strand'] = 1
								tmRNA_coords = tmRNA_coords.replace("[","").replace("]","").split(",")
							indtmRNA['begin'] = int(tmRNA_coords[0])
							indtmRNA['end'] = int(tmRNA_coords[1])
							tmRNAdict[newrecord.id][tRNAseq.id] = indtmRNA
						else:
							indtmRNA['product'] = "tmRNA"
							tmRNA_coords = str(tRNA_information[2])
							Beginningrevcomppat = re.compile("^c")
							if re.match(Beginningrevcomppat, tmRNA_coords):
								indtmRNA['strand'] = -1
								tmRNA_coords = tmRNA_coords.replace("c[","").replace("]","").split(",")
							else:
								indtmRNA['strand'] = 1
								tmRNA_coords = tmRNA_coords.replace("[","").replace("]","").split(",")
							indtmRNA['begin'] = int(tmRNA_coords[0])
							indtmRNA['end'] = int(tmRNA_coords[1])
							tmRNAdict[newrecord.id][tRNAseq.id] = indtmRNA
					elif re.match(tRNApat, tRNA_information[1]):
						indtRNA['product'] = re.sub("\(\w{3}\)", "",  tRNA_information[1])
						tRNA_coords = tRNA_information[2]
						Beginningrevcomppat = re.compile("^c")
						if re.match(Beginningrevcomppat, tRNA_coords):
							indtRNA['strand'] = -1
							tRNA_coords = tRNA_coords.replace("c[","").replace("]","").split(",")
						else:
							indtRNA['strand'] = 1
							tRNA_coords = tRNA_coords.replace("[","").replace("]","").split(",")
						indtRNA['begin'] = int(tRNA_coords[0])
						indtRNA['end'] = int(tRNA_coords[1])
						tRNAdict[newrecord.id][tRNAseq.id] = indtRNA
endtime3 = time.time()
durationtime3 = endtime3 - starttime3
eprint("Done: tRNA and tmRNA detection took %s seconds" % str(durationtime3))

## Predicting CRISPR repeats
starttime4 = time.time()
eprint("\nRunning PILER-CR to predict CRISPR repeats for all contigs")
with open("/dev/null", "w") as apocalypse:
	subprocess.call(["pilercr", "-in", "CONTIGS_ALL.fasta", "-out", "crisprfile.txt", "-noinfo", "-minrepeat", str(args.minrepeat), "-maxrepeat", str(args.maxrepeat), "-minspacer", str(args.minspacer), "-maxspacer", str(args.maxspacer)], stderr=apocalypse)
information_CRISPR = {}
with open("crisprfile.txt", "rU") as crisprfile:
	for line in crisprfile:
		if "SUMMARY BY POSITION" in line:
			for line in crisprfile:
				try:
					patC1 = re.compile('^\s+\d+\s+(.{16})\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d?\s+\w+')
					sequence_id = re.split("\s", re.match(patC1, line).groups()[0])[0]
					information_CRISPR[sequence_id] = {}
				except AttributeError:
					continue
with open("crisprfile.txt", "rU") as crisprfile:
	for line in crisprfile:
		if "SUMMARY BY POSITION" in line:
			for line in crisprfile:
				information_crispr_repeat = {}
				try:
					patC = re.compile('^\s+(\d+)\s+(.{16})\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+\d+\s+\d?\s+(\w+)')
					key, sequence_id, start, length, seq = re.match(patC, line).groups()
					sequence_id = re.split("\s", sequence_id)[0]
				except AttributeError:
					continue
				else:
					information_crispr_repeat['start'] = start
					information_crispr_repeat['end'] = int(start) + int(length)
					information_crispr_repeat['repeatseq'] = seq
					information_crispr_repeat['repeatend'] = int(start) + len(seq)
					information_CRISPR[sequence_id][key] = information_crispr_repeat
for seqid in SeqIO.parse("CONTIGS_ALL.fasta", "fasta"):
	try:
		num_CRISPRs = len(information_CRISPR[seqid.id])
	except KeyError:
		eprint("Detected 0 CRISPR repeats in %s" % seqid.id)
	else:
		eprint("Detected %i CRISPR repeats in %s" % (num_CRISPRs, seqid.id))
endtime4 = time.time()
durationtime4 = endtime4 - starttime4
eprint("Done: CRISPR repeats detection took %s seconds" % str(durationtime4))

## Predicting tandem repeats using TRF
starttime5 = time.time()
eprint("\nPredicting tandem repeats in the sequences using TRF")
with open("/dev/null", "w") as stderr:
	subprocess.call(["trf", "CONTIGS_ALL.fasta", "2", "7", "7", "80", "10", "50", "500", "-h"], stdout=stderr, stderr=stderr)
	os.rename("CONTIGS_ALL.fasta.2.7.7.80.10.50.500.dat", "trf_temp.dat")
information_TRF = {}
with open("trf_temp.dat", "rU") as trfile:
	for line in trfile:
		information_tandem_repeat = {}	
		try:
			patT1 = re.compile('^Sequence:\s+(.+)')
			seqtid = re.split(" ", re.match(patT1, line).groups()[0])[0]
		except AttributeError:
			try:
				patT = re.compile('^(\d+)\s(\d+)\s\d+\s')
				start, end = re.match(patT, line).groups()
			except AttributeError:
				continue
			else:
				information_tandem_repeat['start'] = start
				information_tandem_repeat['end'] = end
				information_TRF[seqtid][count] = information_tandem_repeat
				count += 1
		else:
			information_TRF[seqtid] = {}
			count = 1
for contigid in sorted(information_TRF):
	number_TRs = len(information_TRF[contigid])
	eprint("Detected %i tandem repeats in %s" % (number_TRs, contigid))
endtime5 = time.time()
durationtime5 = endtime5 - starttime5
eprint("Done: tandem repeats detection took %s seconds" % str(durationtime5))

## Predicting inverted repeats using IRF
starttime6 = time.time()
eprint("\nPredicting inverted repeats in the sequences using IRF")
with open("/dev/null", "w") as stderr:
	subprocess.call(["irf", "CONTIGS_ALL.fasta", "2", "3", "5", "80", "10", "40", "500000", "10000", "-d", "-h"], stdout=stderr, stderr=stderr)
	os.rename("CONTIGS_ALL.fasta.2.3.5.80.10.40.500000.10000.dat", "irf_temp.dat")
information_IRF = {}
with open("irf_temp.dat", "rU") as irfile:
	for line in irfile:
		information_inverted_repeat = {}
		try:
			patI1 = re.compile('^Sequence:\s+(.+)')
			seqiid = re.split(" ", re.match(patI1, line).groups()[0])[0]
		except AttributeError:
			try:
				patI = re.compile('^(\d+)\s(\d+)\s\d+\s')
				start, end = re.match(patI, line).groups()
			except AttributeError:
				continue
			else:
				information_inverted_repeat['start'] = start
				information_inverted_repeat['end'] = end
				information_IRF[seqiid][count] = information_inverted_repeat
				count += 1
		else:
			information_IRF[seqiid] = {}
			count = 1
for contigid in sorted(information_IRF):
	number_IRs = len(information_IRF[contigid])
	eprint("Detected %i inverted repeats in %s" % (number_IRs, contigid))
endtime6 = time.time()
durationtime6 = endtime6 - starttime6
eprint("Done: inverted repeats detection took %s seconds" % str(durationtime6))

## Predicting genes using PRODIGAL
starttime7 = time.time()
eprint("\nRunning Prodigal to predict the ORFs in all contigs")
for contigfile in sorted(glob.glob("CONTIG_*.fna")):
	for record in SeqIO.parse(contigfile, "fasta"):
		orffile = "orffile_%s.faa" % record.id
		length_contig = len(record.seq)
		if (length_contig >= 100000):
			if genomeshape[record.id]['genomeshape'] == 'linear':
				subprocess.call(["prodigal", "-a", "pretemp.faa", "-i", contigfile, "-o", "/dev/null", "-g", args.gcode, "-c", "-q"])
			else:
				subprocess.call(["prodigal", "-a", "pretemp.faa", "-i", contigfile, "-o", "/dev/null", "-g", args.gcode, "-q"])
		else:
			if genomeshape[record.id]['genomeshape'] == 'linear':
				subprocess.call(["prodigal", "-a", "pretemp.faa", "-i", contigfile, "-o", "/dev/null", "-p", "meta", "-g", args.gcode, "-c", "-q"])
			else:
				subprocess.call(["prodigal", "-a", "pretemp.faa", "-i", contigfile, "-o", "/dev/null", "-p", "meta", "-g", args.gcode, "-q"])
		with open("pretemp.faa", "rU") as originalfaa, open(orffile, "w") as correctedfaa:
			sequences = SeqIO.parse(originalfaa, "fasta")
			for record in sequences:
				record.seq = record.seq.rstrip("*")
				SeqIO.write(record, correctedfaa, "fasta")
		num_seqs = len(list(SeqIO.parse("pretemp.faa", "fasta")))
		eprint("Detected %i ORFs in %s" % (num_seqs, contigfile))
cat_all(sorted(glob.glob('orffile*.faa')), 'PROTS_FIRST_ROUND.faa')
endtime7 = time.time()
durationtime7 = endtime7 - starttime7
eprint("Done: protein prediction took %s seconds" % str(durationtime7))
 
## Predicting protein function based on homology using DIAMOND
starttime8 = time.time()
eprint("\nRunning DIAMOND to predict the protein function according to homology inference using default parameters")
with open("PROTS_FIRST_ROUND.faa", "r") as inputstep:
	first_line = inputstep.readline()
	if first_line.startswith(">"):
		subprocess.call(['diamond', 'blastp', '-q', 'PROTS_FIRST_ROUND.faa', '-d', args.diamonddatabase, '-e', str(args.diamondevalue), '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qstart', 'qend', 'evalue', 'bitscore', 'stitle', '-o', 'PROTS_FIRST_ROUND.csv', '-k', '10', "-p", str(args.ncpus), '--quiet'])
	else:
		open("PROTS_FIRST_ROUND.csv", 'a').close()
hypotheticalpat = re.compile(r'(((((?i)hypothetical)|(?i)uncharacteri(z|s)ed|(?i)predicted))( phage)?( membrane)? protein)|((?i)ORF|((?i)unnamed protein product|(?i)gp\d+|protein of unknown function|phage protein))')
with open("PROTS_FIRST_ROUND.csv", "rU") as diamondresults:
	reader = csv.DictReader(diamondresults, delimiter='\t', fieldnames=['qseqid','sseqid','pident','length','qlen','slen','qstart','qend','evalue','bitscore','stitle'])
	information_proteins_diamond = {}
	for row in reader:
		perc_cover = round(100.00*(float(row['length'])/float(row['qlen'])),2)
		perc_id = float(row['pident'])
		infoprot_diamond = {}
		infoprot_diamond['sseqid'] = row['sseqid']
		infoprot_diamond['pident'] = perc_id
		infoprot_diamond['pcover'] = perc_cover
		infoprot_diamond['evalue'] = row['evalue']
		infoprot_diamond['descr'] = row['stitle']
		try:
			data = information_proteins_diamond[row['qseqid']]
		except KeyError:
			if not re.search(hypotheticalpat, infoprot_diamond['descr']) and float(perc_id) >= float(args.diamondidthreshold) and float(perc_cover) >= float(args.diamondcovthreshold) and float(row['evalue']) <= float(args.diamondevalue):
				information_proteins_diamond[row['qseqid']] = infoprot_diamond
			else:
				continue
		else:
			if not re.search(hypotheticalpat, infoprot_diamond['descr']) and float(perc_id) >= float(args.diamondidthreshold) and float(perc_id) >= float(infoprot_diamond['pident']) and float(perc_cover) >= float(args.diamondcovthreshold) and float(perc_cover) >= float(infoprot_diamond['pcover']) and float(row['evalue']) <= float(args.diamondevalue):
				information_proteins_diamond[row['qseqid']] = infoprot_diamond
records_in_memory = list(SeqIO.parse(open("PROTS_FIRST_ROUND.faa", "r"),"fasta"))
i = 0
while i < len(records_in_memory):
	dataprot = records_in_memory[i].description.split(' # ')
	try:
		if information_proteins_diamond[dataprot[0]]['descr'] != None:
			pass
		elif re.search(hypotheticalpat, information_proteins_diamond[dataprot[0]]['descr']):
			SeqIO.write(records_in_memory[i], open("PROTS_SECOND_ROUND.faa", "a"), "fasta")
		else:
			SeqIO.write(records_in_memory[i], open("PROTS_SECOND_ROUND.faa", "a"), "fasta")
	except KeyError:
		SeqIO.write(records_in_memory[i], open("PROTS_SECOND_ROUND.faa", "a"), "fasta")
	i += 1
endtime8 = time.time()
durationtime8 = endtime8 - starttime8
eprint("Done: function prediction based on DIAMOND took %s seconds" % str(durationtime8))

## Running BLAST to refine the gene annotation
starttime9 = time.time()
with open("PROTS_SECOND_ROUND.faa", "r") as secondinputstep:
	first_line = secondinputstep.readline()
	if first_line.startswith(">") and args.blastexh==True:
		eprint("\nRunning BLAST to predict the genes according to homology inference using exhaustive mode (see Fozo et al. (2010) Nucleic Acids Res for details)")
		subprocess.call(['blastp', '-query', "PROTS_SECOND_ROUND.faa", '-db', args.blastdatabase, '-evalue', str(args.blastevalue), '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle', '-out', 'PROTS_SECOND_ROUND.csv', '-max_target_seqs', '10', '-word_size', '2', '-gapopen', '8', '-gapextend', '2', '-matrix', 'PAM70', '-comp_based_stats', '"0"', "-num_threads", str(args.ncpus)])
	elif first_line.startswith(">") and args.blastexh==False:
		eprint("\nRunning BLAST to predict the genes according to homology inference using default parameters")
		subprocess.call(['blastp', '-query', "PROTS_SECOND_ROUND.faa", '-db', args.blastdatabase, '-evalue', str(args.blastevalue), '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle', '-out', 'PROTS_SECOND_ROUND.csv', '-max_target_seqs', '10', "-num_threads", str(args.ncpus)])
	else:
		open("PROTS_SECOND_ROUND.csv", 'a').close()
with open("PROTS_SECOND_ROUND.csv", "rU") as blastresults:
	reader = csv.DictReader(blastresults, delimiter='\t', fieldnames=['qseqid','sseqid','pident','length','qlen','slen','qstart','qend','evalue','bitscore','stitle'])
	information_proteins_blast = {}
	for row in reader:
		perc_cover = round(100.00*(float(row['length'])/float(row['qlen'])),2)
		perc_id = float(row['pident'])
		infoprot_blast = {}
		infoprot_blast['sseqid'] = row['sseqid']
		infoprot_blast['pident'] = perc_id
		infoprot_blast['pcover'] = perc_cover
		infoprot_blast['evalue'] = row['evalue']
		infoprot_blast['descr'] = row['stitle']
		try:
			data = information_proteins_blast[row['qseqid']]
		except KeyError:
			if not re.search(hypotheticalpat, infoprot_blast['descr']) and float(perc_id) >= float(args.blastidthreshold) and float(perc_cover) >= float(args.blastcovthreshold) and float(row['evalue']) <= float(args.blastevalue):
				information_proteins_blast[row['qseqid']] = infoprot_blast
			else:
				continue
		else:
			if not re.search(hypotheticalpat, infoprot_blast['descr']) and float(perc_id) >= float(args.blastidthreshold) and float(perc_id) >= float(infoprot_blast['pident']) and float(perc_cover) >= float(args.blastcovthreshold) and float(perc_cover) >= float(infoprot_blast['pcover']) and float(row['evalue']) <= float(args.blastevalue):
				information_proteins_blast[row['qseqid']] = infoprot_blast

## Defining the proteins and preparing the file for the third round (if activated)
records_in_memory = list(SeqIO.parse(open("PROTS_FIRST_ROUND.faa", "r"),"fasta"))
protsdict = {}
i = 0
namelocuspat = re.compile(r'(\S+)\_\d+')
while i < len(records_in_memory):
	dataprot = records_in_memory[i].description.split(' # ')
	namelocus = re.match(namelocuspat, dataprot[0]).groups()[0]
	protsdict[namelocus] = {}
	i += 1
j = 0
while j < len(records_in_memory):
	dataprot = records_in_memory[j].description.split(' # ')
	namelocus = re.match(namelocuspat, dataprot[0]).groups()[0]
	modseq = str(records_in_memory[j].seq).replace("X","")
	analysed_seq = ProteinAnalysis(modseq)
	protsdict[namelocus][dataprot[0]] = {}
	protsdict[namelocus][dataprot[0]]['length'] = len(records_in_memory[j].seq)
	protsdict[namelocus][dataprot[0]]['isoelectricpoint'] = analysed_seq.isoelectric_point()
	protsdict[namelocus][dataprot[0]]['molweightkda'] = analysed_seq.molecular_weight()/1000.00
	protsdict[namelocus][dataprot[0]]['instability'] = analysed_seq.instability_index()
	protsdict[namelocus][dataprot[0]]['protein_id'] = dataprot[0]
	protsdict[namelocus][dataprot[0]]['strand'] = int(dataprot[3])
	protsdict[namelocus][dataprot[0]]['begin'] = int(dataprot[1])-1
	protsdict[namelocus][dataprot[0]]['end'] = int(dataprot[2])
	protsdict[namelocus][dataprot[0]]['translation'] = records_in_memory[j].seq
	try:
		if information_proteins_diamond[dataprot[0]]['descr'] == None and information_proteins_blast[dataprot[0]]['descr'] == None:
			protsdict[namelocus][dataprot[0]]['descr'] = 'Hypothetical protein'
			protsdict[namelocus][dataprot[0]]['source'] = "NO_HIT"
			protsdict[namelocus][dataprot[0]]['pident'] = "NA"
			protsdict[namelocus][dataprot[0]]['pcover'] = "NA"
			protsdict[namelocus][dataprot[0]]['evalue'] = "NA"
		elif information_proteins_diamond[dataprot[0]]['descr'] == None and re.search(Hypotheticalpat, information_proteins_blast[dataprot[0]]['descr']):
			protsdict[namelocus][dataprot[0]]['descr'] = 'Hypothetical protein'
			protsdict[namelocus][dataprot[0]]['source'] = "NO_HIT"
			protsdict[namelocus][dataprot[0]]['pident'] = "NA"
			protsdict[namelocus][dataprot[0]]['pcover'] = "NA"
			protsdict[namelocus][dataprot[0]]['evalue'] = "NA"
		elif information_proteins_diamond[dataprot[0]]['descr'] == None:
			listdescr = information_proteins_blast[dataprot[0]]['descr'].split(" ")
			del listdescr[0]
			protsdict[namelocus][dataprot[0]]['descr'] = " ".join(listdescr)			
			protsdict[namelocus][dataprot[0]]['source'] = "BLAST"
			protsdict[namelocus][dataprot[0]]['pident'] = information_proteins_blast[dataprot[0]]['pident']
			protsdict[namelocus][dataprot[0]]['pcover'] = information_proteins_blast[dataprot[0]]['pcover']
			protsdict[namelocus][dataprot[0]]['evalue'] = information_proteins_blast[dataprot[0]]['evalue']
		elif re.search(hypotheticalpat, information_proteins_diamond[dataprot[0]]['descr']):
			protsdict[namelocus][dataprot[0]]['descr'] = 'Conserved hypothetical protein'
		else:
			listdescr = information_proteins_diamond[dataprot[0]]['descr'].split(" ")
			del listdescr[0]
			protsdict[namelocus][dataprot[0]]['descr'] = " ".join(listdescr)
			protsdict[namelocus][dataprot[0]]['source'] = "DIAMOND"
			protsdict[namelocus][dataprot[0]]['pident'] = information_proteins_diamond[dataprot[0]]['pident']
			protsdict[namelocus][dataprot[0]]['pcover'] = information_proteins_diamond[dataprot[0]]['pcover']
			protsdict[namelocus][dataprot[0]]['evalue'] = information_proteins_diamond[dataprot[0]]['evalue']
	except KeyError:
		protsdict[namelocus][dataprot[0]]['descr'] = 'Hypothetical protein'
		protsdict[namelocus][dataprot[0]]['source'] = "NO_HIT"
		protsdict[namelocus][dataprot[0]]['pident'] = "NA"
		protsdict[namelocus][dataprot[0]]['pcover'] = "NA"
		protsdict[namelocus][dataprot[0]]['evalue'] = "NA"
	j += 1
endtime9 = time.time()
durationtime9 = endtime9 - starttime9
eprint("Done: function prediction based on BLAST took %s seconds" % str(durationtime9))

## Predicting the function of the proteins based on HMM predictions using HMMer 3.0 (Slow step)
if args.nohmmer == False:
	starttime10 = time.time()
	eprint("\nRunning HMMER to enrich the annotations for all viral proteins using PVOGs")
	subprocess.call(['hmmsearch', '--cpu', str(args.ncpus), '--domtblout', 'PROTS_FIRST_ROUND.tbl', '-E', str(args.hmmerevalue), '-o', '/dev/null', args.hmmerdatabase, 'PROTS_FIRST_ROUND.faa'])
	information_proteins_hmmer = {}
	with open("PROTS_FIRST_ROUND.tbl", "rU") as tblfile:
		for line in tblfile:
			if line.startswith("#"):
				continue
			else:
				infoprot_hmmer = {}
				try:
					pat = re.compile("^(\S+)\s+\S\s+\d+\s+(\S+)\s+\S\s+(\d+)\s+((?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?)\s+\S+\s+\S+\s+\S+\s+\S+\s+(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+\S+\s+\S+\s+(\S+)\s+(.*)")
					matchname, lociname, length, evaluehh, start, end, pident, description = re.match(pat, line).groups()
				except AttributeError:
					continue
				else:
					length = float(length)
					pident = 100.00*float(pident)
					protarea = 100.00*(((float(end)-1.00) - (float(start)-1.00))/length)
					try:
						data2 = infoprot_hmmer['lociname']
					except KeyError:
						if float(protarea) >= float(args.hmmercovthreshold) and float(evaluehh) <= float(args.hmmerevalue) and float(pident) >= float(args.hmmeridthreshold):
							infoprot_hmmer['lociname'] = matchname
							infoprot_hmmer['name'] = lociname
							infoprot_hmmer['evalue'] = float(evaluehh)
							infoprot_hmmer['pcover'] = float(protarea)
							infoprot_hmmer['pident'] = float(pident)
							infoprot_hmmer['descr'] = description
						else:
							try:
								if float(protarea) >= float(args.hmmercovthreshold) and float(evaluehh) <= float(args.hmmerevalue) and float(pident) >= float(args.hmmeridthreshold) and float(pident) >= infoprot_hmmer['pident']:
									infoprot_hmmer['lociname'] = matchname
									infoprot_hmmer['name'] = lociname
									infoprot_hmmer['evalue'] = float(evaluehh)
									infoprot_hmmer['pcover'] = float(protarea)
									infoprot_hmmer['pident'] = float(pident)
									infoprot_hmmer['descr'] = description
							except KeyError:
									continue
					else:
						if float(protarea) >= float(args.hmmercovthreshold) and float(evaluehh) <= float(args.hmmerevalue) and float(pident) >= float(args.hmmeridthreshold) and float(pident) >= infoprot_hmmer['pident']:
							infoprot_hmmer['lociname'] = matchname
							infoprot_hmmer['name'] = lociname
							infoprot_hmmer['evalue'] = float(evaluehh)
							infoprot_hmmer['pcover'] = float(protarea)
							infoprot_hmmer['pident'] = float(pident)
							infoprot_hmmer['descr'] = description
					information_proteins_hmmer[matchname] = infoprot_hmmer
	for contig in sorted(protsdict):
		for protein in sorted(protsdict[contig]):
			try:
				data3 = information_proteins_hmmer[protein]['name']
			except KeyError:
				protsdict[contig][protein]['pvog'] = "NA"
				protsdict[contig][protein]['refinement'] = "NO_REFINEMENT"
				protsdict[contig][protein]['pvog_pident'] = "NA"
				protsdict[contig][protein]['pvog_pcover'] = "NA"
				protsdict[contig][protein]['pvog_evalue'] = "NA"
			else:
				protsdict[contig][protein]['pvog'] = information_proteins_hmmer[protein]['name']
				protsdict[contig][protein]['refinement'] = "YES_PVOGS"
				protsdict[contig][protein]['pvog_pident'] = information_proteins_hmmer[protein]['pident']
				protsdict[contig][protein]['pvog_pcover'] = information_proteins_hmmer[protein]['pcover']
				protsdict[contig][protein]['pvog_evalue'] = information_proteins_hmmer[protein]['evalue']

	endtime10 = time.time()
	durationtime10 = endtime10 - starttime10
	eprint("Done: function prediction based on HMMER took %s seconds" % str(durationtime10))

## Creating the CSV table with all protein statistics
eprint("\nCreating all output files")
with open("%s.csv" % root_output, "w") as tablefile:
	print("\t".join(["Contig", "Protein ID", "Start", "Stop", "Strand", "size_aa", "pI", "Mol_weight_kDa", "Instability_index", "Description", "Source", "Perc_ID", "Perc_Cov", "E-value", "Refinement", "PVOG", "Perc_ID", "Perc_Cov", "E-value"]), file=tablefile)
	for contig in sorted(protsdict):
		for locus in sorted(protsdict[contig], key = stringSplitByNumbers):
			try:
				if protsdict[contig][locus]['refinement'] == "YES_PVOGS":
					print("\t".join([contig, locus, str(protsdict[contig][locus]['begin']), str(protsdict[contig][locus]['end']), str(protsdict[contig][locus]['strand']), str(protsdict[contig][locus]['length']), str(protsdict[contig][locus]['isoelectricpoint']), str(protsdict[contig][locus]['molweightkda']), str(protsdict[contig][locus]['instability']), protsdict[contig][locus]['descr'], protsdict[contig][locus]['source'], str(protsdict[contig][locus]['pident']), str(protsdict[contig][locus]['pcover']), str(protsdict[contig][locus]['evalue']), "YES", str(protsdict[contig][locus]['pvog']), str(protsdict[contig][locus]['pvog_pident']), str(protsdict[contig][locus]['pvog_pcover']), str(protsdict[contig][locus]['pvog_evalue'])]), file=tablefile)
				else:
					print("\t".join([contig, locus, str(protsdict[contig][locus]['begin']), str(protsdict[contig][locus]['end']), str(protsdict[contig][locus]['strand']), str(protsdict[contig][locus]['length']), str(protsdict[contig][locus]['isoelectricpoint']), str(protsdict[contig][locus]['molweightkda']), str(protsdict[contig][locus]['instability']), protsdict[contig][locus]['descr'], protsdict[contig][locus]['source'], str(protsdict[contig][locus]['pident']), str(protsdict[contig][locus]['pcover']), str(protsdict[contig][locus]['evalue']), "NO", "NO", "NA", "NA", "NA"]), file=tablefile)
			except KeyError:
				print("\t".join([contig, locus, str(protsdict[contig][locus]['begin']), str(protsdict[contig][locus]['end']), str(protsdict[contig][locus]['strand']), str(protsdict[contig][locus]['length']), str(protsdict[contig][locus]['isoelectricpoint']), str(protsdict[contig][locus]['molweightkda']), str(protsdict[contig][locus]['instability']), protsdict[contig][locus]['descr'], protsdict[contig][locus]['source'], str(protsdict[contig][locus]['pident']), str(protsdict[contig][locus]['pcover']), str(protsdict[contig][locus]['evalue']), "NO", "NO", "NA", "NA", "NA"]), file=tablefile)				

## Creating a new Genbank file
for newfile in sorted(glob.glob("CONTIG_*.fna")):
	newtempgbk = "%s.temp.gbk" % newfile
	with open(newfile, "rU") as basefile, open(newtempgbk, "w"):
		for record in SeqIO.parse(basefile, "fasta", IUPAC.ambiguous_dna):
			whole_sequence = SeqRecord(record.seq)
			whole_sequence.id = str(record.id)
			whole_sequence.annotations['data_file_division'] = args.typedata.upper()
			whole_sequence.annotations['date'] = strftime("%d-%b-%Y").upper()
			for locus_protein in sorted(protsdict):
				if locus_protein == record.id:
					for protein in sorted(protsdict[locus_protein], key = stringSplitByNumbers):
						putative_start = int(protsdict[locus_protein][protein]['begin'])
						start_pos = SeqFeature.ExactPosition(putative_start)
						end_pos = SeqFeature.ExactPosition(protsdict[locus_protein][protein]['end'])
						feature_location = SeqFeature.FeatureLocation(start_pos, end_pos, strand=protsdict[locus_protein][protein]['strand'])
						qualifiersgene = OrderedDict([('locus_tag', protsdict[locus_protein][protein]['protein_id'])])
						new_data_gene = SeqFeature.SeqFeature(feature_location, type = "gene", strand = protsdict[locus_protein][protein]['strand'], qualifiers = qualifiersgene)
						whole_sequence.features.append(new_data_gene)
						if args.nohmmer == False:
							if protsdict[locus_protein][protein]['pvog'] == "NA":
								qualifiers = [('locus_tag', protsdict[locus_protein][protein]['protein_id']), ('product', protsdict[locus_protein][protein]['descr']), ('protein_id', protsdict[locus_protein][protein]['protein_id']), ('translation', protsdict[locus_protein][protein]['translation'])]
							else:
								qualifiers = [('locus_tag', protsdict[locus_protein][protein]['protein_id']), ('product', protsdict[locus_protein][protein]['descr']), ('protein_id', protsdict[locus_protein][protein]['protein_id']), ('note', 'PVOGs: %s' % protsdict[locus_protein][protein]['pvog']), ('translation', protsdict[locus_protein][protein]['translation'])]
						else:
							qualifiers = [('locus_tag', protsdict[locus_protein][protein]['protein_id']), ('product', protsdict[locus_protein][protein]['descr']), ('protein_id', protsdict[locus_protein][protein]['protein_id']), ('translation', protsdict[locus_protein][protein]['translation'])]
						feature_qualifiers = OrderedDict(qualifiers)
						new_data_cds = SeqFeature.SeqFeature(feature_location, type = "CDS", strand = protsdict[locus_protein][protein]['strand'], qualifiers = feature_qualifiers)
						whole_sequence.features.append(new_data_cds)
			for locus_rRNA in sorted(subunits):
				if locus_rRNA == record.id:
					for rRNA in subunits[locus_rRNA].keys():
						count = 0
						try:
							lengthlist = len(subunits[locus_rRNA][rRNA])
						except KeyError:
							continue
						else:
							while count < lengthlist:
								putative_start = int(subunits[locus_rRNA][rRNA][count]['begin'])
								start_pos = SeqFeature.ExactPosition(putative_start)
								end_pos = SeqFeature.ExactPosition(subunits[locus_rRNA][rRNA][count]['end'])
								feature_location = SeqFeature.FeatureLocation(start_pos, end_pos, strand=subunits[locus_rRNA][rRNA][count]['strand'])
								new_data_gene = SeqFeature.SeqFeature(feature_location, type = "gene", strand = subunits[locus_rRNA][rRNA][count]['strand'])
								whole_sequence.features.append(new_data_gene)
								qualifiers = [('product', subunits[locus_rRNA][rRNA][count]['product'])]
								feature_qualifiers = OrderedDict(qualifiers)
								new_data_rRNA = SeqFeature.SeqFeature(feature_location, type = "rRNA", strand = subunits[locus_rRNA][rRNA][count]['strand'], qualifiers = feature_qualifiers)
								whole_sequence.features.append(new_data_rRNA)
								count += 1
			for locus_tRNA in sorted(tRNAdict):
				if locus_tRNA == record.id:
					for tRNA in sorted(tRNAdict[locus_tRNA], key = stringSplitByNumbers):
						putative_start = int(tRNAdict[locus_tRNA][tRNA]['begin'])
						start_pos = SeqFeature.ExactPosition(putative_start)
						end_pos = SeqFeature.ExactPosition(tRNAdict[locus_tRNA][tRNA]['end'])
						feature_location = SeqFeature.FeatureLocation(start_pos, end_pos, strand=tRNAdict[locus_tRNA][tRNA]['strand'])
						new_data_gene = SeqFeature.SeqFeature(feature_location, type = "gene", strand = tRNAdict[locus_tRNA][tRNA]['strand'])
						whole_sequence.features.append(new_data_gene)
						qualifiers = [('product', tRNAdict[locus_tRNA][tRNA]['product'])]
						feature_qualifiers = OrderedDict(qualifiers)
						new_data_tRNA = SeqFeature.SeqFeature(feature_location, type = "tRNA", strand = tRNAdict[locus_tRNA][tRNA]['strand'], qualifiers = feature_qualifiers)
						whole_sequence.features.append(new_data_tRNA)
			for locus_tmRNA in sorted(tmRNAdict):
				if locus_tmRNA == record.id:
					for tmRNA in sorted(tmRNAdict[locus_tmRNA], key = stringSplitByNumbers):
						putative_start = int(tmRNAdict[locus_tmRNA][tmRNA]['begin'])
						start_pos = SeqFeature.ExactPosition(putative_start)
						end_pos = SeqFeature.ExactPosition(tmRNAdict[locus_tmRNA][tmRNA]['end'])
						feature_location = SeqFeature.FeatureLocation(start_pos, end_pos, strand=tmRNAdict[locus_tmRNA][tmRNA]['strand'])
						new_data_gene = SeqFeature.SeqFeature(feature_location, type = "gene", strand = tmRNAdict[locus_tmRNA][tmRNA]['strand'])
						whole_sequence.features.append(new_data_gene)
						qualifiers = [('product', tmRNAdict[locus_tmRNA][tmRNA]['product'])]
						feature_qualifiers = OrderedDict(qualifiers)
						new_data_tmRNA = SeqFeature.SeqFeature(feature_location, type = "tmRNA", strand = tmRNAdict[locus_tmRNA][tmRNA]['strand'], qualifiers = feature_qualifiers)
						whole_sequence.features.append(new_data_tmRNA)
			for locus_CRISPR in sorted(information_CRISPR, key = stringSplitByNumbers):
				if locus_CRISPR == record.id:
					for CRISPR in sorted(information_CRISPR[locus_CRISPR], key = stringSplitByNumbers):
						putative_start = int(information_CRISPR[locus_CRISPR][CRISPR]['start'])
						start_pos = SeqFeature.ExactPosition(putative_start)
						end_pos = SeqFeature.ExactPosition(information_CRISPR[locus_CRISPR][CRISPR]['end'])
						feature_location = SeqFeature.FeatureLocation(start_pos, end_pos)
						qualifiers = [('rpt_family', 'CRISPR'), ('rpt_type', 'direct'), ('rpt_unit_range', "%i..%i" % (int(information_CRISPR[locus_CRISPR][CRISPR]['start']), int(information_CRISPR[locus_CRISPR][CRISPR]['repeatend']))), ('rpt_unit_seq', information_CRISPR[locus_CRISPR][CRISPR]['repeatseq'])]
						feature_qualifiers = OrderedDict(qualifiers)
						new_data_CRISPRrepeat = SeqFeature.SeqFeature(feature_location, type = "repeat_region", qualifiers = feature_qualifiers)
						whole_sequence.features.append(new_data_CRISPRrepeat)
			for locus_tandem in sorted(information_TRF):
				if locus_tandem == record.id:
					for tandem in sorted(information_TRF[locus_tandem]):
						putative_start = int(information_TRF[locus_tandem][tandem]['start'])
						start_pos = SeqFeature.ExactPosition(putative_start)
						end_pos = SeqFeature.ExactPosition(information_TRF[locus_tandem][tandem]['end'])
						feature_location = SeqFeature.FeatureLocation(start_pos, end_pos)
						qualifiers = [('rpt_type', 'direct')]
						feature_qualifiers = OrderedDict(qualifiers)
						new_data_tandemrepeat = SeqFeature.SeqFeature(feature_location, type = "repeat_region", qualifiers = feature_qualifiers)
						whole_sequence.features.append(new_data_tandemrepeat)
			for locus_inverted in sorted(information_IRF):
				if locus_inverted == record.id:
					for inverted in sorted(information_IRF[locus_inverted]):
						putative_start = int(information_IRF[locus_inverted][inverted]['start'])
						start_pos = SeqFeature.ExactPosition(putative_start)
						end_pos = SeqFeature.ExactPosition(information_IRF[locus_inverted][inverted]['end'])
						feature_location = SeqFeature.FeatureLocation(start_pos, end_pos)
						qualifiers = [('rpt_type', 'inverted')]
						feature_qualifiers = OrderedDict(qualifiers)
						new_data_invertedrepeat = SeqFeature.SeqFeature(feature_location, type = "repeat_region", qualifiers = feature_qualifiers)
						whole_sequence.features.append(new_data_invertedrepeat)
			SeqIO.write(whole_sequence, newtempgbk, "genbank")
	newgbk = "%s.gbk" % newfile
	with open(newtempgbk, "rU") as gbktempfile, open(newgbk, "w") as gbkrealfile:
		newpat = re.compile("D|RNA\s+(CON|PHG|VRL|BCT)")
		newpat2 = re.compile("LOCUS\s+(\S+)")
		for line in gbktempfile:
			if line.startswith("LOCUS ") and re.search(newpat, line):
				namelocus = re.match(newpat2, line).groups()[0]
				if genomeshape[namelocus]['genomeshape'] == "linear":
					newline = re.sub("bp    DNA\s+", "bp    DNA     linear   ", line)
				else:
					newline = re.sub("bp    DNA\s+", "bp    DNA     circular ", line)
				gbkrealfile.write(newline)
			else:
				gbkrealfile.write(line)
	for f in glob.glob("*.temp.gbk"):
		os.remove(f)
cat_all(sorted(glob.glob('CONTIG_*.fna.gbk')), "%s.gbk" % root_output)

## Printing the GFF output (optional)
if args.gffprint==True:
	with open(("%s.gff" % root_output), "w") as outgff, open("%s.gbk" % root_output, "rU") as ingbk:
		BCBio.GFF.write(SeqIO.parse(ingbk, "genbank"), outgff)

## Removing all intermediate files
eprint("Cleaning all intermediate files")
os.remove("CONTIGS_ALL.fasta")
os.remove("temporal_circular.fasta")
os.remove("rrnafile.csv")
for splittedtrnafiles in sorted(glob.glob('trnafile_*.fasta')):
	os.remove(splittedtrnafiles)
os.remove("crisprfile.txt")
os.remove("trf_temp.dat")
os.remove("irf_temp.dat")
os.remove("pretemp.faa")
for splittedorffiles in sorted(glob.glob('orffile_*.faa')):
	os.remove(splittedorffiles)
os.remove("PROTS_FIRST_ROUND.faa")
os.remove("PROTS_FIRST_ROUND.csv")
if args.nohmmer == False:
	os.remove("PROTS_FIRST_ROUND.tbl")
os.remove("PROTS_SECOND_ROUND.faa")
os.remove("PROTS_SECOND_ROUND.csv")
for splittedfastafiles in sorted(glob.glob('CONTIG_*.fna')):
	os.remove(splittedfastafiles)
for tempgbk in glob.glob("CONTIG_*.gbk"):
	os.remove(tempgbk)
 
## Preparing sequences for GenBank submission (Original code from Wan Yu's gbk2tbl.py script [https://github.com/wanyuac/BINF_toolkit/blob/master/gbk2tbl.py])
allowed_qualifiers = ['locus_tag', 'gene', 'product', 'pseudo', 'protein_id', 'gene_desc', 'old_locus_tag', 'note', 'inference', 'organism', 'mol_type', 'strain', 'sub_species', 'isolation-source', 'country']
with open(args.modifiers, "rU") as modifiers, open("%s.gbk" % root_output, "r") as genbank_fh, open("%s.fasta" % root_output, "w") as fasta_fh, open("%s.tbl" % root_output, "w") as feature_fh: 
	info = modifiers.readline()
	wholelist = list(SeqIO.parse(genbank_fh, 'genbank'))
	for record in wholelist:
		if len(record) <= args.mincontigsize:
			eprint("WARNING: Skipping small contig %s" % record.id)
			continue
		record.description = "%s %s" % (record.id, info)
		SeqIO.write([record], fasta_fh, 'fasta')
		print('>Feature %s' % (record.name), file=feature_fh)
		for line in record.features:
			if line.strand == 1:
				print('%d\t%d\t%s' % (line.location.nofuzzy_start + 1, line.location.nofuzzy_end, line.type), file=feature_fh)
			else:
				print('%d\t%d\t%s' % (line.location.nofuzzy_end, line.location.nofuzzy_start + 1, line.type), file=feature_fh)
			for (key, values) in line.qualifiers.iteritems():
				if key not in allowed_qualifiers:
					continue
				for v in values:
					print('\t\t\t%s\t%s' % (key, v), file=feature_fh)

## Final statement
eprint("\nGenome annotation done!")
eprint("The GenBank file is %s.gbk" % root_output)
if args.gffprint==True:
	eprint("The GFF3 file is %s.gff" % root_output)
eprint("The table file for GenBank submission is %s.tbl" % root_output)
eprint("The FASTA file for GenBank submission is %s.fasta" % root_output)
eprint("The table file with all protein information is %s.csv" % root_output)
