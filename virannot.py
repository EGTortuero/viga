# -*- coding: utf-8 -*-

#!/usr/bin/env python

# Virannot - De-novo viral genome annotator
#
# Copyright (C) 2017 - Enrique Gonzalez-Tortuero
#                      Vimalkumar Velayudhan
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

# Importing python libraries
import argparse
import csv
import fileinput
import fractions
import glob
import os
import re
import sys
import subprocess
from BCBio import GFF
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import OrderedDict
from time import strftime

# Preparing functions
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

def cmd_exists(cmd):
	return subprocess.call("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

def stringSplitByNumbers(x):
	r = re.compile('(\d+)')
	l = r.split(x)
	return [int(y) if y.isdigit() else y for y in l]

# Defining the program version
version = "0.2.0"

# Processing the parameters
parser = argparse.ArgumentParser(description='Virannot is a automatic de novo viral genome annotator.')
basic_group = parser.add_argument_group('Basic options for virannot [REQUIRED]')

basic_group.add_argument("--input", dest="inputfile", type=str, required=True, help='Input file as a FASTA file', metavar="FASTAFILE")
basic_group.add_argument("--blastdb", dest="blastdatabase", type=str, required=True, help='BLAST Database that will be used for the protein function prediction. The database must be an amino acid one, not  nucleotidic', metavar="BLASTDB")
basic_group.add_argument("--hhblitsdb", dest="hhblitsdatabase", type=str, required=True, help='HHBLITS Database that will be used for the first step of the protein function prediction.In this case, HHBLITSDB should be in the format "/full/path/to/db1/db1 (without the extension _a3m_db)"', metavar="HHBLITSDB")
basic_group.add_argument("--hhsearchdb", dest="hhsearchdatabase", nargs='+', type=str, required=True, help='HHSEARCH Database(s) that will be used for the second step of the protein function prediction. You can use more than a single database for the analysis. In this case, HHSEARCHDB should be in the format -d3 /full/path/to/db1/db1_hhm_db /full/path/to/db2/db2_hhm_db', metavar="HHSEARCHDB")
basic_group.add_argument("--modifiers", dest="modifiers", type=str, required=True, help='Input file as a plain text file with the modifiers per every FASTA header according to SeqIn (https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html). All modifiers must be written in a single line and are separated by a single space character. No space should be placed besides the = sign. For example: [organism=Serratia marcescens subsp. marcescens] [sub-species=marcescens] [strain=AH0650_Sm1] [topology=linear] [moltype=DNA] [tech=wgs] [gcode=11] [country=Australia] [isolation-source=sputum]. This line will be copied and printed along with the record name as the definition line of every contig sequence.', metavar="TEXTFILE")

advanced_group = parser.add_argument_group('Advanced options for virannot [OPTIONAL]')
advanced_group.add_argument("--readlength", dest="read_length", type=int, default=101, help='Read length for the circularity prediction (default: 101 bp)', metavar="INT")
advanced_group.add_argument("--out", dest="rootoutput", type=str, help='Name of the outputs files (without extension)', metavar="OUTPUTNAME")
advanced_group.add_argument("--locus", dest="locus", type=str, default='LOC', help='Name of the sequences. If the input is a multiFASTA file, please put a general name as the program will add the number of the contig at the end of the name (Default: %(default)s)', metavar="STRING")
advanced_group.add_argument("--threads", dest="ncpus", default=1, help='Number of threads/cpus (Default: %(default)s cpu)', metavar="INT")
advanced_group.add_argument("--noparallel", dest="noparallel", action='store_true', default=False, help="Don't use of GNU Parallel to run BLAST and HHSEARCH in parallel jobs (Default=FALSE)")
advanced_group.add_argument("--gff", dest="gffprint", action='store_true', default=False, help='Printing the output as GFF3 file (Default: False)')
advanced_group.add_argument("--blastevalue", dest="blastevalue", default=0.00001, help='Blast e-value threshold (Default: 0.00001)', metavar="FLOAT")
advanced_group.add_argument("--hhsuiteevalue", dest="hhsuiteevalue", default=0.001, help='HHSUITE e-value threshold (Default: 0.001)', metavar="FLOAT")

type_choices = {'CON': 'Contig', 'VRL': 'Eukaryotic/Archaea virus', 'PHG': 'Phages'}
type_help = ('GenBank Division: One of the following codes - {0}. (Default: %(default)s)'.format(', '.join('{0} ({1})'.format(k, v) for k, v in type_choices.items())))
advanced_group.add_argument("--typedata", dest="typedata", type=str, default='CON', help=type_help, metavar="CON|VRL|PHG")

gcode_choices = {'1': 'Standard genetic code [Eukaryotic]', '4': 'Mycoplasma/Spiroplasma', '6': 'Protozoa [nuclear]', '11': 'Bacteria/Archaea/Phages', '25': 'Gracilibacteria & Candidate division SR1'}
gcode_help = ('Number of GenBank translation table. At this moment, the available options are {0}. (Default: %(default)s)'.format('{}'.format(', '.join('{0} ({1})'.format(k, v) for k, v in gcode_choices.items()))))
advanced_group.add_argument("--gcode", dest="gcode", type=str, default='11', help=gcode_help, metavar="NUMBER")

advanced_group.add_argument('--mincontigsize', dest="mincontigsize", type=int, default = 200, help = 'Minimum contig length to be considered in the final files (Default: 200 bp)', metavar="INT")
advanced_group.add_argument("--idthr", dest="idthreshold", default=50.00, help='ID threshold (Default: 50.0)', metavar="FLOAT")
advanced_group.add_argument("--coverthr", dest="covthreshold", default=50.00, help='Coverage threshold (Default: 50.0)', metavar="FLOAT")
advanced_group.add_argument("--maxfilt", dest="maxfilt", default=20000, help='Max number of hits allowed to pass 2nd prefilter in HHSUITE (Default: 20000)', metavar="INT (>100)")
advanced_group.add_argument("--neffmax", dest="neffmax", default=10, help="Skip further search iterations in HHSUITE when diversity of query MSA becomes larger than this value (Default: 10)", metavar="FLOAT [1.00-20.00]")
advanced_group.add_argument("--diffid", dest="diffid", default=5.00, help='Max allowed difference between the ID percentages of BLAST and HHSEARCH. (Default = 5.00; Not recommended to change such value)', metavar="FLOAT (>0.01)")
advanced_group.add_argument("--blastexh", dest="blastexh", action='store_true', default=False, help='Use of exhaustive BLAST to predict the proteins by homology according to Fozo et al. (2010) Nucleic Acids Res (Default=FALSE)')
advanced_group.add_argument("--hhsearchexh", dest="hhsearchexh", action='store_true', default=False, help='Use of exhaustive HHSEARCH to tackle proteins of unknown function according to Fidler et al. (2016) Traffic (Default=FALSE)')
args = parser.parse_args()

root_output = args.rootoutput
if not root_output:
    root_output = '{}_annotated'.format(os.path.splitext(args.inputfile)[0])

hh_search_dbs = '{!r}'.format(' '.join(args.hhsearchdatabase))

# Printing the header of the program 
print "This is VirAnnot ", version
print "Written by Enrique Gonzalez Tortuero & Vimalkumar Velayudhan"
print "Homepage is https://github.com/EGTortuero/virannot"
print "Local time: ", strftime("%a, %d %b %Y %H:%M:%S")
print "\n\n"

# checking the presence of the programs in the system
if not cmd_exists("lastz")==True:
	sys.exit("You must install LASTZ before using this script")
elif not cmd_exists("prodigal")==True:
	sys.exit("You must install prodigal before using this script")
elif not cmd_exists("parallel")==True and args.noparallel==False:
	sys.exit("You must install GNU Parallel before using this script")
elif not cmd_exists("blastp")==True:
	sys.exit("You must install BLAST before using this script")
elif not cmd_exists("hhblits")==True and not cmd_exists("hhsearch")==True:
	sys.exit("You must install HHsuite before using this script")
elif not cmd_exists("aragorn")==True:
	sys.exit("You must install ARAGORN before using this script")
elif not cmd_exists("trf")==True:
	sys.exit("You must install Tandem Repeats Finder before using this script")
elif not cmd_exists("irf")==True:
	sys.exit("You must install Inverted Repeats Finder before using this script")
print "\nData type is {0} and GenBank translation table no is {1}\n".format(args.typedata, args.gcode)

# Correcting the original file (long headers problem + multiple FASTA files)
record_iter = SeqIO.parse(open(args.inputfile, "rU"),"fasta")
counter = 1
for i, batch in enumerate(batch_iterator(record_iter, 1)):
	seq_index = "CONTIG_%i" % (i + 1)
	filename = "%s.temp.fna" % seq_index
	newfilename = "%s.fna" % seq_index
	with open(filename, "w") as handle:
		count = SeqIO.write(batch, filename, "fasta")

	with open(filename, "rU") as original, open(newfilename, "w") as corrected:
		sequences = SeqIO.parse(original, "fasta", IUPAC.ambiguous_dna)
		for record in sequences:
			original_name = record.id
			record.id = "%s_%i" % (args.locus, counter)
			record.description = record.description
			counter += 1
			print "WARNING: The name of the sequence %s was corrected as %s" % (original_name, record.id)
		SeqIO.write(record, corrected, "fasta")
	os.remove(filename)

for newfile in sorted(glob.glob("CONTIG_*.fna")):

	# Predicting the shape of the contig (code based on Alex Crits-Christoph's find_circular.py script [https://github.com/alexcritschristoph/VICA/blob/master/find_circular.py])
	print "Predicting the shape of the contig using LASTZ"
	genomeshape = {}
	with open(newfile, "rU") as targetfasta:
		Sequence = SeqIO.parse(newfile, "fasta")
		for record in Sequence:
			seq_beginning = str(record.seq[0:args.read_length])
			seq_ending = str(record.seq[len(record.seq)-args.read_length:len(record.seq)])
			combined_seqs = SeqRecord(Seq(seq_beginning + seq_ending, IUPAC.ambiguous_dna), id = "test")
			SeqIO.write(combined_seqs, "temporal_circular.fasta", "fasta")
			outputlastz = subprocess.check_output(["lastz", "temporal_circular.fasta", "--self", "--notrivial", "--nomirror", "--format=general-:start1,end1,start2,end2,score,strand1,strand2,identity,length1"]);
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
						genomeshape['genomeshape'] = "circular"
						try:
							if genomeshape['identity'] >= float(fractions.Fraction(identity)):
								genomeshape['identity'] = float(fractions.Fraction(identity))
								genomeshape['length'] = length
						except KeyError:
							genomeshape['identity'] = float(fractions.Fraction(identity))
							genomeshape['length'] = length
					else:
						continue
			try:
				if genomeshape['genomeshape'] == "":
						genomeshape['genomeshape'] = "linear"
			except KeyError:
				genomeshape['genomeshape'] = "linear"
			else:
				genomeshape['genomeshape'] = "circular"
				with open("temp.fasta", "w") as correctedcircular:
					Corrseq = str(record.seq[int(genomeshape['length'])//2:-int(genomeshape['length'])//2])
					Newseq = SeqRecord(Seq(Corrseq, IUPAC.ambiguous_dna), id = record.id)
					SeqIO.write(Newseq, correctedcircular, "fasta")
				os.rename("temp.fasta", newfile)
		print "Done. LASTZ predicted that the contig is %s\n" % genomeshape['genomeshape']
	
	# Predicting genes using PRODIGAL
	print "Running Prodigal to predict the genes in %s" % newfile
	for record in SeqIO.parse(newfile, "fasta"):
		length_contig = len(record.seq)
		if (length_contig >= 100000):
			if genomeshape['genomeshape'] == 'linear':
				subprocess.call(["prodigal", "-a", "pretemp.faa", "-i", newfile, "-o", "/dev/null", "-g", args.gcode, "-c", "-q"])
			else:
				subprocess.call(["prodigal", "-a", "pretemp.faa", "-i", newfile, "-o", "/dev/null", "-g", args.gcode, "-q"])
		else:
			if genomeshape['genomeshape'] == 'linear':
				subprocess.call(["prodigal", "-a", "pretemp.faa", "-i", newfile, "-o", "/dev/null", "-p", "meta", "-g", args.gcode, "-c", "-q"])
			else:
				subprocess.call(["prodigal", "-a", "pretemp.faa", "-i", newfile, "-o", "/dev/null", "-p", "meta", "-g", args.gcode, "-q"])
	num_seqs = len(list(SeqIO.parse("pretemp.faa", "fasta")))
	print "Done. Prodigal was able to predict %i genes\n" % num_seqs
	
	with open("pretemp.faa", "rU") as originalfaa, open("temp.faa", "w") as correctedfaa:
		sequences = SeqIO.parse(originalfaa, "fasta")
		for record in sequences:
			record.seq = record.seq.rstrip("*")
			SeqIO.write(record, correctedfaa, "fasta")
	os.remove("pretemp.faa")

	# Predicting the function of the proteins based on homology using BLAST
	equivalence = {}
	record_iter = SeqIO.parse(open("temp.faa"),"fasta")
	for i, batch in enumerate(batch_iterator(record_iter, 1)):
		seq_index = "SEQ_%i" % (i + 1)
		filename = "%s.faa" % seq_index
		with open(filename, "w") as handle:
			count = SeqIO.write(batch, handle, "fasta")
			equivalence[seq_index] = batch[0].id

	if args.noparallel==True:
		if args.blastexh==True:
			print "Running BLAST to predict the genes according to homology inference in %s using exhaustive mode (see Fozo et al. (2010) Nucleic Acids Res for details)" % newfile
			subprocess.call(['blastp', '-query', "temp.faa", '-db', args.blastdatabase, '-evalue', str(args.blastevalue), '-outfmt', '"6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle', '-out', 'temp_blast.csv', '-max_target_seqs', '10', '-word_size', '2', '-gapopen', '8', '-gapextend', '2', '-matrix', '"PAM70"', '-comp_based_stats', '"0"', "-num_threads", str(args.ncpus)])
		else:
			print "Running BLAST to predict the genes according to homology inference in %s using default parameters" % newfile
			subprocess.call(['blastp', '-query', "temp.faa", '-db', args.blastdatabase, '-evalue', str(args.blastevalue), '-outfmt', '"6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle"', '-out', 'temp_blast.csv', '-max_target_seqs', '10', "-num_threads", str(args.ncpus)])
	else:
		with open("commands.sh", "w") as commands:
			for j in sorted(glob.glob("SEQ_*.faa")):
				jout = "%s.blast.csv" % j
				if args.blastexh==True:
					print "Creating file to run parallel BLAST: adding %s using exhaustive mode  (see Fozo et al. (2010) Nucleic Acids Res for details)" % j
					line = ['blastp', '-query', j, '-db', args.blastdatabase, '-evalue', str(args.blastevalue), '-outfmt', '"6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle"', '-out', jout, '-max_target_seqs', '10', '-word_size', '2', '-gapopen', '8', '-gapextend', '2', '-matrix', '"PAM70"', '-comp_based_stats', '"0"', '\n']
					line2write = ' '.join(line)
					commands.write(line2write)
				else:
					print "Creating file to run parallel BLAST: adding %s using default parameters" % j
					line = ['blastp', '-query', j, '-db', args.blastdatabase, '-evalue', str(args.blastevalue), '-outfmt', '"6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle"', '-out', jout, '-max_target_seqs', '10', '\n']
					line2write = ' '.join(line)
					commands.write(line2write)
		print "Running parallel BLAST"
		subprocess.call(['parallel', '-j', str(args.ncpus)], stdin=open('commands.sh'))
		listblastcsv = sorted(glob.glob('SEQ_*.csv'))
		with open('temp_blast.csv', 'w') as finalblastcsv:
			for fname in listblastcsv:
				with open(fname) as infile:
					for line in infile:
						finalblastcsv.write(line)
		os.remove("commands.sh")
		for fname in listblastcsv:
			os.remove(fname)
	print "Done. BLASTp was done to predict the genes by homology\n"

	# Parsing the results from BLAST
	with open("temp_blast.csv", "rU") as blastresults:
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
				if (float(perc_id) >= float(args.idthreshold)) and (float(perc_cover) >= float(args.covthreshold)) and (float(row['evalue']) <= float(args.blastevalue)):
					information_proteins_blast[row['qseqid']] = infoprot_blast
			else:
				if (float(perc_cover) > float(data['pcover'])) and (float(perc_cover) >= float(args.covthreshold)) and (float(row['evalue']) <= float(args.blastevalue) and (float(row['evalue']) <= float(data['evalue']))):
					information_proteins_blast[row['qseqid']] = infoprot_blast
	
	## Predicting the function of the proteins based on HMM-HMM comparisons using HH-SUITE

	with open("commandsB.sh", "w") as commandsB, open("commandsC.sh", "w") as commandsC:
		for singleprot in sorted(glob.glob("SEQ_*.faa")):
			hhtempout = "%s.a3m" % singleprot
			hhout = "%s.hhr" % singleprot
			if args.noparallel==True:
				if args.hhsearchexh==True:
					print "Running HHblits to predict the proteins according to HMM-HMM comparisons in %s using exhaustive mode (see Fidler et al. (2016) Traffic for details)." % singleprot
					subprocess.call(['hhblits', '-i', singleprot, '-d', args.hhblitsdatabase, '-o' , '/dev/null', '-oa3m', hhtempout, '-n', '8', '-e', '0.001', '-E', '0.01', '-maxfilt', str(args.maxfilt), '-neffmax', str(args.neffmax), '-v', "0", "-cpu", str(args.ncpus)])
					print "Done. HHblits was done to predict the function of the genes by HMM-HMM comparisons\n\n"
					print "Running HHpred to predict the proteins according to HMM-HMM comparisons in %s using exhaustive mode (see Fidler et al. (2016) Traffic for details)." % singleprot
					subprocess.call(['hhsearch', '-i', hhtempout, '-d', hh_search_dbs, '-o', hhout, '-mac', '-e', '0.01', '-v', "0", "-cpu", str(args.ncpus)])
				else:
					print "Running HHblits to predict the proteins according to HMM-HMM comparisons in %s using default parameters." % singleprot
					subprocess.call(['hhblits', '-i', singleprot, '-d', args.hhblitsdatabase, '-o' , '/dev/null', '-oa3m', hhtempout, '-n', '3', '-E', '0.5', '-maxfilt', str(args.maxfilt), '-neffmax', str(args.neffmax), '-v', "0", "-cpu", str(args.ncpus)])
					print "Done. HHblits was done to predict the function of the genes by HMM-HMM comparisons\n\n"
					print "Running HHpred to predict the proteins according to HMM-HMM comparisons in %s using default parameters." % singleprot
					subprocess.call(['hhsearch', '-i', hhtempout, '-d', hh_search_dbs, '-o', hhout, '-v', "0", "-cpu", str(args.ncpus)])
			else:
				if args.hhsearchexh==True:
					print "Creating file to run parallel HHblits: adding %s using exhaustive mode (see Fidler et al. (2016) Traffic for details)." % singleprot
					lineB = ['hhblits', '-i', singleprot, '-d', args.hhblitsdatabase, '-o' , '/dev/null', '-oa3m', hhtempout, '-n', '8', '-e', '0.001', '-E', '0.01', '-maxfilt', str(args.maxfilt), '-neffmax', str(args.neffmax), '-v', "0", '\n']
					line2writeB = ' '.join(lineB)
					commandsB.write(line2writeB)
					print "Creating file to run parallel HHpred: adding %s using exhaustive mode (see Fidler et al. (2016) Traffic for details)." % singleprot
					lineC = ['hhsearch', '-i', hhtempout, '-d', hh_search_dbs, '-o', hhout, '-mac', '-e', '0.01', '-v', "0", '\n']
					line2writeC = ' '.join(lineC)
					commandsC.write(line2writeC)
				else:
					print "Creating file to run parallel HHblits: adding %s using default parameters." % singleprot
					lineB = ['hhblits', '-i', singleprot, '-d', args.hhblitsdatabase, '-o' , '/dev/null', '-oa3m', hhtempout, '-n', '3', '-E', '0.5', '-maxfilt', str(args.maxfilt), '-neffmax', str(args.neffmax), '-v', "0", '\n']
					line2writeB = ' '.join(lineB)
					commandsB.write(line2writeB)
					print "Creating file to run parallel HHpred: adding %s using default parameters." % singleprot
					lineC = ['hhsearch', '-i', hhtempout, '-d', hh_search_dbs, '-o', hhout, '-v', "0", '\n']
					line2writeC = ' '.join(lineC)
					commandsC.write(line2writeC)

	if args.noparallel==False:
		print "Running parallel HHblits"
		subprocess.call(['parallel', '-j', str(args.ncpus)], stdin=open('commandsB.sh'))
		print "Running parallel HHpred"		
		subprocess.call(['parallel', '-j', str(args.ncpus)], stdin=open('commandsC.sh'))
		os.remove("commandsB.sh")
		os.remove("commandsC.sh")
		for j in sorted(glob.glob("*a3m")):
			os.remove(j)
	print "Done. HHpred was done to predict the function of the genes by HMM-HMM comparisons\n\n"

	# Parsing the results from HH-SUITE
	information_proteins_hhsuite = {}
	for singlehhr in sorted(glob.glob("*.hhr")):
		rootname = singlehhr.replace(".faa.hhr", "")
		with open(singlehhr) as hhrfile:
			infoprot_hhsuite = {}
			match_found = False
			for line in hhrfile:
				try:
					pat = re.compile('^(\s+)?(\d+)\s(.{30})\s+(\d+\.\d+)\s+((?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?)\s+(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?\s+\d+\.\d+\s+\d+\.\d+\s+\d+\s+(\d+-\d+)')
					elm2rem, key, name, prob, evaluehh, coords = re.match(pat, line).groups()
				except AttributeError:
					try:
						data = infoprot_hhsuite['key']
					except KeyError:
						Matchcolsentence = re.compile("Match_columns ")
						if re.match(Matchcolsentence, line):
							length = float(line.replace("Match_columns ",""))
						else:
							continue
					else:
						Infoprothhsuitepat = re.compile(r"No "+infoprot_hhsuite['key']+"\s+$")
						if re.match(Infoprothhsuitepat, line):
							match_found = True
							continue
						if match_found:
							if line.startswith('>'):
								pat2 = re.compile('^>\S+\s(.+)$')
								infoprot_hhsuite['descr'] = re.match(pat2, line).groups()[0]
							if line.startswith('Probab'):
								pat3 = re.compile('Identities=(\d+)\%')
								infoprot_hhsuite['pident'] = float(re.search(pat3, line).groups()[0])
								match_found = False
						else:
							pass
				else:
					if float(prob) >= 50.00:
						start, end = coords.split('-')
						protarea = 100.00*((float(end)-1.00) - (float(start)-1.00))/length
						try:
							data2 = infoprot_hhsuite['key']
						except KeyError:
							if float(protarea) >= float(args.covthreshold) and (float(evaluehh) <= float(args.hhsuiteevalue)):
								infoprot_hhsuite['key'] = key
								infoprot_hhsuite['prob'] = float(prob)
								infoprot_hhsuite['evalue'] = float(evaluehh)
								infoprot_hhsuite['pcover'] = float(protarea)
								infoprot_hhsuite['name'] = name
						else:
							if (float(protarea) >= float(args.covthreshold)) and (float(protarea) >= float(infoprot_hhsuite['pcover'])) and (float(prob) >= float(infoprot_hhsuite['prob'])) and (float(evaluehh) <= float(args.hhsuiteevalue)) and (float(evaluehh) <= float(infoprot_hhsuite['evalue'])):
								infoprot_hhsuite['key'] = key
								infoprot_hhsuite['prob'] = float(prob)
								infoprot_hhsuite['evalue'] = float(evaluehh)
								infoprot_hhsuite['pcover'] = float(protarea)
								infoprot_hhsuite['name'] = name
			information_proteins_hhsuite[rootname] = infoprot_hhsuite

	#Storing protein information in memory
	with open("temp.faa", "rU") as protsfile:
		tempprotsdict = {}
		for protseq in SeqIO.parse(protsfile, "fasta"):
			tempindprot = {}
			dataprot = protseq.description.split(' # ')
			modseq = str(protseq.seq).replace("X","") # Removing all ambiguous amino acids to avoid problems with ProteinAnalysis module
			analysed_seq = ProteinAnalysis(modseq)
			tempindprot['length'] = len(protseq.seq)
			tempindprot['isoelectricpoint'] = analysed_seq.isoelectric_point()
			tempindprot['molweightkda'] = analysed_seq.molecular_weight()/1000.00
			tempindprot['instability'] = analysed_seq.instability_index()
			tempindprot['protein_id'] = dataprot[0]
			tempindprot['strand'] = int(dataprot[3])
			tempindprot['begin'] = int(dataprot[1])-1
			tempindprot['end'] = int(dataprot[2])
			tempprotsdict[dataprot[0]] = tempindprot

	# Creation of table
	debugtable = "%s.csv" % newfile
	with open(debugtable, "w") as tablefile:
		print >>tablefile, "\t".join(["Identifier", "Start", "Stop", "Strand", "size_aa", "pI", "Mol_weight_kDa", "Instability_index", "ID_BLAST", "Descr_BLAST", "evalue_BLAST", "%ID_BLAST", "%Cover_BLAST", "ID_HHPRED", "Descr_HHPRED", "Prob_HHPRED", "evalue_HHPRED", "%ID_HHPRED", "%Cover_HHPRED"])
		keylist = information_proteins_hhsuite.keys()
		keylist.sort()
		for keyB in keylist:
			keyB = keyB.replace(".hhr", "")
			try:
				print >>tablefile, "\t".join([equivalence[keyB], str(tempprotsdict[equivalence[keyB]]['begin']), str(tempprotsdict[equivalence[keyB]]['end']), str(tempprotsdict[equivalence[keyB]]['strand']), str(tempprotsdict[equivalence[keyB]]['length']), str(tempprotsdict[equivalence[keyB]]['isoelectricpoint']), str(tempprotsdict[equivalence[keyB]]['molweightkda']), str(tempprotsdict[equivalence[keyB]]['instability']), information_proteins_blast[equivalence[keyB]]['sseqid'], information_proteins_blast[equivalence[keyB]]['descr'], str(information_proteins_blast[equivalence[keyB]]['evalue']), str(information_proteins_blast[equivalence[keyB]]['pident']), str(information_proteins_blast[equivalence[keyB]]['pcover']), information_proteins_hhsuite[keyB]['name'], information_proteins_hhsuite[keyB]['descr'], str(information_proteins_hhsuite[keyB]['prob']), str(information_proteins_hhsuite[keyB]['evalue']), str(information_proteins_hhsuite[keyB]['pident']), str(information_proteins_hhsuite[keyB]['pcover'])])
			except KeyError:
				try:
					print >>tablefile, "\t".join([equivalence[keyB], str(tempprotsdict[equivalence[keyB]]['begin']), str(tempprotsdict[equivalence[keyB]]['end']), str(tempprotsdict[equivalence[keyB]]['strand']), str(tempprotsdict[equivalence[keyB]]['length']), str(tempprotsdict[equivalence[keyB]]['isoelectricpoint']), str(tempprotsdict[equivalence[keyB]]['molweightkda']), str(tempprotsdict[equivalence[keyB]]['instability']), information_proteins_blast[equivalence[keyB]]['sseqid'], information_proteins_blast[equivalence[keyB]]['descr'], str(information_proteins_blast[equivalence[keyB]]['evalue']), str(information_proteins_blast[equivalence[keyB]]['pident']), str(information_proteins_blast[equivalence[keyB]]['pcover']), "None", "None", "NA", "NA", "NA", "NA"])
				except KeyError:
					try:
						print >>tablefile, "\t".join([equivalence[keyB], str(tempprotsdict[equivalence[keyB]]['begin']), str(tempprotsdict[equivalence[keyB]]['end']), str(tempprotsdict[equivalence[keyB]]['strand']), str(tempprotsdict[equivalence[keyB]]['length']), str(tempprotsdict[equivalence[keyB]]['isoelectricpoint']), str(tempprotsdict[equivalence[keyB]]['molweightkda']), str(tempprotsdict[equivalence[keyB]]['instability']), "None", "None", "NA", "NA", "NA", information_proteins_hhsuite[keyB]['name'], information_proteins_hhsuite[keyB]['descr'], str(information_proteins_hhsuite[keyB]['prob']), str(information_proteins_hhsuite[keyB]['evalue']), str(information_proteins_hhsuite[keyB]['pident']), str(information_proteins_hhsuite[keyB]['pcover'])])
					except KeyError:
						print >>tablefile, "\t".join([equivalence[keyB], str(tempprotsdict[equivalence[keyB]]['begin']), str(tempprotsdict[equivalence[keyB]]['end']), str(tempprotsdict[equivalence[keyB]]['strand']), str(tempprotsdict[equivalence[keyB]]['length']), str(tempprotsdict[equivalence[keyB]]['isoelectricpoint']), str(tempprotsdict[equivalence[keyB]]['molweightkda']), str(tempprotsdict[equivalence[keyB]]['instability']), "None", "None", "NA", "NA", "NA",  "None", "None", "NA", "NA", "NA", "NA"])

	# Algorithm of decisions (which one: BLAST/HHSUITE?)
	multipleprots = {}
	keylist = information_proteins_hhsuite.keys()
	keylist.sort()
	Hypotheticalpat = re.compile(r'([H|h]ypothetical|[U|u]ncharacteri[z|s]ed) protein')
	for keyB in keylist:
		keyB = keyB.replace(".hhr", "")
		singleprot = {}
		singleprot['name'] = equivalence[keyB]
		if (equivalence[keyB] in information_proteins_blast) and (keyB in information_proteins_hhsuite):
			if re.match(Hypotheticalpat, information_proteins_blast[equivalence[keyB]]['descr']):
				try:
					if re.match(Hypotheticalpat, information_proteins_hhsuite[keyB]['descr']):
						singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
					else:
						singleprot['descr'] = ' '.join(["Putative", information_proteins_hhsuite[keyB]['descr']])
				except KeyError:
					singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
			else:
				try:
					if (float(information_proteins_blast[equivalence[keyB]]['pident'])>float(information_proteins_hhsuite[keyB]['pident'])) and (float(information_proteins_blast[equivalence[keyB]]['pcover'])>float(information_proteins_hhsuite[keyB]['pcover'])):
						singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
					elif (float(information_proteins_blast[equivalence[keyB]]['pident'])<float(information_proteins_hhsuite[keyB]['pident'])) and (float(information_proteins_blast[equivalence[keyB]]['pcover'])<float(information_proteins_hhsuite[keyB]['pcover'])):
						singleprot['descr'] = information_proteins_hhsuite[keyB]['descr']
					elif (float(information_proteins_blast[equivalence[keyB]]['pident'])>float(information_proteins_hhsuite[keyB]['pident'])) and (float(information_proteins_blast[equivalence[keyB]]['pcover'])<float(information_proteins_hhsuite[keyB]['pcover'])):
						if (float(information_proteins_blast[equivalence[keyB]]['pident'])-float(information_proteins_hhsuite[keyB]['pident']) >= args.diffid):
							singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
						else:
							singleprot['descr'] = information_proteins_hhsuite[keyB]['descr']
					else:
						if (float(information_proteins_hhsuite[keyB]['pident'])-float(information_proteins_blast[equivalence[keyB]]['pident']) >= args.diffid):
							singleprot['descr'] = information_proteins_hhsuite[keyB]['descr']
						else:
							singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
				except KeyError:
					try:
						if (float(information_proteins_blast[equivalence[keyB]]['pcover'])>float(information_proteins_hhsuite[keyB]['pcover'])):
							singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
						else:
							singleprot['descr'] = information_proteins_hhsuite[keyB]['descr']
					except KeyError:
						singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
		elif equivalence[keyB] in information_proteins_blast:
			singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
		elif keyB in information_proteins_hhsuite:
			try:
				singleprot['descr'] = information_proteins_hhsuite[keyB]['descr']
			except KeyError:
				singleprot['descr'] = 'Hypothetical protein'
		else:
			singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
		multipleprots[keyB] = singleprot

	# (DEBUG) Solution of the decision algorithm
	#print "NAME\tDESCRIPTION"
	#for keyB in sorted(multipleprots):
	#	print "\t".join([multipleprots[keyB]['name'], multipleprots[keyB]['descr']])

	#Storing protein information in memory
	with open("temp.faa", "rU") as protsfile:
		protsdict = {}
		for protseq in SeqIO.parse(protsfile, "fasta"):
			indprot = {}
			dataprot = protseq.description.split(' # ')
			indprot['translation'] = protseq.seq
			indprot['protein_id'] = dataprot[0]
			indprot['strand'] = int(dataprot[3])
			indprot['begin'] = int(dataprot[1])-1
			indprot['end'] = int(dataprot[2])
			for keyOmega in sorted(multipleprots):
				if multipleprots[keyOmega]['name'] == dataprot[0]:
					indprot['product'] = multipleprots[keyOmega]['descr']
			protsdict[dataprot[0]] = indprot

	# Predicting the tRNA sequences
	print "Running ARAGORN to predict tRNA-like sequences in %s" % newfile
	genetictable = "-gc%s" % str(args.gcode)
	with open("trnafile.fasta", "w") as trnafile:
	    subprocess.call(["aragorn", "-l", "-fon", genetictable, newfile], stdout=trnafile)
	print "DONE. ARAGORN was done to predict tRNA sequences\n\n"

	#Storing tRNA information in memory
	with open("trnafile.fasta", "rU") as trnafile:
		tRNAdict = {}
		for tRNAseq in SeqIO.parse(trnafile, "fasta"):
			indtRNA = {}
			tRNA_information = tRNAseq.description.split(" ")
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
			tRNAdict[tRNAseq.id] = indtRNA

	# Predicting repeats
	print "Predicting repeats in the sequences using TRF and IRF"
	with open("/dev/null", "w") as stderr:
		subprocess.call(["trf", newfile, "2", "7", "7", "80", "10", "50", "500", "-h"], stderr=stderr)
		os.rename("%s.2.7.7.80.10.50.500.dat" % newfile, "trf_temp.dat")
	with open("/dev/null", "w") as stderr:
		subprocess.call(["irf", newfile, "2", "3", "5", "80", "10", "40", "500000", "10000", "-d", "-h"], stderr=stderr)
		os.rename("%s.2.3.5.80.10.40.500000.10000.dat" % newfile, "irf_temp.dat")
	print "DONE. TRF and IRF were launched to predict repeats"

	# Storing repeats information
	information_TRF = {}
	with open("trf_temp.dat", "rU") as trfile:
		information_tandem_repeat = {}
		for line in trfile:
			try:
				patT = re.compile('^(\d+)\s(\d+)\s\d+\s\d+\.\d+\s')
				start, end = re.match(patT, line).groups()
				combinedinfo = "%s_%s" % (str(start), str(end))
			except AttributeError:
				continue
			else:
				information_tandem_repeat['start'] = start
				information_tandem_repeat['end'] = end
			information_TRF[combinedinfo] = information_tandem_repeat

	information_IRF = {}
	with open("irf_temp.dat", "rU") as irfile:
		information_inverted_repeat = {}
		for line in irfile:
			try:
				patI = re.compile('^(\d+)\s(\d+)\s\d+\s\d+\s\d+')
				start, end = re.match(patI, line).groups()
				combinedinfo = "%s_%s" % (str(start), str(end))
			except AttributeError:
				continue
			else:
				information_inverted_repeat['start'] = start
				information_inverted_repeat['end'] = end
			information_IRF[combinedinfo] = information_inverted_repeat

	# Creating a new Genbank and GFF file
	newgbk = "%s.gbk" % newfile
	with open(newfile, "rU") as basefile, open(newgbk, "w"):
		for record in SeqIO.parse(basefile, "fasta", IUPAC.ambiguous_dna):
			whole_sequence = SeqRecord(record.seq)
			whole_sequence.id = re.sub("\_\d+$", "", str(record.id))
			whole_sequence.annotations['data_file_division'] = args.typedata.upper()
			whole_sequence.annotations['date'] = strftime("%d-%b-%Y").upper()
			whole_sequence.annotations['topology'] = genomeshape['genomeshape']
			for protein in sorted(protsdict, key = stringSplitByNumbers):
				start_pos = SeqFeature.ExactPosition(protsdict[protein]['begin'])
				end_pos = SeqFeature.ExactPosition(protsdict[protein]['end'])
				feature_location = SeqFeature.FeatureLocation(start_pos, end_pos, strand=protsdict[protein]['strand'])
				new_data_gene = SeqFeature.SeqFeature(feature_location, type = "gene", strand = protsdict[protein]['strand'])
				whole_sequence.features.append(new_data_gene)
				qualifiers = [('product', protsdict[protein]['product']), ('protein_id', protsdict[protein]['protein_id']), ('translation', protsdict[protein]['translation'])]
				feature_qualifiers = OrderedDict(qualifiers)
				new_data_cds = SeqFeature.SeqFeature(feature_location, type = "CDS", strand = protsdict[protein]['strand'], qualifiers = feature_qualifiers)
				whole_sequence.features.append(new_data_cds)
			for tRNA in sorted(tRNAdict, key = stringSplitByNumbers):
				start_pos = SeqFeature.ExactPosition(tRNAdict[tRNA]['begin'])
				end_pos = SeqFeature.ExactPosition(tRNAdict[tRNA]['end'])
				feature_location = SeqFeature.FeatureLocation(start_pos, end_pos, strand=tRNAdict[tRNA]['strand'])
				new_data_gene = SeqFeature.SeqFeature(feature_location, type = "gene", strand = tRNAdict[tRNA]['strand'])
				whole_sequence.features.append(new_data_gene)
				qualifiers = [('product', tRNAdict[tRNA]['product'])]
				feature_qualifiers = OrderedDict(qualifiers)
				new_data_tRNA = SeqFeature.SeqFeature(feature_location, type = "tRNA", strand = tRNAdict[tRNA]['strand'], qualifiers = feature_qualifiers)
				whole_sequence.features.append(new_data_tRNA)
			for tandem in sorted(information_TRF, key = stringSplitByNumbers):
				start_pos = SeqFeature.ExactPosition(information_TRF[tandem]['start'])
				end_pos = SeqFeature.ExactPosition(information_TRF[tandem]['end'])
				feature_location = SeqFeature.FeatureLocation(start_pos, end_pos)
				qualifiers = [('rpt_type', 'direct')]
				feature_qualifiers = OrderedDict(qualifiers)
				new_data_tandemrepeat = SeqFeature.SeqFeature(feature_location, type = "repeat_region", qualifiers = feature_qualifiers)
				whole_sequence.features.append(new_data_tandemrepeat)
			for inverted in sorted(information_IRF, key = stringSplitByNumbers):
				start_pos = SeqFeature.ExactPosition(information_IRF[inverted]['start'])
				end_pos = SeqFeature.ExactPosition(information_IRF[inverted]['end'])
				feature_location = SeqFeature.FeatureLocation(start_pos, end_pos)
				qualifiers = [('rpt_type', 'inverted')]
				feature_qualifiers = OrderedDict(qualifiers)
				new_data_invertedrepeat = SeqFeature.SeqFeature(feature_location, type = "repeat_region", qualifiers = feature_qualifiers)
				whole_sequence.features.append(new_data_invertedrepeat)
			SeqIO.write(whole_sequence, newgbk, "genbank")

	if args.gffprint==True:
		newgff = "%s.gff" % newfile
		with open(newgff, "w") as outgff, open(newgbk, "rU") as ingbk:
			GFF.write(SeqIO.parse(ingbk, "genbank"), outgff)

	# Removing intermediate files
	os.remove(newfile)
	os.remove("temporal_circular.fasta")
	os.remove("temp.faa")
	os.remove("temp_blast.csv")
	os.remove("trnafile.fasta")
	os.remove("trf_temp.dat")
	os.remove("irf_temp.dat")
	for f in glob.glob("SEQ*"):
		os.remove(f)

# Joining all GENBANK files into one
listgbk = sorted(glob.glob('CONTIG_*.gbk'))
gbkoutputfile = "%s.gbk" % root_output
with open(gbkoutputfile, 'w') as finalgbk:
	for fname in listgbk:
		with open(fname) as infile:
			for line in infile:
				finalgbk.write(line)

for tempgbk in glob.glob("CONTIG_*.gbk"):
	os.remove(tempgbk)

# Joining all GFF files into one
if args.gffprint==True:
	listgff = sorted(glob.glob('CONTIG_*.gff'))
	gffoutputfile = "%s.gff" % root_output
	with open(gffoutputfile, 'w') as finalgff:
		for fname in listgff:
			with open(fname) as infile:
				for line in infile:
					finalgff.write(line)
	for tempgff in glob.glob("CONTIG_*.gff"):
		os.remove(tempgff)

# Joining all TABLE files into one
listcsv = sorted(glob.glob('CONTIG_*.csv'))
tbloutputfile = "%s.csv" % root_output
with open(tbloutputfile, 'w') as finaltable:
	for fname in listcsv:
		with open(fname) as infile:
			for line in infile:
				finaltable.write(line)

for temptbl in glob.glob("CONTIG_*.csv"):
	os.remove(temptbl)

# Preparing sequences for GenBank submission (Original code from Wan Yu's gbk2tbl.py script [https://github.com/wanyuac/BINF_toolkit/blob/master/gbk2tbl.py])

allowed_qualifiers = ['locus_tag', 'gene', 'product', 'pseudo', 'protein_id', 'gene_desc', 'old_locus_tag', 'note', 'inference', 'organism', 'mol_type', 'strain', 'sub_species', 'isolation-source', 'country']
newfastafile = "%s.fasta" % root_output
newtablefile = "%s.tbl" % root_output
with open(args.modifiers, "rU") as modifiers, open(gbkoutputfile, "r") as genbank_fh, open(newfastafile, "w") as fasta_fh, open(newtablefile, "w") as feature_fh: 
	info = modifiers.readline()
	wholelist = list(SeqIO.parse(genbank_fh, 'genbank'))
	for record in wholelist:
		if len(record) <= args.mincontigsize:
			print "WARNING: Skipping small contig %s" % rec.id
			continue
		record.description = "%s %s" % (record.id, info)
		SeqIO.write([record], fasta_fh, 'fasta')
		print >>feature_fh, '>Feature %s' % (record.name)
		for line in record.features:
			if line.strand == 1:
				print >>feature_fh, '%d\t%d\t%s' % (line.location.nofuzzy_start + 1, line.location.nofuzzy_end, line.type)
			else:
				print >>feature_fh, '%d\t%d\t%s' % (line.location.nofuzzy_end, line.location.nofuzzy_start + 1, line.type)
			for (key, values) in line.qualifiers.iteritems():
				if key not in allowed_qualifiers:
					continue
				for v in values:
					print >>feature_fh, '\t\t\t%s\t%s' % (key, v)

# Final statement
print "Genome annotation done!"
print "The GenBank file is %s" % gbkoutputfile
if args.gffprint==True:
	print "The GFF3 file is %s" % gffoutputfile
print "The table file for GenBank submission is %s" % tbloutputfile
print "The FASTA file for GenBank submission is %s" % newfastafile
print "The table file with all protein information is %s" % newtablefile
