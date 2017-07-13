#!/usr/bin/env python

# -*- coding: utf-8 -*-

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
from __future__ import print_function
import argparse
import csv
import fileinput
import fractions
import glob
import os
import re
import sys
import shutil
import subprocess
from BCBio import GFF
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import OrderedDict, defaultdict
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

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

#def find(name, path):
#	for root, dirs, files in os.walk(path):
#		if name in files:
#			return os.path.join(root, name)

def stringSplitByNumbers(x):
	r = re.compile('(\d+)')
	l = r.split(x)
	return [int(y) if y.isdigit() else y for y in l]

# Defining the program version
version = "0.8.1"

# Processing the parameters
parser = argparse.ArgumentParser(description='Virannot is a automatic de novo viral genome annotator.')
basic_group = parser.add_argument_group('Basic options for virannot [REQUIRED]')

basic_group.add_argument("--input", dest="inputfile", type=str, required=True, help='Input file as a FASTA file', metavar="FASTAFILE")
basic_group.add_argument("--rfamdb", dest="rfamdatabase", type=str, required=True, help='RFAM Database that will be used for the ribosomal RNA prediction. RFAMDB should be in the format "/full/path/to/rfamdb/Rfam.cm" and must be compressed accordingly (see INFERNAL manual) before running the script.', metavar="RFAMDB")
basic_group.add_argument("--modifiers", dest="modifiers", type=str, required=True, help='Input file as a plain text file with the modifiers per every FASTA header according to SeqIn (https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html). All modifiers must be written in a single line and are separated by a single space character. No space should be placed besides the = sign. For example: [organism=Serratia marcescens subsp. marcescens] [sub-species=marcescens] [strain=AH0650_Sm1] [topology=linear] [moltype=DNA] [tech=wgs] [gcode=11] [country=Australia] [isolation-source=sputum]. This line will be copied and printed along with the record name as the definition line of every contig sequence.', metavar="TEXTFILE")

advanced_group = parser.add_argument_group('Advanced options for virannot [OPTIONAL]')
advanced_group.add_argument("--readlength", dest="read_length", type=int, default=101, help='Read length for the circularity prediction (default: 101 bp)', metavar="INT")
advanced_group.add_argument("--out", dest="rootoutput", type=str, help='Name of the outputs files (without extension)', metavar="OUTPUTNAME")
advanced_group.add_argument("--locus", dest="locus", type=str, default='LOC', help='Name of the sequences. If the input is a multiFASTA file, please put a general name as the program will add the number of the contig at the end of the name (Default: %(default)s)', metavar="STRING")
advanced_group.add_argument("--gff", dest="gffprint", action='store_true', default=False, help='Printing the output as GFF3 file (Default: False)')
advanced_group.add_argument("--threads", dest="ncpus", default=1, help='Number of threads/cpus (Default: %(default)s cpu)', metavar="INT")
advanced_group.add_argument("--nohmmer", dest="nohmmer", action='store_true', default=False, help='Running only BLAST to predict protein function. (Default: False)')
advanced_group.add_argument("--noblast", dest="noblast", action='store_true', default=False, help='Running DIAMOND instead of BLAST to predict protein function according to homology. This will be less sensitive but faster than BLAST. (Default: False)')
advanced_group.add_argument("--blastevalue", dest="blastevalue", default=0.00001, help='Blast e-value threshold (Default: 0.00001)', metavar="FLOAT")
basic_group.add_argument("--blastdb", dest="blastdatabase", type=str, help='BLAST Database that will be used for the protein function prediction. The database must be an amino acid one, not  nucleotidic', metavar="BLASTDB")
basic_group.add_argument("--diamonddb", dest="diamonddatabase", type=str, help='DIAMOND Database that will be used for the protein function prediction. The database must be created from a amino acid FASTA file as indicated in https://github.com/bbuchfink/diamond. This argument is mandatory when --ultrafast option is enabled', metavar="DIAMONDDB")
advanced_group.add_argument("--hmmdb", dest="hmmdatabase", type=str, help='PHMMER Database that will be used for the protein function prediction according to Hidden Markov Models. In this case, HMMDB must be in FASTA format (e.g. UniProt: "', metavar="HMMDB")
advanced_group.add_argument("--hmmerevalue", dest="hmmerevalue", default=0.001, help='PHMMER e-value threshold (Default: 0.001)', metavar="FLOAT")

type_choices = {'BCT': 'Prokaryotic chromosome', 'CON': 'Contig', 'PHG': 'Phages', 'VRL': 'Eukaryotic/Archaea virus'}
type_help = ('GenBank Division: One of the following codes - {0}. (Default: %(default)s)'.format(', '.join('{0} ({1})'.format(k, v) for k, v in type_choices.items())))
advanced_group.add_argument("--typedata", dest="typedata", type=str, default='CON', help=type_help, metavar="BCT|CON|VRL|PHG")

gcode_choices = {'1': 'Standard genetic code [Eukaryotic]', '2': 'Vertebrate mitochondrial code', '3': 'Yeast mitochondrial code', '4': 'Mycoplasma/Spiroplasma and Protozoan/mold/coelenterate mitochondrial code', '5': 'Invertebrate mitochondrial code', '6': 'Ciliate, dasycladacean and hexamita nuclear code', '9': 'Echinoderm/flatworm mitochondrial code', '10': 'Euplotid nuclear code', '11': 'Bacteria/Archaea/Phages/Plant plastid', '12': 'Alternative yeast nuclear code', '13': 'Ascidian mitochondrial code', '14': 'Alternative flatworm mitochondrial code', '16': 'Chlorophycean mitochondrial code', '21': 'Trematode mitochondrial code', '22': 'Scedenesmus obliquus mitochondrial code', '23': 'Thraustochytrium mitochondrial code', '24': 'Pterobranquia mitochondrial code', '25': 'Gracilibacteria & Candidate division SR1', '26': 'Pachysolen tannophilus nuclear code', '27': 'Karyorelict nuclear code', '28': 'Condylostoma nuclear code', '29': 'Mesodinium nuclear code', '30': 'Peritrich nuclear code', '31': 'Blastocrithidia nuclear code'}
gcode_help = ('Number of GenBank translation table. At this moment, the available options are {0}. (Default: %(default)s)'.format('{}'.format(', '.join('{0} ({1})'.format(k, v) for k, v in gcode_choices.items()))))
advanced_group.add_argument("--gcode", dest="gcode", type=str, default='11', help=gcode_help, metavar="NUMBER")

advanced_group.add_argument('--mincontigsize', dest="mincontigsize", type=int, default = 200, help = 'Minimum contig length to be considered in the final files (Default: 200 bp)', metavar="INT")
advanced_group.add_argument("--idthr", dest="idthreshold", default=50.00, help='ID threshold (Default: 50.0)', metavar="FLOAT")
advanced_group.add_argument("--coverthr", dest="covthreshold", default=50.00, help='Coverage threshold (Default: 50.0)', metavar="FLOAT")
advanced_group.add_argument("--diffid", dest="diffid", default=5.00, help='Max allowed difference between the ID percentages of BLAST and HMMER. (Default = 5.00; Not recommended to change such value)', metavar="FLOAT (>0.01)")
advanced_group.add_argument("--minrepeat", dest="minrepeat", type=int, default=16, help="Minimum repeat length for CRISPR detection (Default: 16)", metavar="INT")
advanced_group.add_argument("--maxrepeat", dest="maxrepeat", type=int, default=64, help="Maximum repeat length for CRISPR detection (Default: 64)")
advanced_group.add_argument("--minspacer", dest="minspacer", type=int, default=8, help="Minimum spacer length for CRISPR detection (Default: 8)")
advanced_group.add_argument("--maxspacer", dest="maxspacer", type=int, default=64, help="Maximum spacer length for CRISPR detection (Default: 64)")
advanced_group.add_argument("--blastexh", dest="blastexh", action='store_true', default=False, help='Use of exhaustive BLAST to predict the proteins by homology according to Fozo et al. (2010) Nucleic Acids Res (Default=FALSE)')

args = parser.parse_args()

root_output = args.rootoutput
if not root_output:
	root_output = '{}_annotated'.format(os.path.splitext(args.inputfile)[0])

if args.noblast == False and args.blastdatabase == None:
    sys.exit('You MUST specify BLAST database using the parameter --blastdb if you are not using --noblast option')

if args.noblast == True and args.diamonddatabase == None:
    sys.exit('You MUST specify DIAMOND database using the parameter --diamonddb if you are using --noblast option')

if args.nohmmer == False and args.noblast == False and args.hmmdatabase == None:
		sys.exit('You MUST specify HMMER database using the parameter --hmmdb if you are not using --nohmmer option')

# Printing the header of the program 
eprint("This is VirAnnot %s" % str(version))
eprint("Written by Enrique Gonzalez Tortuero & Vimalkumar Velayudhan")
eprint("Homepage is https://github.com/EGTortuero/virannot")
eprint("Local time: ", strftime("%a, %d %b %Y %H:%M:%S"))
eprint("\n\n")

# checking the presence of the programs in the system

if not cmd_exists("lastz")==True:
	sys.exit("You must install LASTZ before using this script")
elif not cmd_exists("cmscan")==True:
	sys.exit("You must install INFERNAL before using this script")
elif not cmd_exists("prodigal")==True:
	sys.exit("You must install prodigal before using this script")
elif not cmd_exists("parallel")==True:
	sys.exit("You must install GNU Parallel before using this script")
elif not cmd_exists("blastp")==True:
	sys.exit("You must install BLAST before using this script")
elif not cmd_exists("diamond")==True:
	sys.exit("You must install DIAMOND before using this script")
elif not cmd_exists("phmmer")==True:
	sys.exit("You must install HMMER 3 before using this script")
elif not cmd_exists("aragorn")==True:
	sys.exit("You must install ARAGORN before using this script")
elif not cmd_exists("pilercr")==True:
	sys.exit("You must install PILERCR before using this script")
elif not cmd_exists("trf")==True:
	sys.exit("You must install Tandem Repeats Finder before using this script")
elif not cmd_exists("irf")==True:
	sys.exit("You must install Inverted Repeats Finder before using this script")

eprint("Data type is {0} and GenBank translation table no is {1}\n".format(args.typedata, args.gcode))

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
			eprint("WARNING: The name of the sequence %s was corrected as %s" % (original_name, record.id))
		SeqIO.write(record, corrected, "fasta")
	os.remove(filename)

for newfile in sorted(glob.glob("CONTIG_*.fna")):

	# Predicting the shape of the contig (code based on Alex Crits-Christoph's find_circular.py script [https://github.com/alexcritschristoph/VICA/blob/master/find_circular.py])
	eprint("Predicting the shape of the contig using LASTZ")
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
		eprint("LASTZ predicted that %s is %s\n" % (newfile, genomeshape['genomeshape']))

	# Predicting genes using PRODIGAL
	eprint("\nRunning Prodigal to predict the genes in %s" % newfile)
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

	if args.blastexh==True:
		eprint("Running BLAST to predict the genes according to homology inference in %s using exhaustive mode (see Fozo et al. (2010) Nucleic Acids Res for details)" % newfile)
		subprocess.call(['blastp', '-query', "temp.faa", '-db', args.blastdatabase, '-evalue', str(args.blastevalue), '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle', '-out', 'temp_blast.csv', '-max_target_seqs', '10', '-word_size', '2', '-gapopen', '8', '-gapextend', '2', '-matrix', '"PAM70"', '-comp_based_stats', '"0"', "-num_threads", str(args.ncpus)])
	elif args.noblast==True:
		eprint("Running DIAMOND to predict the genes according to homology inference in %s using default parameters" % newfile)
		subprocess.call(['diamond', 'blastp', '-q', "temp.faa", '-d', args.diamonddatabase, '-e', str(args.blastevalue), '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qstart', 'qend', 'evalue', 'bitscore', 'stitle', '-o', 'temp_blast.csv', '-k', '10', "-p", str(args.ncpus), '--quiet'])
	else:
		eprint("Running BLAST to predict the genes according to homology inference in %s using default parameters" % newfile)
		subprocess.call(['blastp', '-query', "temp.faa", '-db', args.blastdatabase, '-evalue', str(args.blastevalue), '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle', '-out', 'temp_blast.csv', '-max_target_seqs', '10', "-num_threads", str(args.ncpus)])
	eprint("Done. BLASTp was done to predict the genes by homology\n")

	blast_log = "%s.blast.log" % newfile
	shutil.copyfile("temp_blast.csv", blast_log)

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

	## Predicting the function of the proteins based on HMM predictions using phmmer
	if args.nohmmer == False and args.noblast == False:
		with open("commands.sh", "w") as commands:
			for singleprot in sorted(glob.glob("SEQ_*.faa")):
				hhmtable = "%s.tbl" % singleprot
				eprint("Creating file to run parallel PHMMER")
				eprint("Adding %s to run PHMMER." % singleprot)
				line2write = ' '.join(["phmmer", "--cpu", "1", "--domtblout", hhmtable, "-E", str(args.hmmerevalue), "-o", "/dev/null", singleprot, args.hmmdatabase, '\n'])
				commands.write(line2write)

		eprint("Running parallel PHMMER")
		subprocess.call(['parallel', '-j', str(args.ncpus)], stdin=open('commands.sh'))
		eprint("Done. PHMMER was done to predict the function of the genes according to Hidden Markov Models\n")
		os.remove("commands.sh")

		# Parsing the results from HMMER
		information_proteins_hmmer = {}
		for singletbl in sorted(glob.glob("*.faa.tbl")):
			rootname = singletbl.replace(".faa.tbl", "")
			with open(singletbl) as tblfile:
				infoprot_hmmer = {}
				for line in tblfile:
					if line.startswith("#"):
						continue
					else:
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
								if float(protarea) >= float(args.covthreshold) and (float(evaluehh) <= float(args.hmmerevalue)) and (float(pident) >= 50.00):
									infoprot_hmmer['lociname'] = lociname
									infoprot_hmmer['name'] = matchname
									infoprot_hmmer['evalue'] = float(evaluehh)
									infoprot_hmmer['pcover'] = float(protarea)
									infoprot_hmmer['pident'] = float(pident)
									infoprot_hmmer['descr'] = description
								else:
									try:
										if (float(protarea) >= float(args.covthreshold)) and (float(protarea) >= float(infoprot_hmmer['pcover'])) and (float(evaluehh) <= float(args.hmmerevalue)) and (float(evaluehh) <= float(infoprot_hmmer['evalue'])) and (float(pident) >= 50.00) and (float(pident) >= infoprot_hmmer['pident']) :
											infoprot_hmmer['lociname'] = lociname
											infoprot_hmmer['name'] = matchname
											infoprot_hmmer['evalue'] = float(evaluehh)
											infoprot_hmmer['pcover'] = float(protarea)
											infoprot_hmmer['pident'] = float(pident)
											infoprot_hmmer['descr'] = description
									except KeyError:
											continue
							else:
								if (float(protarea) >= float(args.covthreshold)) and (float(protarea) >= float(infoprot_hmmer['pcover'])) and (float(evaluehh) <= float(args.hmmerevalue)) and (float(evaluehh) <= float(infoprot_hmmer['evalue'])) and (float(pident) >= 50.00):
									infoprot_hmmer['lociname'] = lociname
									infoprot_hmmer['name'] = matchname
									infoprot_hmmer['evalue'] = float(evaluehh)
									infoprot_hmmer['pcover'] = float(protarea)
									infoprot_hmmer['pident'] = float(pident)
									infoprot_hmmer['descr'] = description
				information_proteins_hmmer[rootname] = infoprot_hmmer

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
		if args.noblast == False and args.nohmmer == False:
			print("\t".join(["Identifier", "Start", "Stop", "Strand", "size_aa", "pI", "Mol_weight_kDa", "Instability_index", "ID_BLAST", "Descr_BLAST", "evalue_BLAST", "%ID_BLAST", "%Cover_BLAST", "ID_HMMER", "Descr_HMMER", "evalue_HMMER", "%ID_HMMER", "%Cover_HMMER"]), file=tablefile)
			keylist = information_proteins_hmmer.keys()
			keylist.sort()
			for keyB in keylist:
				keyB = keyB.replace(".faa.tbl", "")
				try:
					print("\t".join([equivalence[keyB], str(tempprotsdict[equivalence[keyB]]['begin']), str(tempprotsdict[equivalence[keyB]]['end']), str(tempprotsdict[equivalence[keyB]]['strand']), str(tempprotsdict[equivalence[keyB]]['length']), str(tempprotsdict[equivalence[keyB]]['isoelectricpoint']), str(tempprotsdict[equivalence[keyB]]['molweightkda']), str(tempprotsdict[equivalence[keyB]]['instability']), information_proteins_blast[equivalence[keyB]]['sseqid'], information_proteins_blast[equivalence[keyB]]['descr'], str(information_proteins_blast[equivalence[keyB]]['evalue']), str(information_proteins_blast[equivalence[keyB]]['pident']), str(information_proteins_blast[equivalence[keyB]]['pcover']), information_proteins_hmmer[keyB]['name'], information_proteins_hmmer[keyB]['descr'], str(information_proteins_hmmer[keyB]['evalue']), str(information_proteins_hmmer[keyB]['pident']), str(information_proteins_hmmer[keyB]['pcover'])]), file=tablefile)
				except KeyError:
					try:
						print("\t".join([equivalence[keyB], str(tempprotsdict[equivalence[keyB]]['begin']), str(tempprotsdict[equivalence[keyB]]['end']), str(tempprotsdict[equivalence[keyB]]['strand']), str(tempprotsdict[equivalence[keyB]]['length']), str(tempprotsdict[equivalence[keyB]]['isoelectricpoint']), str(tempprotsdict[equivalence[keyB]]['molweightkda']), str(tempprotsdict[equivalence[keyB]]['instability']), information_proteins_blast[equivalence[keyB]]['sseqid'], information_proteins_blast[equivalence[keyB]]['descr'], str(information_proteins_blast[equivalence[keyB]]['evalue']), str(information_proteins_blast[equivalence[keyB]]['pident']), str(information_proteins_blast[equivalence[keyB]]['pcover']), "None", "None", "NA", "NA", "NA"]), file=tablefile)
					except KeyError:
						try:
							print("\t".join([equivalence[keyB], str(tempprotsdict[equivalence[keyB]]['begin']), str(tempprotsdict[equivalence[keyB]]['end']), str(tempprotsdict[equivalence[keyB]]['strand']), str(tempprotsdict[equivalence[keyB]]['length']), str(tempprotsdict[equivalence[keyB]]['isoelectricpoint']), str(tempprotsdict[equivalence[keyB]]['molweightkda']), str(tempprotsdict[equivalence[keyB]]['instability']), "None", "None", "NA", "NA", "NA", information_proteins_hmmer[keyB]['name'], information_proteins_hmmer[keyB]['descr'], str(information_proteins_hmmer[keyB]['evalue']), str(information_proteins_hmmer[keyB]['pident']), str(information_proteins_hmmer[keyB]['pcover'])]), file=tablefile)
						except KeyError:
							print("\t".join([equivalence[keyB], str(tempprotsdict[equivalence[keyB]]['begin']), str(tempprotsdict[equivalence[keyB]]['end']), str(tempprotsdict[equivalence[keyB]]['strand']), str(tempprotsdict[equivalence[keyB]]['length']), str(tempprotsdict[equivalence[keyB]]['isoelectricpoint']), str(tempprotsdict[equivalence[keyB]]['molweightkda']), str(tempprotsdict[equivalence[keyB]]['instability']), "None", "None", "NA", "NA", "NA",  "None", "None", "NA", "NA", "NA"]), file=tablefile)
		else:
			print("\t".join(["Identifier", "Start", "Stop", "Strand", "size_aa", "pI", "Mol_weight_kDa", "Instability_index", "ID_BLAST", "Descr_BLAST", "evalue_BLAST", "%ID_BLAST", "%Cover_BLAST"]), file=tablefile)
			keylist = equivalence.values()
			keylist.sort()
			for keyB in keylist:
				try:
					print("\t".join([keyB, str(tempprotsdict[keyB]['begin']), str(tempprotsdict[keyB]['end']), str(tempprotsdict[keyB]['strand']), str(tempprotsdict[keyB]['length']), str(tempprotsdict[keyB]['isoelectricpoint']), str(tempprotsdict[keyB]['molweightkda']), str(tempprotsdict[keyB]['instability']), information_proteins_blast[keyB]['sseqid'], information_proteins_blast[keyB]['descr'], str(information_proteins_blast[keyB]['evalue']), str(information_proteins_blast[keyB]['pident']), str(information_proteins_blast[keyB]['pcover'])]), file=tablefile)
				except KeyError:
					print("\t".join([keyB, str(tempprotsdict[keyB]['begin']), str(tempprotsdict[keyB]['end']), str(tempprotsdict[keyB]['strand']), str(tempprotsdict[keyB]['length']), str(tempprotsdict[keyB]['isoelectricpoint']), str(tempprotsdict[keyB]['molweightkda']), str(tempprotsdict[keyB]['instability']), "None", "None", "NA", "NA", "NA"]), file=tablefile)

	# Algorithm of decisions (which one: BLAST/HMMER?)
	multipleprots = {}
	Hypotheticalpat = re.compile(r'(((H|h)ypothetical)|((U|u)ncharacteri(z|s)ed)) protein')
	if args.nohmmer == False and args.noblast == False:
		keylist = information_proteins_hmmer.keys()
		keylist.sort()
		for keyB in keylist:
			keyB = keyB.replace(".faa.tbl", "")
			singleprot = {}
			singleprot['name'] = equivalence[keyB]
			if (equivalence[keyB] in information_proteins_blast) and (keyB in information_proteins_hmmer):
				if re.match(Hypotheticalpat, information_proteins_blast[equivalence[keyB]]['descr']):
					try:
						if re.match(Hypotheticalpat, information_proteins_hmmer[keyB]['descr']):
							singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
						else:
							singleprot['descr'] = information_proteins_hmmer[keyB]['descr']
					except KeyError:
						singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
				else:
					try:
						if (float(information_proteins_blast[equivalence[keyB]]['pident'])>float(information_proteins_hmmer[keyB]['pident'])) and (float(information_proteins_blast[equivalence[keyB]]['pcover'])>float(information_proteins_hmmer[keyB]['pcover'])):
							singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
						elif (float(information_proteins_blast[equivalence[keyB]]['pident'])<float(information_proteins_hmmer[keyB]['pident'])) and (float(information_proteins_blast[equivalence[keyB]]['pcover'])<float(information_proteins_hmmer[keyB]['pcover'])):
							singleprot['descr'] = information_proteins_hmmer[keyB]['descr']
						elif (float(information_proteins_blast[equivalence[keyB]]['pident'])>float(information_proteins_hmmer[keyB]['pident'])) and (float(information_proteins_blast[equivalence[keyB]]['pcover'])<float(information_proteins_hmmer[keyB]['pcover'])):
							if (float(information_proteins_blast[equivalence[keyB]]['pident'])-float(information_proteins_hmmer[keyB]['pident']) >= args.diffid):
								singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
							else:
								singleprot['descr'] = information_proteins_hmmer[keyB]['descr']
						else:
							if (float(information_proteins_hmmer[keyB]['pident'])-float(information_proteins_blast[equivalence[keyB]]['pident']) >= args.diffid):
								singleprot['descr'] = information_proteins_hmmer[keyB]['descr']
							else:
								singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
					except KeyError:
						try:
							if (float(information_proteins_blast[equivalence[keyB]]['pcover'])>float(information_proteins_hmmer[keyB]['pcover'])):
								singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
							else:
								singleprot['descr'] = information_proteins_hmmer[keyB]['descr']
						except KeyError:
							singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
			elif equivalence[keyB] in information_proteins_blast:
				singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
			elif keyB in information_proteins_hmmer:
				try:
					singleprot['descr'] = information_proteins_hmmer[keyB]['descr']
				except KeyError:
					singleprot['descr'] = 'Hypothetical protein'
			else:
				singleprot['descr'] = information_proteins_blast[equivalence[keyB]]['descr']
			multipleprots[keyB] = singleprot
	else:
		keylist = equivalence.values()
		keylist.sort()
		for keyB in keylist:
			singleprot = {}
			singleprot['name'] = keyB
			try:
				if information_proteins_blast[keyB]['descr'] == None:
					singleprot['descr'] = 'Hypothetical protein'
				else:
					singleprot['descr'] = information_proteins_blast[keyB]['descr']
			except KeyError:
				singleprot['descr'] = 'Hypothetical protein'
			multipleprots[keyB] = singleprot

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

	# Predicting the rRNA sequences
	with open(newfile, "rU") as targetfasta, open("/dev/null", "w") as apocalypse:
		eprint("Running INFERNAL+RFAM to predict rRNA-like sequences in %s" % newfile)
		subprocess.call(["cmscan", "--rfam", "--cut_ga", "--nohmmonly", "--tblout", "rrnafile.csv", "--cpu", args.ncpus, args.rfamdatabase, newfile], stdout=apocalypse)

	#Storing rRNA information in memory
	with open("rrnafile.csv", "rU") as rrnafile:
		namedict = {"SSU_rRNA_archaea": "16s_rRNA", "SSU_rRNA_bacteria": "16s_rRNA", "SSU_rRNA_eukarya": "18s_rRNA", "SSU_rRNA_microsporidia": "16s_rRNA", "LSU_rRNA_archaea": "23s_rRNA", "LSU_rRNA_bacteria": "23s_rRNA", "LSU_rRNA_eukarya": "28s_rRNA", "5S_rRNA": "5s_rRNA"}
		rRNAdict = defaultdict(list)
		for line in rrnafile:
			if not line.startswith("#"):
				InfoLINE = re.sub("\s+", ",", line)
				line_splitted = InfoLINE.split(",")
				item_type = line_splitted[0]
				if item_type.startswith(('LSU', 'SSU', '5S')):
					strand = 0
					if line_splitted[9] == "+":
						strand = 1
					else:
						strand = -1
					rRNAdict[item_type].append({'score': float(line_splitted[14]), 'product': namedict[line_splitted[0]], 'begin': int(line_splitted[7]), 'end': int(line_splitted[8]), 'strand': strand})

		subunits = {'LSU': {'max_score': 0 }, 'SSU': {'max_score': 0 }, '5S': {'max_score': 0 }}
		for type_rRNA, rRNA_data in rRNAdict.items():
			max_score = max([item['score'] for item in rRNAdict[type_rRNA]])
			for item in ('LSU', 'SSU'):
				if type_rRNA.startswith(item):
					if max_score > subunits[item]['max_score']:
						subunits[item]['listdata'] = rRNA_data
						subunits[item]['max_score'] = max_score
			if type_rRNA.startswith('5S'):
				subunits['5S']['listdata'] = rRNA_data			
				subunits['5S']['max_score'] = max_score
		
		for rRNA in subunits:
			i = 0
			try:
				lengthlist = len(subunits[rRNA]['listdata'])
			except KeyError:
				continue
			else:
				while i < lengthlist:
					eprint("%s harbours a %s from %i to %i" % (newfile, subunits[rRNA]['listdata'][i]['product'], int(subunits[rRNA]['listdata'][i]['begin']), int(subunits[rRNA]['listdata'][i]['end'])))
					i += 1

	# Predicting the tRNA sequences
	eprint("Running ARAGORN to predict tRNA-like sequences in %s" % newfile)
	genetictable = "-gc%s" % str(args.gcode)
	with open("trnafile.fasta", "w") as trnafile:
		if genomeshape['genomeshape'] == "circular":
			subprocess.call(["aragorn", "-c", "-fon", genetictable, newfile], stdout=trnafile)
		else:
			subprocess.call(["aragorn", "-l", "-fon", genetictable, newfile], stdout=trnafile)
	num_tRNA = len(list(SeqIO.parse("trnafile.fasta", "fasta")))
	eprint("ARAGORN was able to predict %i tRNAs in %s\n" % (num_tRNA, newfile))

	#Storing tRNA and tmRNA information in memory
	with open("trnafile.fasta", "rU") as trnafile:
		tRNAdict = {}
		tmRNAdict = {}
		for tRNAseq in SeqIO.parse(trnafile, "fasta"):
			indtRNA = {}
			indtmRNA = {}
			tRNA_information = tRNAseq.description.split(" ")
			tRNApat = re.compile("^tRNA-")
			if tRNA_information[1] == "tmRNA":
				if tRNA_information[2] == "(Permuted)":
					indtmRNA['product'] = "tmRNA"
					tmRNA_coords = tRNA_information[3]
					Beginningrevcomppat = re.compile("^c")
					if re.match(Beginningrevcomppat, tRNA_coords):
						indtRNA['strand'] = -1
						tmRNA_coords = tmRNA_coords.replace("c[","").replace("]","").split(",")
					else:
						indtRNA['strand'] = 1
						tmRNA_coords = tmRNA_coords.replace("[","").replace("]","").split(",")
					indtmRNA['begin'] = int(tmRNA_coords[0])
					indtmRNA['end'] = int(tmRNA_coords[1])
					tmRNAdict[tRNAseq.id] = indtmRNA
				else:
					indtmRNA['product'] = "tmRNA"
					tmRNA_coords = tRNA_information[2]
					Beginningrevcomppat = re.compile("^c")
					if re.match(Beginningrevcomppat, tRNA_coords):
						indtRNA['strand'] = -1
						tmRNA_coords = tmRNA_coords.replace("c[","").replace("]","").split(",")
					else:
						indtRNA['strand'] = 1
						tmRNA_coords = tmRNA_coords.replace("[","").replace("]","").split(",")
					indtmRNA['begin'] = int(tmRNA_coords[0])
					indtmRNA['end'] = int(tmRNA_coords[1])
					tmRNAdict[tRNAseq.id] = indtmRNA
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
				tRNAdict[tRNAseq.id] = indtRNA

	#Predicting CRISPR repeats and others
	eprint("Running PILERCR to predict CRISPR repeats in %s" % newfile)
	subprocess.call(["pilercr", "-in", newfile, "-out", "crisprfile.txt", "-noinfo", "-minrepeat", str(args.minrepeat), "-maxrepeat", str(args.maxrepeat), "-minspacer", str(args.minspacer), "-maxspacer", str(args.maxspacer)])
	eprint("Predicting repeats in the sequences using TRF and IRF")
	with open("/dev/null", "w") as stderr:
		subprocess.call(["trf", newfile, "2", "7", "7", "80", "10", "50", "500", "-h"], stderr=stderr)
		os.rename("%s.2.7.7.80.10.50.500.dat" % newfile, "trf_temp.dat")
	with open("/dev/null", "w") as stderr:
		subprocess.call(["irf", newfile, "2", "3", "5", "80", "10", "40", "500000", "10000", "-d", "-h"], stderr=stderr)
		os.rename("%s.2.3.5.80.10.40.500000.10000.dat" % newfile, "irf_temp.dat")

	# Storing CRISPR repeats information
	information_CRISPR = {}
	with open("crisprfile.txt", "rU") as crisprfile:
		for line in crisprfile:
			if "SUMMARY BY POSITION" in line:
				for line in crisprfile:
					information_crispr_repeat = {}
					try:
						patC = re.compile('^\s+(\d+)\s+.{16}\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+\d+\s+\d?\s+(\w+)')
						key, start, length, seq = re.match(patC, line).groups()
					except AttributeError:
						continue
					else:
						information_crispr_repeat['start'] = start
						information_crispr_repeat['end'] = int(start) + int(length)
						information_crispr_repeat['repeatseq'] = seq
						information_crispr_repeat['repeatend'] = int(start) + len(seq)
						information_CRISPR[key] = information_crispr_repeat

	# Storing tandem repeats information
	information_TRF = {}
	count = 1
	with open("trf_temp.dat", "rU") as trfile:
		for line in trfile:
			information_tandem_repeat = {}
			try:
				patT = re.compile('^(\d+)\s(\d+)\s\d+\s\d+\.\d+\s')
				start, end = re.match(patT, line).groups()
			except AttributeError:
				continue
			else:
				information_tandem_repeat['start'] = start
				information_tandem_repeat['end'] = end
				information_TRF[count] = information_tandem_repeat
				count += 1

	# Storing inverted repeats information
	information_IRF = {}
	count = 1
	with open("irf_temp.dat", "rU") as irfile:
		for line in irfile:
			information_inverted_repeat = {}
			try:
				patI = re.compile('^(\d+)\s(\d+)\s\d+\s\d+\s\d+')
				start, end = re.match(patI, line).groups()
			except AttributeError:
				continue
			else:
				information_inverted_repeat['start'] = start
				information_inverted_repeat['end'] = end
				information_IRF[count] = information_inverted_repeat
				count += 1

	# Creating a new Genbank and GFF file
	eprint("Creating the output files")
	newtempgbk = "%s.temp.gbk" % newfile
	with open(newfile, "rU") as basefile, open(newtempgbk, "w"):
		for record in SeqIO.parse(basefile, "fasta", IUPAC.ambiguous_dna):
			whole_sequence = SeqRecord(record.seq)
			whole_sequence.id = str(record.id)
			whole_sequence.annotations['data_file_division'] = args.typedata.upper()
			whole_sequence.annotations['date'] = strftime("%d-%b-%Y").upper()
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
			for rRNA in sorted(subunits, key = stringSplitByNumbers):
				i = 0
				try:
					lengthlist = len(subunits[rRNA]['listdata'])
				except KeyError:
					continue
				else:
					while i < lengthlist:
						start_pos = SeqFeature.ExactPosition(subunits[rRNA]['listdata'][i]['begin'])
						end_pos = SeqFeature.ExactPosition(subunits[rRNA]['listdata'][i]['end'])
						feature_location = SeqFeature.FeatureLocation(start_pos, end_pos, strand=subunits[rRNA]['listdata'][i]['strand'])
						new_data_gene = SeqFeature.SeqFeature(feature_location, type = "gene", strand = subunits[rRNA]['listdata'][i]['strand'])
						whole_sequence.features.append(new_data_gene)
						qualifiers = [('product', subunits[rRNA]['listdata'][i]['product'])]
						feature_qualifiers = OrderedDict(qualifiers)
						new_data_rRNA = SeqFeature.SeqFeature(feature_location, type = "rRNA", strand = subunits[rRNA]['listdata'][i]['strand'], qualifiers = feature_qualifiers)
						whole_sequence.features.append(new_data_rRNA)
						i += 1
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
			for tmRNA in sorted(tmRNAdict, key = stringSplitByNumbers):
				start_pos = SeqFeature.ExactPosition(tmRNAdict[tmRNA]['begin'])
				end_pos = SeqFeature.ExactPosition(tmRNAdict[tmRNA]['end'])
				feature_location = SeqFeature.FeatureLocation(start_pos, end_pos, strand=tmRNAdict[tRNA]['strand'])
				new_data_gene = SeqFeature.SeqFeature(feature_location, type = "gene", strand = tmRNAdict[tRNA]['strand'])
				whole_sequence.features.append(new_data_gene)
				qualifiers = [('product', tmRNAdict[tmRNA]['product'])]
				feature_qualifiers = OrderedDict(qualifiers)
				new_data_tmRNA = SeqFeature.SeqFeature(feature_location, type = "tmRNA", strand = tmRNAdict[tmRNA]['strand'], qualifiers = feature_qualifiers)
				whole_sequence.features.append(new_data_tmRNA)
			for CRISPR in sorted(information_CRISPR, key = stringSplitByNumbers):
				start_pos = SeqFeature.ExactPosition(information_CRISPR[CRISPR]['start'])
				end_pos = SeqFeature.ExactPosition(information_CRISPR[CRISPR]['end'])
				feature_location = SeqFeature.FeatureLocation(start_pos, end_pos)
				qualifiers = [('rpt_family', 'CRISPR'), ('rpt_type', 'direct'), ('rpt_unit_range', "%i..%i" % (int(information_CRISPR[CRISPR]['start']), int(information_CRISPR[CRISPR]['repeatend']))), ('rpt_unit_seq', information_CRISPR[CRISPR]['repeatseq'])]
				feature_qualifiers = OrderedDict(qualifiers)
				new_data_CRISPRrepeat = SeqFeature.SeqFeature(feature_location, type = "repeat_region", qualifiers = feature_qualifiers)
				whole_sequence.features.append(new_data_CRISPRrepeat)
			for tandem in sorted(information_TRF):
				start_pos = SeqFeature.ExactPosition(information_TRF[tandem]['start'])
				end_pos = SeqFeature.ExactPosition(information_TRF[tandem]['end'])
				feature_location = SeqFeature.FeatureLocation(start_pos, end_pos)
				qualifiers = [('rpt_type', 'direct')]
				feature_qualifiers = OrderedDict(qualifiers)
				new_data_tandemrepeat = SeqFeature.SeqFeature(feature_location, type = "repeat_region", qualifiers = feature_qualifiers)
				whole_sequence.features.append(new_data_tandemrepeat)
			for inverted in sorted(information_IRF):
				start_pos = SeqFeature.ExactPosition(information_IRF[inverted]['start'])
				end_pos = SeqFeature.ExactPosition(information_IRF[inverted]['end'])
				feature_location = SeqFeature.FeatureLocation(start_pos, end_pos)
				qualifiers = [('rpt_type', 'inverted')]
				feature_qualifiers = OrderedDict(qualifiers)
				new_data_invertedrepeat = SeqFeature.SeqFeature(feature_location, type = "repeat_region", qualifiers = feature_qualifiers)
				whole_sequence.features.append(new_data_invertedrepeat)
			SeqIO.write(whole_sequence, newtempgbk, "genbank")

	newgbk = "%s.gbk" % newfile
	with open(newtempgbk, "rU") as gbktempfile, open(newgbk, "w") as gbkrealfile:
		newpat = re.compile("D|RNA\s+(CON|PHG|VRL|BCT)")
		for line in gbktempfile:
			if line.startswith("LOCUS ") and re.search(newpat, line):
				if genomeshape['genomeshape'] == "linear":
					newline = re.sub("bp    DNA\s+", "bp    DNA     linear   ", line)
				else:
					newline = re.sub("bp    DNA\s+", "bp    DNA     circular ", line)
				gbkrealfile.write(newline)
			else:
				gbkrealfile.write(line)

	for f in glob.glob("*.temp.gbk"):
		os.remove(f)

	if args.gffprint==True:
		newgff = "%s.gff" % newfile
		with open(newgff, "w") as outgff, open(newgbk, "rU") as ingbk:
			GFF.write(SeqIO.parse(ingbk, "genbank"), outgff)

	# Removing intermediate files
	os.remove(newfile)
	os.remove("temporal_circular.fasta")
	os.remove("temp.faa")
	os.remove("temp_blast.csv")
	os.remove("crisprfile.txt")
	os.remove("trnafile.fasta")
	os.remove("rrnafile.csv")
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
			eprint("WARNING: Skipping small contig %s" % rec.id)
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

# Final statement
eprint("Genome annotation done!")
eprint("The GenBank file is %s" % gbkoutputfile)
if args.gffprint==True:
	eprint("The GFF3 file is %s" % gffoutputfile)
eprint("The table file for GenBank submission is %s" % tbloutputfile)
eprint("The FASTA file for GenBank submission is %s" % newfastafile)
eprint("The table file with all protein information is %s" % newtablefile)
