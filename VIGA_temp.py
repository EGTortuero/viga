#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# VIGA - Automatic de-novo VIral Genome Annotator
#
# Copyright (C) 2019 - Enrique Gonzalez-Tortuero
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
import csv
import fractions
import glob
import multiprocessing
import numpy
import os
import os.path
import re
import shutil
import sys
import subprocess
import time
from BCBio import GFF
from Bio import SeqIO, SeqFeature
try:
	from Bio.Alphabet import IUPAC
except ImportError:
	IUPAC = None
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, ExactPosition, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import OrderedDict, defaultdict
from io import StringIO
from itertools import product
from pathlib import Path
import scipy.signal as signal
from time import strftime
from typing import Dict

## Defining the program version
version = "0.11.2"

## Preparing functions
# A batch iterator
def batch_iterator(iterator, batch_size):
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch
       
# Function equivalent to 'cat' command in UNIX systems (used to concatenate all files in a single one)
def cat_all(listinputfiles, outputfile):
    with open(outputfile, 'wb') as outfile:
        for fname in listinputfiles:
            with open(fname, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
                    
# Function equivalent to 'whereis'/'which' commands in UNIX systems (used to find a command in your computer)
def cmd_exists(cmd):
	return subprocess.call(["which", cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

# Function to calculate the GC skew and determine the ori and ter regions
def GCskew(name, length, seq, window, slide):
    replacements = defaultdict(int, {'G': 1, 'C': -1, 'A': 0, 'T': 0, 'N': 0})
    gmc = [replacements[base] for base in seq]
    gpc = [abs(i) for i in gmc]
    weights = numpy.ones(window) / window
    gmc = [[i, c] for i, c in enumerate(signal.fftconvolve(gmc, weights, 'same').tolist())]
    gpc = [[i, c] for i, c in enumerate(signal.fftconvolve(gpc, weights, 'same').tolist())]
    skew = [[], []]
    c_skew = [[], []]
    cs = 0
    for i, m in gmc[0::slide]:
        p = gpc[i][1]
        gcs = m/p if p != 0 else 0
        cs += gcs
        skew[0].append(i)
        c_skew[0].append(i)
        skew[1].append(gcs)
        c_skew[1].append(cs)

    c_skew_min = signal.argrelextrema(numpy.asarray(c_skew[1]), numpy.less, order = 1)[0].tolist()
    c_skew_max = signal.argrelextrema(numpy.asarray(c_skew[1]), numpy.greater, order = 1)[0].tolist()
    if len(c_skew_min) == 0 or len(c_skew_min) == 0:
        return [False, False], skew, c_skew
    else:
        c_skew_min = [[c_skew[0][i], c_skew[1][i]] for i in c_skew_min]
        c_skew_max = [[c_skew[0][i], c_skew[1][i]] for i in c_skew_max]
        closest, farthest = int(length * 0.45), int(length * 0.55)
        pairs = []
        for pair in product(*[[c_skew_min, c_skew_max]]):
            tr, pk = sorted(list(pair), key = lambda x: x[1], reverse = False)
            a = min(tr[0] - pk[0], pk[0] - tr[0])
            pt = abs(tr[1] - pk[1])
            if closest <= a <= farthest:
                pairs.append([pt, tr, pk])
        if not pairs:
            return [False, False], skew, c_skew
        pt, tr, pk = max(pairs, key = lambda x: x[0])
        return [tr[0], pk[0]], skew, c_skew

# Function to print all messages in standard error instead of standard output
def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Function to split a string numerically and not by alphabetic order
def stringSplitByNumbers(x):
	r = re.compile('(\d+)')
	l = r.split(x)
	return [int(y) if y.isdigit() else y for y in l]

# Functions to process the arguments
def create_parser():
    parser = argparse.ArgumentParser(description='VIGA is an automatic de novo VIral Genome Annotator.')
    return parser

def add_basic_group(parser):
    basic_group = parser.add_argument_group('Basic options for VIGA [REQUIRED]')
    basic_group.add_argument("--input", dest="inputfile", type=str, required=True, help='Input file as a FASTA file', metavar="FASTAFILE")
    basic_group.add_argument("--modifiers", dest="modifiers", type=str, required=True, help='Input file as a plain text file with the modifiers per every FASTA header according to SeqIn (https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html). All modifiers must be written in a single line and are separated by a single space character. No space should be placed besides the = sign. For example: [organism=Serratia marcescens subsp. marcescens] [sub-species=marcescens] [strain=AH0650_Sm1] [topology=linear] [moltype=DNA] [tech=wgs] [gcode=11] [country=Australia] [isolation-source=sputum]. This line will be copied and printed along with the record name as the definition line of every contig sequence.', metavar="TEXTFILE")
    return parser

def add_advanced_general_group(parser):
    advanced_General_group = parser.add_argument_group('Advanced general options for VIGA [OPTIONAL]')
    advanced_General_group.add_argument("--out", dest="root_output", type=str, help='Name of the outputs files (without extension)', metavar="OUTPUTNAME")
    advanced_General_group.add_argument("--locus", dest="locus", type=str, default='LOC', help='Name of the sequences. If the input is a multiFASTA file, please put a general name as the program will add the number of the contig at the end of the name (Default: %(default)s)', metavar="STRING")
    advanced_General_group.add_argument("--threads", dest="ncpus", default=multiprocessing.cpu_count(), help='Number of threads/cpus (Default: %(default)s cpu [Max. number of cpus])', metavar="INT")
    advanced_General_group.add_argument('--mincontigsize', dest="mincontigsize", type=int, default = 200, help = 'Minimum contig length to be considered in the final files (Default: 200 bp)', metavar="INT")
    advanced_General_group.add_argument("--blast", dest="blastswitch", action='store_true', default=False, help='Using BLAST to predict protein function based on homology. Alternatively, DIAMOND will be used for such purpose (Default: False)')
    return parser

def add_advanced_circularity_group(parser):
    advanced_circularity_group = parser.add_argument_group('Advanced options for contig shape prediction in VIGA [OPTIONAL]')
    advanced_circularity_group.add_argument("--readlength", dest="read_length", type=int, default=101, help='Read length for the circularity prediction (default: 101 bp)', metavar="INT")
    advanced_circularity_group.add_argument("--windowsize", dest="window", type=int, default=100, help='Window size used to determine the origin of replication in circular contigs according to the cumulative GC skew (default: 100 bp)', metavar="INT")
    advanced_circularity_group.add_argument("--slidingsize", dest="slide", type=int, default=10, help='Window size used to determine the origin of replication in circular contigs according to the cumulative GC skew (default: 10 bp)', metavar="INT")
    return advanced_circularity_group

def add_advanced_repeats_group(parser):
    advanced_repeats_group = parser.add_argument_group('Advanced options for detection of repeats in VIGA [OPTIONAL]')
    advanced_repeats_group.add_argument("--minrepeat", dest="minrepeat", type=int, default=16, help="Minimum repeat length for CRISPR detection (Default: 16)", metavar="INT")
    advanced_repeats_group.add_argument("--maxrepeat", dest="maxrepeat", type=int, default=64, help="Maximum repeat length for CRISPR detection (Default: 64)")
    advanced_repeats_group.add_argument("--minspacer", dest="minspacer", type=int, default=8, help="Minimum spacer length for CRISPR detection (Default: 8)")
    advanced_repeats_group.add_argument("--maxspacer", dest="maxspacer", type=int, default=64, help="Maximum spacer length for CRISPR detection (Default: 64)")
    return advanced_repeats_group

def add_gcode_group(parser):
    advanced_gcode_group = parser.add_argument_group('Advanced options for genetic code in VIGA [OPTIONAL]')
    gcode_choices = {'1': 'Standard genetic code [Eukaryotic]', '2': 'Vertebrate mitochondrial code', '3': 'Yeast mitochondrial code', '4': 'Mycoplasma/Spiroplasma and Protozoan/mold/coelenterate mitochondrial code', '5': 'Invertebrate mitochondrial code', '6': 'Ciliate, dasycladacean and hexamita nuclear code', '9': 'Echinoderm/flatworm mitochondrial code', '10': 'Euplotid nuclear code', '11': 'Bacteria/Archaea/Phages/Plant plastid', '12': 'Alternative yeast nuclear code', '13': 'Ascidian mitochondrial code', '14': 'Alternative flatworm mitochondrial code', '16': 'Chlorophycean mitochondrial code', '21': 'Trematode mitochondrial code', '22': 'Scedenesmus obliquus mitochondrial code', '23': 'Thraustochytrium mitochondrial code', '24': 'Pterobranquia mitochondrial code', '25': 'Gracilibacteria & Candidate division SR1', '26': 'Pachysolen tannophilus nuclear code', '27': 'Karyorelict nuclear code', '28': 'Condylostoma nuclear code', '29': 'Mesodinium nuclear code', '30': 'Peritrich nuclear code', '31': 'Blastocrithidia nuclear code', '33': 'Cephalodiscidae mitochondrial code'}
    gcode_help = ('Number of GenBank translation table. At this moment, the available options are {0}. (Default: %(default)s)'.format('{}'.format(', '.join('{0} ({1})'.format(k, v) for k, v in sorted(gcode_choices.items())))))
    advanced_gcode_group.add_argument("--gcode", dest="gcode", type=str, default='11', help=gcode_help, metavar="NUMBER")

def add_type_group(parser):
    advanced_type_group = parser.add_argument_group('Advanced options for GenBank division in VIGA [OPTIONAL]')
    type_choices = {'BCT': 'Prokaryotic chromosome', 'CON': 'Contig', 'PHG': 'Phages', 'VRL': 'Eukaryotic/Archaea virus'}
    type_help = ('GenBank Division: One of the following codes - {0}. (Default: %(default)s)'.format(', '.join('{0} ({1})'.format(k, v) for k, v in type_choices.items())))
    advanced_type_group.add_argument("--typedata", dest="typedata", type=str, default='CON', help=type_help, metavar="BCT|CON|VRL|PHG")

def add_rfam_group(parser):
    advanced_rfam_group = parser.add_argument_group('Advanced options for the ncRNA prediction based on Covariance Models in VIGA [OPTIONAL]')
    advanced_rfam_group.add_argument("--norfam", dest="norfam", action='store_true', default=False, help="Don't run RFAM to predict other ncRNAs, apart of rRNAs and tRNAs. (Default: False)")
    advanced_rfam_group.add_argument("--rfamdb", dest="rfamdatabase", type=str, help='RFAM Database that will be used for the ncRNA prediction. RFAMDB should be in the format "/full/path/to/rfamdb/Rfam.cm" and must be compressed accordingly (see INFERNAL manual) before running the script. By default, the program will try to search Rfam inside the folder database/ (after running the Create_databases.sh script)', metavar="RFAMDB")
    return parser

def add_diamond_group(parser):
    advanced_diamond_group = parser.add_argument_group('Advanced options for the first protein function prediction based on homology in VIGA [OPTIONAL]')
    advanced_diamond_group.add_argument("--diamonddb", dest="diamonddatabase", type=str, help='DIAMOND Database that will be used for the protein function prediction. The database must be created from a amino acid FASTA file as indicated in https://github.com/bbuchfink/diamond. By default, the program will try to search the RefSeq Viral Protein DB formatted for its use in Diamond inside the folder database/ only if --blast parameter is disabled and after running the Create_databases.sh script', metavar="DIAMONDDB")
    advanced_diamond_group.add_argument("--diamondevalue", dest="diamondevalue", default=0.00001, help='DIAMOND e-value threshold (Default: 0.00001)', metavar="FLOAT")
    advanced_diamond_group.add_argument("--diamondwidthr", dest="diamondwidthreshold", default=50.00, help='DIAMOND ID threshold (Default: 50.0)', metavar="FLOAT")
    advanced_diamond_group.add_argument("--diamondcoverthr", dest="diamondcovthreshold", default=50.00, help='DIAMOND Coverage threshold (Default: 50.0)', metavar="FLOAT")
    return parser
    
def add_advanced_blast_group(parser):
    advanced_blast_group = parser.add_argument_group('Advanced options for the second protein function prediction based on homology in VIGA [OPTIONAL]')
    advanced_blast_group.add_argument("--blastdb", dest="blastdatabase", type=str, help='BLAST Database that will be used to refine the protein function prediction in hypothetical proteins. The database must be an amino acid one, not nucleotidic. By default, the program will try to search the RefSeq Viral Protein DB formatted for its use in BLAST inside the folder database/ only if --blast parameter is ensabled and after running the Create_databases.sh script', metavar="BLASTDB")
    advanced_blast_group.add_argument("--blastexh", dest="blastexh", action='store_true', default=False, help='Use of exhaustive BLAST to predict the proteins by homology according to Fozo et al. (2010) Nucleic Acids Res (Default=FALSE)')
    advanced_blast_group.add_argument("--blastevalue", dest="blastevalue", default=0.00001, help='BLAST e-value threshold (Default: 0.00001)', metavar="FLOAT")
    advanced_blast_group.add_argument("--blastwidthr", dest="blastwidthreshold", default=50.00, help='BLAST ID threshold (Default: 50.0)', metavar="FLOAT")
    advanced_blast_group.add_argument("--blastcoverthr", dest="blastcovthreshold", default=50.00, help='BLAST Coverage threshold (Default: 50.0)', metavar="FLOAT")

def add_advanced_hmm_group(parser):
    advanced_hmm_group = parser.add_argument_group('Advanced options for the third protein function prediction based on HMM in VIGA [OPTIONAL]')
    advanced_hmm_group.add_argument("--nohmmer", dest="nohmmer", action='store_true', default=False, help='Running only BLAST to predict protein function. (Default: False)')
    advanced_hmm_group.add_argument("--novogs", dest="novogs", action='store_true', default=False, help='HMMer analyses will not consider the VOGs database. (Default: False)')
    advanced_hmm_group.add_argument("--norvdb", dest="norvdb", action='store_true', default=False, help='HMMer analyses will not consider the RVDB database. (Default: False)')
    advanced_hmm_group.add_argument("--nophrogs", dest="nophrogs", action='store_true', default=False, help='HMMer analyses will not consider the PHROGs database (Default: False)')
    advanced_hmm_group.add_argument("--vogsdb", dest="vogsdb", type=str, help='VOG Database that will be used to refine the protein function prediction in hypothetical proteins. By default, the program will try to search the Viral Orthologous Genes DB formatted for its use in HMMer inside the folder database/ only if --nohmmer parameter is disabled and after running the Create_databases.sh script.')
    advanced_hmm_group.add_argument("--rvdbdb", dest="rvdbdb", type=str, help='RVDB Database that will be used to refine the protein function prediction in hypothetical proteins. By default, the program will try to search the Reference Virus DataBase formatted for its use in HMMer inside the folder database/ only if --nohmmer parameter is disabled and after running the Create_databases.sh script.')
    advanced_hmm_group.add_argument("--phrogsdb", dest="phrogsdb", type=str, help='PHROGs Database that will be used to refine the protein function prediction in hypothetical proteins. By default, the program will try to search the Prokaryotic Virus Remote Homologous Groups formatted for its use in HMMer inside the folder database/ only if --nohmmer parameter is disabled and after running the Create_databases.sh script.')
    advanced_hmm_group.add_argument("--hmmerevalue", dest="hmmerevalue", default=0.001, help='HMMER e-value threshold (Default: 0.001)', metavar="FLOAT")
    advanced_hmm_group.add_argument("--hmmeridthr", dest="hmmeridthreshold", default=50.00, help='HMMER ID threshold (Default: 50.0)', metavar="FLOAT")
    advanced_hmm_group.add_argument("--hmmercoverthr", dest="hmmercovthreshold", default=50.00, help='HMMER Coverage threshold (Default: 50.0)', metavar="FLOAT")

# Interpreting the parameters
def check_and_set_default_databases(args):
    root_output = args.root_output
    if not root_output:
        root_output = '{}_annotated'.format(os.path.splitext(args.inputfile)[0])

    if args.norfam == False and args.rfamdatabase == None:
        my_file = Path(os.path.dirname(os.path.abspath(__file__)) + '/databases/rfam/Rfam.cm')
        try:
            my_abs_path = my_file.resolve(strict=True)
        except FileNotFoundError:
            sys.exit('You MUST specify RFAM database using the parameter --rfamdb if you are not using --norfam option')
        else:
            args.rfamdatabase = my_file

    if args.blastswitch == False and args.diamonddatabase == None:
        my_file = Path(os.path.dirname(os.path.abspath(__file__)) + '/databases/RefSeq_Viral_DIAMOND/refseq_viral_proteins.dmnd')
        try:
            my_abs_path = my_file.resolve(strict=True)
        except FileNotFoundError:
            sys.exit('You MUST specify RefSeq Viral Protein Database formatted for Diamond using the parameter --diamonddb if you are not using --blast option')
        else:
            args.diamonddatabase = my_file

    if args.blastswitch == True and args.blastdatabase == None:
        my_file = Path(os.path.dirname(os.path.abspath(__file__)) + '/databases/RefSeq_Viral_BLAST/refseq_viral_proteins.pdb')
        try:
            my_abs_path = my_file.resolve(strict=True)
        except FileNotFoundError:
            sys.exit('You MUST specify RefSeq Viral Protein Database formatted for BLAST using the parameter --blastdb if you are using the --blast option')
        else:
            args.blastdatabase = os.path.dirname(os.path.abspath(__file__)) + '/databases/RefSeq_Viral_BLAST/refseq_viral_proteins'

    if args.nohmmer == False:
        if args.novogs == True and args.norvdb == True and args.nophrogs == True:
            sys.exit('You MUST use one HMMer database, at least, to run this step! ')
        if args.novogs == False and args.vogsdb == None:
            vogs_path = Path(os.path.dirname(os.path.abspath(__file__)) + '/databases/vogs/vog_latest.hmm')
            try:
                my_abs_path = vogs_path.resolve(strict=True)
            except FileNotFoundError:
                sys.exit('You do not have installed VOGS database')
            else:
                args.vogsdb = vogs_path
        if args.norvdb == False and args.rvdbdb == None:
            rvdb_path = Path(os.path.dirname(os.path.abspath(__file__)) + '/databases/rvdb/RVDB_26.0.hmm')
            try:
                my_abs_path = rvdb_path.resolve(strict=True)
            except FileNotFoundError:
                sys.exit('You do not have installed RVDB database')
            else:
                args.rvdbdb = rvdb_path
        if args.nophrogs == False and args.phrogsdb == None:
            phrogs_path = Path(os.path.dirname(os.path.abspath(__file__)) + '/databases/phrogs/phrogs_v4.hmm')
            try:
                args.phrogsdb = phrogs_path.resolve(strict=True)
            except FileNotFoundError:
                sys.exit('You do not have installed PHROGS database')
            else:
                args.phrogsdb = phrogs_path

    return root_output

def check_required_cmds(args):
    required_cmds = ["lastz", "aragorn", "pilercr", "prodigal", "diamond"]

    if args.norfam == False:
        required_cmds.append("cmscan")

    if args.blastswitch == True:
        required_cmds.append("blastp")

    if args.nohmmer == False:
        required_cmds.extend(["hmmsearch", "hmmbuild"])

    missing_cmds = [cmd for cmd in required_cmds if not cmd_exists(cmd)]

    if missing_cmds:
        sys.exit(f"The following commands are missing: {', '.join(missing_cmds)}. You need to run the installer.sh script before running this pipeline")

## Function to correct the original file (long headers problem + multiple FASTA files)
def rename_sequences(args, record_iter):
    newnamessequences = {}
    counter = 1
    for i, batch in enumerate(batch_iterator(record_iter, 1)):
        seq_index = "LOC_%i" % (i + 1)
        with open("%s.temp.fna" % seq_index, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        with open("%s.temp.fna" % seq_index, "r") as original, open("%s.fna" % seq_index, "w") as corrected:
            if IUPAC:
                sequences = SeqIO.parse(original, "fasta", IUPAC.ambiguous_dna)
            else:
                sequences = SeqIO.parse(original, "fasta")
            for record in sequences:
                original_name = record.id
                record.id = "%s_%i" % (args.locus, counter)
                record.description = record.description
                counter += 1
                newnamessequences[record.id] = original_name
                print("WARNING: %s was renamed as %s" % (original_name, record.id))
            SeqIO.write(record, corrected, "fasta")
        with open("logfile.txt", "w") as logfile:
            logfile.write("#Original\tNew\n")
            for oldname in sorted(newnamessequences, key = stringSplitByNumbers):
                logfile.write("%s\t%s\n" % (newnamessequences[oldname], oldname))
        os.remove("%s.temp.fna" % seq_index)
    return newnamessequences

# Determination of the circular shape contigs based on Alex Crits-Christoph's find_circular.py script [https://github.com/alexcritschristoph/VRCA/blob/master/find_circular.py])
def process_resultlastz(resultlastz, record, genomeshape, args):
    if resultlastz != '':
        start1, end1, start2, end2, _, strand1, strand2, identity, _, length = resultlastz.split()
        length = int(length)
        identity = float(fractions.Fraction(identity))
        if (strand1 == strand2 and length > 0.4 * args.read_length and identity > 0.95 and int(start1) < 5 and int(start2) > args.read_length and int(end1) < args.read_length and int(end2) > args.read_length * 2 * 0.9):
            genomeshape[record.id]["genomeshape"] = "circular"
            genomeshape[record.id]["identity"] = min(identity, genomeshape[record.id].get("identity", identity))
            genomeshape[record.id]["length"] = max(length, genomeshape[record.id].get("length", length))

def determine_genome_shape(record, genomeshape, args):
    try:
        if genomeshape[record.id]["genomeshape"] == "":
            genomeshape[record.id]["genomeshape"] = "linear"
    except KeyError:
        genomeshape[record.id]["genomeshape"] = "linear"

    if genomeshape[record.id]["genomeshape"] == "circular":
        with open("temp.fasta", "w") as correctedcircular:
            corr_seq = str(record.seq[: int(genomeshape[record.id]["length"]) // 2:-int(genomeshape[record.id]["length"]) // 2])
            if IUPAC:
                new_seq = SeqRecord(Seq(corr_seq, IUPAC.ambiguous_dna), id=record.description)
            else:
                new_seq = SeqRecord(Seq(corr_seq), id=record.description)
            SeqIO.write(new_seq, correctedcircular, "fasta")
        os.rename("temp.fasta", "temp2.fasta")

## Calculate the cumulative GCskew in circular contigs to determine the origin of replication (Based on iRep -- Brown CT, Olm MR, Thomas BC, Banfield JF (2016) Measurement of bacterial replication rates in microbial communities. Nature Biotechnology 34: 1256-63.)
def circularize_sequence(newfile, genomeshape, record, args, IUPAC):
    if genomeshape[record.id]['genomeshape'] != "circular":
        return
    eprint("Determining the origin of replication of %s according to the GC Skew" % newfile)
    for gotocircularize in SeqIO.parse("temp2.fasta", "fasta"):
        oric, term, _, _ = GCskew(gotocircularize.id, len(gotocircularize.seq), gotocircularize.seq, args.window, args.slide)
        if oric is False:
            oric, term = 'n/a', 'n/a'
        else:
            oric, term = '{:,}'.format(oric), '{:,}'.format(term)
        eprint('%s -> Origin: %s Terminus: %s' % (gotocircularize.id, oric, term))
        if oric is not False:
            firstpartseq = str(gotocircularize.seq[oric:-1])
            secondpartseq = str(gotocircularize.seq[0:(oric-1)])
            combinedcorrectedseq = SeqRecord(Seq(firstpartseq + secondpartseq, IUPAC.ambiguous_dna) if IUPAC else Seq(firstpartseq + secondpartseq), id=gotocircularize.description)
            SeqIO.write(combinedcorrectedseq, newfile, "fasta")
        else:
            eprint("VIGA was unable to predict the origin of replication: %s was not modified!" % record.id)
            os.rename("temp2.fasta", newfile)
    os.remove("temp2.fasta") if os.path.isfile("temp2.fasta") else None

def get_sequence(newfile):
    with open(newfile, "r") as targetfile:
        seq_string = targetfile.read()
        return SeqIO.parse(StringIO(seq_string), "fasta")

def get_combined_seqs(record, args, IUPAC):
    seq_beginning = str(record.seq[0:args.read_length])
    seq_ending = str(record.seq[len(record.seq)-args.read_length:len(record.seq)])
    if IUPAC:
        return SeqRecord(Seq(seq_beginning + seq_ending, IUPAC.ambiguous_dna), id = record.description)
    else:
        return SeqRecord(Seq(seq_beginning + seq_ending), id = record.description)

def write_temp_file(combined_seqs):
    SeqIO.write(combined_seqs, "temporal_circular.fasta", "fasta")

def run_lastz(combined_seqs):
    return subprocess.check_output(["lastz", "temporal_circular.fasta", "--self", "--notrivial", "--nomirror", "--ambiguous=iupac", "--format=general-:start1,end1,start2,end2,score,strand1,strand2,identity,length1"])

def parse_lastz_output(outputlastz):
    return outputlastz.decode().split("\n")

def process_resultlastz(resultlastz, record, genomeshape, args):
    if isinstance(resultlastz, str) and resultlastz != '':
        start1, end1, start2, end2, _, strand1, strand2, identity, _, length = resultlastz.split()
        length = int(length)
        identity = float(fractions.Fraction(identity))
        if (strand1 == strand2 and length > 0.4 * args.read_length and identity > 0.95 and int(start1) < 5 and int(start2) > args.read_length and int(end1) < args.read_length and int(end2) > args.read_length * 2 * 0.9):
            genomeshape.setdefault(record.id, {})  # Initialize the dictionary if it doesn't exist
            genomeshape[record.id]["genomeshape"] = "circular"
            genomeshape[record.id]["identity"] = min(identity, genomeshape[record.id].get("identity", identity))
            genomeshape[record.id]["length"] = max(length, genomeshape[record.id].get("length", length))
        else:
            genomeshape.setdefault(record.id, {})  # Initialize the dictionary if it doesn't exist
            genomeshape[record.id]["genomeshape"] = "linear"
            genomeshape[record.id]["identity"] = identity
            genomeshape[record.id]["length"] = length
            
def determine_genome_shape(record, genomeshape, args):
    genomeshape.setdefault(record.id, {})  # Initialize the dictionary if it doesn't exist
    if "genomeshape" not in genomeshape[record.id]:
        genomeshape[record.id]["genomeshape"] = "linear"

    if genomeshape[record.id]["genomeshape"] == "circular":
        with open("temp.fasta", "w") as correctedcircular:
            corr_seq = str(record.seq[: int(genomeshape[record.id]["length"]) // 2:-int(genomeshape[record.id]["length"]) // 2])
            if IUPAC:
                new_seq = SeqRecord(Seq(corr_seq, IUPAC.ambiguous_dna), id=record.description)
            else:
                new_seq = SeqRecord(Seq(corr_seq), id=record.description)
            SeqIO.write(new_seq, correctedcircular, "fasta")
        os.rename("temp.fasta", "temp2.fasta")

def log_genome_shape(record, genomeshape, logfile):
    eprint("%s seems to be a %s contig according to LASTZ" % (record.id, genomeshape[record.id]['genomeshape']))
    logfile.write("%s\t%s\n" % (record.id, genomeshape[record.id]['genomeshape']))

def predict_shape(newfile, args, IUPAC, logfile):
    Sequence = get_sequence(newfile)
    for record in Sequence:
        genomeshape[record.id] = {}
        combined_seqs = get_combined_seqs(record, args, IUPAC)
        write_temp_file(combined_seqs)
        outputlastz = run_lastz(combined_seqs)
        resultslastz = parse_lastz_output(outputlastz)
        process_resultlastz(resultslastz, record, genomeshape, args)
        determine_genome_shape(record, genomeshape, args)
        log_genome_shape(record, genomeshape, logfile)
    return Sequence, genomeshape

def process_trna_file(newrecord, putativetrnafile, tRNAdict, tmRNAdict):
    with open(putativetrnafile, "r") as trnafile:
        tRNAdict[newrecord.id] = {}
        tmRNAdict[newrecord.id] = {}
        for tRNAseq in SeqIO.parse(trnafile, "fasta"):
            tRNA_information = tRNAseq.description.split(" ")
            if tRNA_information[1] == "tmRNA":
                indtmRNA = process_tmRNA(tRNA_information)
                tmRNAdict[newrecord.id][tRNAseq.id] = indtmRNA
            elif re.match("^tRNA-", tRNA_information[1]):
                indtRNA = process_tRNA(tRNA_information)
                tRNAdict[newrecord.id][tRNAseq.id] = indtRNA

def process_tmRNA(tRNA_information):
    indtmRNA = {}
    indtmRNA['product'] = "tmRNA"
    tmRNA_coords = str(tRNA_information[2] if tRNA_information[2] != "(Permuted)" else tRNA_information[3])
    indtmRNA['strand'] = -1 if re.match("^c", tmRNA_coords) else 1
    tmRNA_coords = tmRNA_coords.replace("c[","").replace("[","").replace("]","").split(",")
    indtmRNA['begin'] = int(tmRNA_coords[0])
    indtmRNA['end'] = int(tmRNA_coords[1])
    return indtmRNA

def process_tRNA(tRNA_information):
    indtRNA = {}
    indtRNA['product'] = re.sub("\(\w{3}\)", "", tRNA_information[1])
    tRNA_coords = tRNA_information[2]
    indtRNA['strand'] = -1 if re.match("^c", tRNA_coords) else 1
    tRNA_coords = tRNA_coords.replace("c[","").replace("[","").replace("]","").split(",")
    indtRNA['begin'] = int(tRNA_coords[0])
    indtRNA['end'] = int(tRNA_coords[1])
    return indtRNA

def process_contig_for_tRNAs(contigfile, tRNAdict, tmRNAdict, genetictable):
    with open(contigfile, "r") as targetfasta:
        for newrecord in SeqIO.parse(targetfasta, "fasta"):
            putativetrnafile = "trnafile_%s.fasta" % newrecord.id
            with open(putativetrnafile, "w") as trnafile:
                if genomeshape[newrecord.id]['genomeshape'] == "circular":
                    subprocess.call(["aragorn", "-c", "-fon", genetictable, contigfile], stdout=trnafile)
                else:
                    subprocess.call(["aragorn", "-l", "-fon", genetictable, contigfile], stdout=trnafile)
            num_tRNA = len(list(SeqIO.parse(putativetrnafile, "fasta")))
            process_trna_file(newrecord, putativetrnafile, tRNAdict, tmRNAdict)
    return num_tRNA

def run_cmscan(rfamdatabase, ncpus, fasta_file, output_file):
    with open("/dev/null", "w") as stderr:
        subprocess.call(["cmscan", "--rfam", "--cut_ga", "--nohmmonly", "--tblout", output_file, "--cpu", str(ncpus), rfamdatabase, fasta_file], stdout=stderr)
    return output_file

def parse_and_read_ncrnafile(filename):
    elementsncRNA = {}
    with open(filename, "r") as ncrnafile:
        count = 0
        for line in ncrnafile:
            if not line.startswith("#"):
                InfoLINE = re.sub("\s{2,}", ",", line)
                line_splitted = InfoLINE.split(",")
                item_type = line_splitted[0]
                contig_id = line_splitted[2]
                if contig_id not in elementsncRNA:
                    elementsncRNA[contig_id] = {'Other': {}}
                if item_type.startswith(('LSU', 'SSU', '5S', '5_8S', 'tRNA')):
                    next
                elif item_type.endswith('tmRNA'):
                    next
                else:
                    count += 1
                    if len(line_splitted) == 16:
                        if line_splitted[9] == "+":
                            elementsncRNA[contig_id]['Other'][count] = {'type': item_type, 'product': line_splitted[15].replace("\n", ""), 'begin': int(line_splitted[7]), 'end': int(line_splitted[8]), 'score': float(line_splitted[14].replace(" !", "")), 'rfamcode': line_splitted[1], 'strand': 1}
                        else:
                            elementsncRNA[contig_id]['Other'][count] = {'type': item_type, 'product': line_splitted[15].replace("\n", ""), 'begin': int(line_splitted[8]), 'end': int(line_splitted[7]), 'score': float(line_splitted[14].replace(" !", "")), 'rfamcode': line_splitted[1], 'strand': -1}
                    elif len(line_splitted) == 15:
                        if line_splitted[9] == "+":
                            elementsncRNA[contig_id]['Other'][count] = {'type': item_type, 'product': line_splitted[14].replace("\n", ""), 'begin': int(line_splitted[7]), 'end': int(line_splitted[8]), 'score': float(line_splitted[13].replace(" !", "")), 'rfamcode': line_splitted[1], 'strand': 1}
                        else:
                            elementsncRNA[contig_id]['Other'][count] = {'type': item_type, 'product': line_splitted[14].replace("\n", ""), 'begin': int(line_splitted[8]), 'end': int(line_splitted[7]), 'score': float(line_splitted[13].replace(" !", "")), 'rfamcode': line_splitted[1], 'strand': -1}
    return elementsncRNA

def extract_crispr_data(filename):
    information_CRISPR = {}
    with open(filename, "r") as crisprfile:
        for line in crisprfile:
            if "SUMMARY BY POSITION" in line:
                for line in crisprfile:
                    try:
                        patC1 = re.compile('^\s+\d+\s+(.{16})\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d?\s+\w+')
                        sequence_id = re.split("\s", re.match(patC1, line).groups()[0])[0]
                        information_CRISPR[sequence_id] = {}
                    except AttributeError:
                        continue
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
    return information_CRISPR

def predict_genes(contigfile, record, orffile, orffile2, length_contig, genomeshape, args):
    cmd = ["prodigal", "-a", "pretemp.faa", "-d", orffile2, "-i", contigfile, "-o", "/dev/null", "-g", args.gcode, "-q"]
    if length_contig >= 100000:
        if genomeshape[record.id]['genomeshape'] == 'linear':
            cmd += ["-c"]
    else:
        cmd += ["-p", "meta"]
        if genomeshape[record.id]['genomeshape'] == 'linear':
            cmd += ["-c"]
    subprocess.call(cmd)
    with open("pretemp.faa", "r") as originalfaa, open(orffile, "w") as correctedfaa:
        sequences = SeqIO.parse(originalfaa, "fasta")
        for record in sequences:
            record.seq = record.seq.rstrip("*")
            SeqIO.write(record, correctedfaa, "fasta")
    num_seqs = len(list(SeqIO.parse("pretemp.faa", "fasta")))
    eprint("Detected %i ORFs in %s" % (num_seqs, contigfile))

def run_blast(file_name, blast_database, blast_evalue, ncpus, blastexh=False):
    with open(file_name, "r") as input_file:
        first_line = input_file.readline()
        if first_line.startswith(">"):
            if blastexh:
                print("\nRunning BLAST to predict the genes according to homology inference using exhaustive mode (see Fozo et al. (2010) Nucleic Acids Res for details)")
                subprocess.call(['blastp', '-query', file_name, '-db', blast_database, '-evalue', str(blast_evalue), '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle', '-out', file_name + ".csv", '-word_size', '2', '-gapopen', '8', '-gapextend', '2', '-matrix', 'PAM70', '-comp_based_stats', '"0"', "-num_threads", str(ncpus)])
            else:
                print("\nRunning BLAST to predict the genes according to homology inference using default parameters")
                subprocess.call(['blastp', '-query', file_name, '-db', blast_database, '-evalue', str(blast_evalue), '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle', '-out', file_name + ".csv", "-num_threads", str(ncpus)])
        else:
            open(file_name + ".csv", 'a').close()

def run_diamond_blastp(inputfile, outputfile, diamonddatabase, diamondevalue, ncpus):
    with open(inputfile, "r") as inputstep:
        first_line = inputstep.readline()
        if first_line.startswith(">"):
            subprocess.call(['diamond', 'blastp', '-q', inputfile, '-d', diamonddatabase, '-e', str(diamondevalue), '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qstart', 'qend', 'evalue', 'bitscore', 'stitle', '-o', outputfile, "-p", str(ncpus), '--quiet'])
        else:
            open(outputfile, 'a').close()

def parse_blast_diamond_results(file_path, threshold_width, threshold_cov, threshold_evalue, program):
    hypotheticalpat = re.compile(r'(?i)(((hypothetical|uncharacteri[z|s]ed|predicted)( phage)?( membrane)? protein)|(ORF|(unnamed protein product|gp\d+|protein of unknown function|phage protein)))')
    information_proteins = {}
    with open(file_path, "r") as results:
        reader = csv.DictReader(results, delimiter='\t', fieldnames=['qseqid','sseqid','pident','length','qlen','slen','qstart','qend','evalue','bitscore','stitle'])
        for row in reader:
            perc_cover = round(100.00 * (float(row['length']) / float(row['qlen'])), 2)
            perc_id = float(row['pident'])
            infoprot = {}
            infoprot['sseqid'] = row['sseqid']
            infoprot['pident'] = perc_id
            infoprot['pcover'] = perc_cover
            infoprot['evalue'] = row['evalue']
            infoprot['descr'] = row['stitle']
            try:
                data = information_proteins[row['qseqid']]
            except KeyError:
                if not re.search(hypotheticalpat, infoprot['descr']) and perc_id >= threshold_width and perc_cover >= threshold_cov and float(row['evalue']) <= threshold_evalue:
                    information_proteins[row['qseqid']] = infoprot
                else:
                    continue
            else:
                if program == "blast":
                    if not re.search(hypotheticalpat, infoprot['descr']) and perc_id >= threshold_width and perc_id >= data['pident'] and perc_cover >= threshold_cov and perc_cover >= data['pcover'] and float(row['evalue']) <= threshold_evalue:
                        information_proteins[row['qseqid']] = infoprot
                elif program == "diamond":
                    if not re.search(hypotheticalpat, infoprot['descr']) and perc_id >= threshold_width and perc_id >= float(data['pident']) and perc_cover >= threshold_cov and perc_cover >= float(data['pcover']) and float(row['evalue']) <= threshold_evalue:
                        information_proteins[row['qseqid']] = infoprot
    return information_proteins

def parse_hmmer_file(tbl_file, hmmercovthreshold, hmmerevalue, hmmeridthreshold, db_name):
    information_proteins_hmmer = {}
    with open(tbl_file, "r") as tblfile:
        for line in tblfile:
            if not line.startswith("#"):
                infoprot_hmmer = {}
                match = re.match("^(\S+)\s+\S\s+\d+\s+(\S+)\s+\S\s+(\d+)\s+((?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?)\s+\S+\s+\S+\s+\S+\s+\S+\s+(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+\S+\s+\S+\s+(\S+)\s+(.*)", line)
                if match:
                    matchname, lociname, length, evaluehh, start, end, pident, description = match.groups()
                    length = float(length)
                    pident = 100.00 * float(pident)
                    protarea = float(end) - float(start)
                    if protarea < 0:
                        protarea = protarea + length
                    covprot = (protarea / length) * 100
                    if pident >= hmmeridthreshold and covprot >= hmmercovthreshold and float(evaluehh) <= hmmerevalue:
                        infoprot_hmmer['db_name'] = db_name
                        infoprot_hmmer['idhmm'] = lociname
                        infoprot_hmmer['name'] = matchname
                        infoprot_hmmer['pident'] = pident
                        infoprot_hmmer['pcover'] = covprot
                        infoprot_hmmer['evalue'] = evaluehh
                        infoprot_hmmer['description'] = description
                        if matchname not in information_proteins_hmmer:
                            information_proteins_hmmer[matchname] = {}
                        information_proteins_hmmer[matchname][db_name] = infoprot_hmmer
    return {k: v for k, v in information_proteins_hmmer.items() if len(v) == 1}

#def translate_PHROGS_idhmm(idhmm):
#    # Read the PHROGS table into a dictionary mapping IDHMMs to descriptions
#    idhmm_to_desc = {}
#    phrogstable = Path(os.path.dirname(os.path.abspath(__file__)) + '/databases/phrogs/phrog_annot_v4.tsv')
#    with open(phrogstable) as f:
#        for i, line in enumerate(f):
#            if i == 0:
#                continue # skip first line
#            fields = line.strip().split('\t')
#            idhmm = int(fields[0])
#            desc = fields[3] if fields[2] == 'NA' else fields[2]
#            idhmm_to_desc[idhmm] = desc
#    
#    # Look up the description for the given IDHMM
#    if idhmm in idhmm_to_desc and idhmm_to_desc[idhmm] != 'NA':
#        return f"{idhmm_to_desc[idhmm]} ({idhmm})"
#    else:
#        return 'unknown function ({idhmm})'

def add_hmmer_info(protsdict, information_proteins_hmmer, hmmer_databases):
    for contig in sorted(protsdict):
        for protein in sorted(protsdict[contig]):
            identifier = "NA"
            refinement = "NO_REFINEMENT"
            hmm_id     = "NA"
            hmm_pident = "NA"
            hmm_pcover = "NA"
            hmm_evalue = "NA"
            hmmer_info_dict = {}
            for database in information_proteins_hmmer:
                for identity in information_proteins_hmmer[database]:
                    if information_proteins_hmmer[database][identity][database]['name'] == protein:
                        hmmer_info_dict[database] = information_proteins_hmmer[database][identity][database]
            if len(hmmer_info_dict) > 0:
                identifier = "|".join([hmmer_info_dict[db]['name'] for db in hmmer_info_dict])
                refinement = "|".join([hmmer_info_dict[db]['db_name'] for db in hmmer_info_dict])
                hmm_id = "|".join("%s: %s" % (hmmer_info_dict[db]['db_name'], hmmer_info_dict[db]['idhmm']) for db in hmmer_info_dict)
                hmm_pident = "|".join([str(hmmer_info_dict[db]['pident']) for db in hmmer_info_dict])
                hmm_pcover = "|".join([str(hmmer_info_dict[db]['pcover']) for db in hmmer_info_dict])
                hmm_evalue = "|".join([str(hmmer_info_dict[db]['evalue']) for db in hmmer_info_dict])
            protsdict[contig][protein]['identifier'] = identifier
            protsdict[contig][protein]['refinement'] = refinement
            protsdict[contig][protein]['hmm_id'] = hmm_id
            protsdict[contig][protein]['hmm_pident'] = hmm_pident
            protsdict[contig][protein]['hmm_pcover'] = hmm_pcover
            protsdict[contig][protein]['hmm_evalue'] = hmm_evalue
    return protsdict

def create_protsdict(records_in_memory):
    protsdict = {}
    i = 0
    namelocuspat = re.compile(r'(\S+)\_\d+')
    for record in records_in_memory:
        dataprot = record.description.split(' # ')
        namelocus = re.match(namelocuspat, dataprot[0]).groups()[0]
        protsdict[namelocus] = {}
    return protsdict

def get_protein_info(record, information_proteins):
    dataprot = record.description.split(' # ')
    namelocuspat = re.compile(r'(\S+)\_\d+')
    namelocus = re.match(namelocuspat, dataprot[0]).groups()[0]
    modseq = str(record.seq).replace("X","")
    analysed_seq = ProteinAnalysis(modseq)
    protinfo = {}
    protinfo['length'] = len(record.seq)
    protinfo['isoelectricpoint'] = analysed_seq.isoelectric_point()
    protinfo['molweightkda'] = analysed_seq.molecular_weight()/1000.00
    protinfo['instability'] = analysed_seq.instability_index()
    protinfo['protein_id'] = dataprot[0]
    protinfo['strand'] = int(dataprot[3])
    protinfo['begin'] = int(dataprot[1])-1
    protinfo['end'] = int(dataprot[2])
    protinfo['translation'] = record.seq
    hypotheticalpat = re.compile(r'(?i)(((hypothetical|uncharacteri(z|s)ed|predicted))( phage)?( membrane)? protein|(ORF|(unnamed protein product|gp\d+|protein of unknown function|phage protein)))')
    try:
        if information_proteins[dataprot[0]]['descr'] == None:
            protinfo['descr'] = 'Hypothetical protein'
            protinfo['source'] = "NO_HIT"
            protinfo['pident'] = "NA"
            protinfo['pcover'] = "NA"            
            protinfo['evalue'] = "NA"            
        elif re.search(hypotheticalpat, information_proteins[dataprot[0]]['descr']):
            protinfo['descr'] = 'Conserved hypothetical protein'
            protinfo['source'] = "Homology"
            protinfo['pident'] = information_proteins[dataprot[0]]['pident']
            protinfo['pcover'] = information_proteins[dataprot[0]]['pcover']
            protinfo['evalue'] = information_proteins[dataprot[0]]['evalue']
        else:
            listdescr = information_proteins[dataprot[0]]['descr'].split(" ")
            del listdescr[0]
            protinfo['descr'] = " ".join(listdescr)
            protinfo['source'] = "Homology"
            protinfo['pident'] = information_proteins[dataprot[0]]['pident']
            protinfo['pcover'] = information_proteins[dataprot[0]]['pcover']
            protinfo['evalue'] = information_proteins[dataprot[0]]['evalue']
    except KeyError:
        protinfo['descr'] = 'Hypothetical protein'
        protinfo['source'] = "NO_HIT"
        protinfo['pident'] = "NA"
        protinfo['pcover'] = "NA"            
        protinfo['evalue'] = "NA"
    return namelocus, dataprot[0], protinfo

def update_protsdict(protsdict, records_in_memory, information_proteins):
    for record in records_in_memory:
        namelocus, dataprot, protinfo = get_protein_info(record, information_proteins)
        protsdict[namelocus][dataprot] = protinfo
    return protsdict

def write_protein_properties_table(protsdict, root_output):
    header = ["Contig", "Protein ID", "Start", "Stop", "Strand", "size_aa", "pI", "Mol_weight_kDa", "Instability_index", "Description", "Source", "Perc_ID", "Perc_Cov", "E-value", "HMMer", "Perc_ID", "Perc_Cov", "E-value"]
    with open(f"{root_output}.csv", "w") as tablefile:
        print("\t".join(header), file=tablefile)
        for contig, loci in sorted(protsdict.items()):
            for locus in sorted(loci, key=stringSplitByNumbers):
                refinement = protsdict[contig][locus].get('refinement', 'NO_REFINEMENT')
                data = [contig, locus, str(protsdict[contig][locus]['begin']), str(protsdict[contig][locus]['end']), str(protsdict[contig][locus]['strand']), str(protsdict[contig][locus]['length']), str(protsdict[contig][locus]['isoelectricpoint']), str(protsdict[contig][locus]['molweightkda']), str(protsdict[contig][locus]['instability']), protsdict[contig][locus]['descr'], protsdict[contig][locus]['source'], str(protsdict[contig][locus]['pident']), str(protsdict[contig][locus]['pcover']), str(protsdict[contig][locus]['evalue'])]
                if refinement != "NO_REFINEMENT":
                    data += [str(protsdict[contig][locus]['hmm_id']), str(protsdict[contig][locus]['hmm_pident']), str(protsdict[contig][locus]['hmm_pcover']), str(protsdict[contig][locus]['hmm_evalue'])]
                else:
                    data += ["NA", "NA", "NA", "NA"]
                print("\t".join(data), file=tablefile)

def add_CDS_features(whole_sequence, protsdict, record, args):
    for locus_protein in sorted(protsdict):
        if locus_protein == record.id:
            for protein in sorted(protsdict[locus_protein], key=stringSplitByNumbers):
                putative_start = int(protsdict[locus_protein][protein]['begin'])
                start_pos = ExactPosition(putative_start)
                end_pos = ExactPosition(protsdict[locus_protein][protein]['end'])
                feature_location = FeatureLocation(start_pos, end_pos, strand=protsdict[locus_protein][protein]['strand'])
                qualifiersgene = OrderedDict([('locus_tag', protsdict[locus_protein][protein]['protein_id'])])
                new_data_gene = SeqFeature(feature_location, type="gene", strand=protsdict[locus_protein][protein]['strand'], qualifiers=qualifiersgene)
                whole_sequence.features.append(new_data_gene)
                if not args.nohmmer:
                    if protsdict[locus_protein][protein]['identifier'] == "NA":
                        qualifiers = [('locus_tag', protsdict[locus_protein][protein]['descr']), ('product', protsdict[locus_protein][protein]['descr']), ('protein_id', protsdict[locus_protein][protein]['protein_id']), ('translation', protsdict[locus_protein][protein]['translation'])]
                    else:
                        qualifiers = [('locus_tag', protsdict[locus_protein][protein]['descr']), ('product', protsdict[locus_protein][protein]['descr']), ('protein_id', protsdict[locus_protein][protein]['protein_id']), ('note', '%s' % protsdict[locus_protein][protein]['hmm_id']), ('translation', protsdict[locus_protein][protein]['translation'])]
                else:
                    qualifiers = [('locus_tag', protsdict[locus_protein][protein]['descr']), ('product', protsdict[locus_protein][protein]['descr']), ('protein_id', protsdict[locus_protein][protein]['protein_id']), ('translation', protsdict[locus_protein][protein]['translation'])]
                feature_qualifiers = OrderedDict(qualifiers)
                new_data_cds = SeqFeature(feature_location, type="CDS", strand=protsdict[locus_protein][protein]['strand'], qualifiers=feature_qualifiers)
                whole_sequence.features.append(new_data_cds)

def add_tRNA_features(whole_sequence, tRNAdict, record):
    for locus_tRNA in sorted(tRNAdict):
        if locus_tRNA == record.id:
            for tRNA in sorted(tRNAdict[locus_tRNA], key=stringSplitByNumbers):
                start_pos = ExactPosition(tRNAdict[locus_tRNA][tRNA]['begin'])
                end_pos = ExactPosition(tRNAdict[locus_tRNA][tRNA]['end'])
                feature_location = FeatureLocation(start_pos, end_pos, strand=tRNAdict[locus_tRNA][tRNA]['strand'])
                new_data_gene = SeqFeature(feature_location, type="gene", strand=tRNAdict[locus_tRNA][tRNA]['strand'])
                whole_sequence.features.append(new_data_gene)
                qualifiers = [('product', tRNAdict[locus_tRNA][tRNA]['product'])]
                feature_qualifiers = OrderedDict(qualifiers)
                new_data_tRNA = SeqFeature(feature_location, type="tRNA", strand=tRNAdict[locus_tRNA][tRNA]['strand'], qualifiers=feature_qualifiers)
                whole_sequence.features.append(new_data_tRNA)

def add_tmRNA_features(whole_sequence, tmRNAdict, record):
    for locus_tmRNA in sorted(tmRNAdict):
        if locus_tmRNA == record.id:
            for tmRNA in sorted(tmRNAdict[locus_tmRNA], key=stringSplitByNumbers):
                start_pos = ExactPosition(tmRNAdict[locus_tmRNA][tmRNA]['begin'])
                end_pos = ExactPosition(tmRNAdict[locus_tmRNA][tmRNA]['end'])
                feature_location = FeatureLocation(start_pos, end_pos, strand=tmRNAdict[locus_tmRNA][tmRNA]['strand'])
                new_data_gene = SeqFeature(feature_location, type="gene", strand=tmRNAdict[locus_tmRNA][tmRNA]['strand'])
                whole_sequence.features.append(new_data_gene)
                qualifiers = [('product', tmRNAdict[locus_tmRNA][tmRNA]['product'])]
                feature_qualifiers = OrderedDict(qualifiers)
                new_data_tmRNA = SeqFeature(feature_location, type="tmRNA", strand=tmRNAdict[locus_tmRNA][tmRNA]['strand'], qualifiers=feature_qualifiers)
                whole_sequence.features.append(new_data_tmRNA)

def add_ncRNA_features(whole_sequence, elementsncRNA, record):
    for locus_ncRNA in sorted(elementsncRNA):
        if locus_ncRNA == record.id:
            for count, data in elementsncRNA[locus_ncRNA]['Other'].items():
                putative_start = data.get('begin')
                if putative_start is not None:
                    start_pos = ExactPosition(int(putative_start))
                    end_pos = ExactPosition(data['end'])
                    feature_location = FeatureLocation(start_pos, end_pos, strand=data['strand'])
                    new_data_gene = SeqFeature(feature_location, type="gene", strand=data['strand'])
                    whole_sequence.features.append(new_data_gene)
                    qualifiers = [('product', data['product']), ('note', 'RFAM: %s' % data['rfamcode'])]
                    feature_qualifiers = OrderedDict(qualifiers)
                    new_data_ncRNA = SeqFeature(feature_location, type="ncRNA", strand=data['strand'], qualifiers=feature_qualifiers)
                    whole_sequence.features.append(new_data_ncRNA)

def add_CRISPR_features(whole_sequence, information_CRISPR, record):
    for locus_CRISPR in sorted(information_CRISPR, key=stringSplitByNumbers):
        if locus_CRISPR == record.id:
            for CRISPR in sorted(information_CRISPR[locus_CRISPR], key=stringSplitByNumbers):
                putative_start = int(information_CRISPR[locus_CRISPR][CRISPR]['start'])
                start_pos = ExactPosition(putative_start)
                end_pos = ExactPosition(int(information_CRISPR[locus_CRISPR][CRISPR]['end']))
                feature_location = FeatureLocation(start_pos, end_pos)
                qualifiers = [('rpt_family', 'CRISPR'), ('rpt_type', 'direct'), ('rpt_unit_range', '%i..%i' % (putative_start, int(information_CRISPR[locus_CRISPR][CRISPR]['repeatend']))), ('rpt_unit_seq', information_CRISPR[locus_CRISPR][CRISPR]['repeatseq'])]
                feature_qualifiers = OrderedDict(qualifiers)
                new_data_CRISPRrepeat = SeqFeature(feature_location, type='repeat_region', qualifiers=feature_qualifiers)
                whole_sequence.features.append(new_data_CRISPRrepeat)

def write_genbank_file(filename, protsdict, tRNAdict, tmRNAdict, elementsncRNA, information_CRISPR, args, IUPAC=None):
    newgbk = "%s.gbk" % filename
    with open(filename, "r") as basefile, open(newgbk, "w") as outfile:
        if IUPAC:
            fileprocessed = SeqIO.parse(basefile, "fasta", IUPAC.ambiguous_dna)
        else:
            fileprocessed = SeqIO.parse(basefile, "fasta")
        for record in fileprocessed:
            whole_sequence = SeqRecord(record.seq)
            whole_sequence.id = str(record.id)
            whole_sequence.annotations['molecule_topology'] = genomeshape[record.id]['genomeshape']
            whole_sequence.annotations['molecule_type'] = "DNA"
            whole_sequence.annotations['data_file_division'] = args.typedata.upper()
            whole_sequence.annotations['date'] = strftime("%d-%b-%Y").upper()
            add_CDS_features(whole_sequence, protsdict, record, args)
            add_tRNA_features(whole_sequence, tRNAdict, record)
            add_tmRNA_features(whole_sequence, tmRNAdict, record)
            if args.norfam == False:
                add_ncRNA_features(whole_sequence, elementsncRNA, record)      
            add_CRISPR_features(whole_sequence, information_CRISPR, record)
            SeqIO.write(whole_sequence, outfile, "genbank")

def write_ptt_file(newfile):
    newgbk = "%s.gbk" % newfile
    for record in SeqIO.parse(newgbk, "gb"):
        record.features = [ptt for ptt in record.features if ptt.type == "CDS"]
        pttout = open("%s.ptt" % record.id, "w")
        pttout.write("{0} - 0..{1}\n".format(record.id, len(record)))
        pttout.write("{0} proteins\n".format(len(record.features)))
        pttout.write("Location\tStrand\tLength\tPID\tGene\tSynonym Code\tCOG\tProduct\n")
        strand = {1:'+', -1:'-'}
        for ptt in record.features:
            pttout.write("{0}\n".format("\t".join([str(ptt.location.start)+".."+str(ptt.location.end),strand[ptt.location.strand],str(abs(ptt.location.start-ptt.location.end)),'-',ptt.qualifiers["locus_tag"][0],ptt.qualifiers["locus_tag"][0],"-",ptt.qualifiers["product"][0]])))
        pttout.close()

## Preparing sequences for GenBank submission (Original code from Wan Yu's gbk2tbl.py script [https://github.com/wanyuac/BINF_toolkit/blob/master/gbk2tbl.py])
def create_genbank_submission_files(args, root_output):
    allowed_qualifiers = ['locus_tag', 'gene', 'product', 'pseudo', 'protein_id', 'gene_desc', 'old_locus_tag', 'note', 'inference', 'organism', 'mol_type', 'strain', 'sub_species', 'isolation-source', 'country']
    with open(args.modifiers, "r") as modifiers, open(f"{root_output}.gbk", "r") as genbank_fh, open(f"{root_output}.fasta", "w") as fasta_fh, open(f"{root_output}.tbl", "w") as feature_fh:
        info = modifiers.readline()
        records = list(SeqIO.parse(genbank_fh, 'genbank'))
        for record in records:
            if len(record) <= args.mincontigsize:
                print(f"WARNING: Skipping small contig {record.id}")
                continue
            record.description = f"{record.id} {info}"
            SeqIO.write([record], fasta_fh, 'fasta')
            print(f'>Feature {record.name}', file=feature_fh)
            for line in record.features:
                start, end = line.location.nofuzzy_start + 1, line.location.nofuzzy_end
                if line.strand == -1:
                    start, end = end, start + 1
                print(f'{start}\t{end}\t{line.type}', file=feature_fh)
                for key, values in line.qualifiers.items():
                    if key not in allowed_qualifiers:
                        continue
                    for v in values:
                        print(f'\t\t\t{key}\t{v}', file=feature_fh)

def remove_files():
    files_to_remove = ["temp.fasta", "CONTIGS_ALL.fasta", "temporal_circular.fasta", "crisprfile.txt", "pretemp.faa", "PROTS_FIRST_ROUND.faa", "PROTS_FIRST_ROUND.faa.csv"]
    if not args.norfam:
        files_to_remove.append("ncrnafile.csv")
    if not args.nohmmer:
        for tblfiles in sorted(glob.glob('PROTS_*.tbl')):
            files_to_remove.append(tblfiles)

    for splittedtrnafiles in sorted(glob.glob('trnafile_*.fasta')):
        files_to_remove.append(splittedtrnafiles)

    for splittedorffiles in sorted(glob.glob('orffile_*.faa')):
        files_to_remove.append(splittedorffiles)

    for splittedorf2files in sorted(glob.glob('orffile_*.fna')):
        files_to_remove.append(splittedorf2files)

    for tempfna in glob.glob("LOC_*.fna"):
        files_to_remove.append(tempfna)

    for tempgbk in glob.glob("LOC_*.gbk"):
        files_to_remove.append(tempgbk)

    for tempptt in glob.glob("LOC_*.ptt"):
        files_to_remove.append(tempptt)

    for file in files_to_remove:
        if os.path.isfile(file):
            os.remove(file)

if __name__ == "__main__":
    # Print the program header
    eprint("This is VIGA %s" % str(version))
    eprint("Written by Enrique Gonzalez Tortuero & Vimalkumar Velayudhan")
    eprint("Homepage is https://github.com/EGTortuero/viga")
    eprint("Local time: ", strftime("%a, %d %b %Y %H:%M:%S"))    

    # Interpret the main parameters
    parser = create_parser()
    add_basic_group(parser)
    add_advanced_general_group(parser)
    add_advanced_circularity_group(parser)
    add_advanced_repeats_group(parser)
    add_gcode_group(parser)
    add_type_group(parser)
    add_rfam_group(parser)
    add_diamond_group(parser)
    add_advanced_blast_group(parser)
    add_advanced_hmm_group(parser)
    args = parser.parse_args()
    
    # Load the program and evaluate the inputfiles
    root_output = check_and_set_default_databases(args)
    check_required_cmds(args)
    eprint("Data type is {0} and GenBank translation table no is {1}\n".format(args.typedata, args.gcode))
    record_iter = SeqIO.parse(open(args.inputfile, "r"),"fasta")
    rename_sequences(args, record_iter)

    # First step: predicting the contig shape
    starttime1 = time.time()
    eprint("\nPredicting the shape for all contigs using LASTZ")
    genomeshape = {}  # Initialize genomeshape for each new file
    with open("logfile.txt", "a") as logfile:
        logfile.write("\n#Contig\tShape\n")
        for newfile in sorted(glob.glob("LOC_*.fna")):
            Sequence, genomeshape = predict_shape(newfile, args, IUPAC, logfile)
            for record in Sequence:
                if genomeshape[record.id]['genomeshape'] == "circular":
                    circularize_sequence(record, newfile, args, IUPAC=IUPAC)
        cat_all(sorted(glob.glob('LOC_*.fna')), 'CONTIGS_ALL.fasta')
    endtime1 = time.time()
    durationtime1 = endtime1 - starttime1
    eprint("Done: shape prediction took %s seconds" % str(durationtime1))

    # Second step: identification of tRNAs
    starttime2 = time.time()
    tRNAdict = dict()
    tmRNAdict = dict()
    genetictable = "-gc%s" % str(args.gcode)
    eprint("\nRunning ARAGORN to predict tRNA-like sequences for all contigs")
    total_tRNA = 0
    for contigfile in sorted(glob.glob("LOC_*.fna")):
        num_tRNA = process_contig_for_tRNAs(contigfile, tRNAdict, tmRNAdict, genetictable)
        total_tRNA += num_tRNA
        eprint("Detected %i tRNAs in %s" % (num_tRNA, contigfile))
    endtime2 = time.time()
    durationtime2 = endtime2 - starttime2
    eprint("Done: tRNA and tmRNA detection took %s seconds" % str(durationtime2))

    # Third step: Predicting all ncRNA sequences (except rRNAs and tRNAs) (OPTIONAL: very SLOW step for genomic data) - Needed some optimisation?
    if args.norfam == False:
        starttime3 = time.time()
        eprint("\nIdentifying all other ncRNA (except rRNAs and tRNAs) for all contigs")
        output_file = run_cmscan(args.rfamdatabase, args.ncpus, "CONTIGS_ALL.fasta", "ncrnafile.csv")
        elementsncRNA = parse_and_read_ncrnafile("ncrnafile.csv")
        endtime3 = time.time()
        durationtime3 = endtime3 - starttime3
        eprint("Done: ncRNA detection took %s seconds" % str(durationtime3))

    # Fourth step: Predicting CRISPR repeats
    starttime4 = time.time()
    eprint("\nRunning PILER-CR to predict CRISPR repeats for all contigs")
    with open("/dev/null", "w") as apocalypse:
        subprocess.call(["pilercr", "-in", "CONTIGS_ALL.fasta", "-out", "crisprfile.txt", "-noinfo", "-minrepeat", str(args.minrepeat), "-maxrepeat", str(args.maxrepeat), "-minspacer", str(args.minspacer), "-maxspacer", str(args.maxspacer)], stderr=apocalypse)
    information_CRISPR = extract_crispr_data("crisprfile.txt")
    for seqid in SeqIO.parse("CONTIGS_ALL.fasta", "fasta"):
        num_CRISPRs = len(information_CRISPR.get(seqid.id, {}))
        eprint("Detected %i CRISPR repeats in %s" % (num_CRISPRs, seqid.id))
    endtime4 = time.time()
    durationtime4 = endtime4 - starttime4
    eprint("Done: CRISPR repeats detection took %s seconds" % str(durationtime4))

    # Fifth step: Predicting genes using PRODIGAL (Future step: adding Pyrodigal as an alternative)
    starttime5 = time.time()
    eprint("\nRunning Prodigal to predict the ORFs in all contigs")
    for contigfile in sorted(glob.glob("LOC_*.fna")):
        for record in SeqIO.parse(contigfile, "fasta"):
            orffile = "orffile_%s.faa" % record.id
            orffile2 = "orffile_%s.fna" % record.id
            length_contig = len(record.seq)
            predict_genes(contigfile, record, orffile, orffile2, length_contig, genomeshape, args)
    cat_all(sorted(glob.glob('orffile*.faa')), 'PROTS_FIRST_ROUND.faa')
    cat_all(sorted(glob.glob('orffile*.faa')), "%s.genes.faa" % root_output)
    cat_all(sorted(glob.glob('orffile*.fna')), "%s.genes.fna" % root_output)
    endtime5 = time.time()
    durationtime5 = endtime5 - starttime5
    eprint("Done: protein prediction took %s seconds" % str(durationtime5))

    # Sixth step: Predicting protein function based on homology using DIAMOND or BLAST
    starttime6 = time.time()
    if args.blastswitch == True:
        run_blast("PROTS_FIRST_ROUND.faa", args.blastdatabase, args.blastevalue, args.ncpus, args.blastexh)
        information_proteins = parse_blast_diamond_results("PROTS_FIRST_ROUND.faa.csv", args.blastwidthreshold, args.blastcovthreshold, args.blastevalue, "blast")
    else:
        eprint("\nRunning DIAMOND to predict the protein function according to homology inference using default parameters")
        run_diamond_blastp("PROTS_FIRST_ROUND.faa", "PROTS_FIRST_ROUND.faa.csv", args.diamonddatabase, args.diamondevalue, args.ncpus)
        information_proteins = parse_blast_diamond_results("PROTS_FIRST_ROUND.faa.csv", args.diamondwidthreshold, args.diamondcovthreshold, args.diamondevalue, "diamond")
    endtime6 = time.time()
    durationtime6 = endtime6 - starttime6
    eprint("Done: function prediction based on homology took %s seconds" % str(durationtime6))    

    # Sixth-extra step: Predicting the function of the proteins based on HMM predictions using HMMer 3.0 (Optional: SLOW step)
    if not args.nohmmer:
        hmmer_dbs = []
        if not args.norvdb:
            hmmer_dbs.append(('RVDB', args.rvdbdb, 'PROTS_RVDB.tbl'))
        if not args.novogs:
            hmmer_dbs.append(('VOGS', args.vogsdb, 'PROTS_VOGS.tbl'))
        if not args.nophrogs:
            hmmer_dbs.append(('PHROGS', args.phrogsdb, 'PROTS_PHROGS.tbl'))

        # Run HMMER and parse output files
        information_proteins_hmmer = {}
        for db_name, db_path, tbl_file in hmmer_dbs:
            starttime6a = time.time()
            eprint(f"\nRunning HMMER to enrich the annotations for all viral proteins using {db_name}")
            subprocess.call(['hmmsearch', '--cpu', str(args.ncpus), '--domtblout', tbl_file, '-E', str(args.hmmerevalue), '-o', '/dev/null', db_path, 'PROTS_FIRST_ROUND.faa'])
            information_proteins_hmmer[db_name] = parse_hmmer_file(tbl_file, args.hmmercovthreshold, args.hmmerevalue, args.hmmeridthreshold, db_name)
            endtime6a = time.time()
            durationtime6a = endtime6a - starttime6a
            eprint(f"Done: {db_name} prediction based on HMMER took {durationtime6a:.2f} seconds")

    # Seven step: Retrieving all the information for every protein
    records_in_memory = list(SeqIO.parse(open("PROTS_FIRST_ROUND.faa", "r"), "fasta"))
    protsdict = create_protsdict(records_in_memory)
    protsdict = update_protsdict(protsdict, records_in_memory, information_proteins)
    if args.nohmmer == False:
        hmmer_databases = {}
        for db_name, _, tbl_file in hmmer_dbs:
            hmmer_databases[db_name] = tbl_file
            protsdict = add_hmmer_info(protsdict, information_proteins_hmmer, hmmer_databases)

    # Eight step: Creating output files
    eprint("\nCreating all output files")
    write_protein_properties_table(protsdict, root_output) # Creating the CSV table with all protein statistics
    for newfile in sorted(glob.glob("LOC_*.fna")):
        if args.norfam == True:
            elementsncRNA = ""
        write_genbank_file(newfile, protsdict, tRNAdict, tmRNAdict, elementsncRNA, information_CRISPR, args, IUPAC)
        write_ptt_file(newfile)
    cat_all(sorted(glob.glob('LOC_*.fna.gbk')), "%s.gbk" % root_output) # Preparing the Genbank files
    cat_all(sorted(glob.glob('LOC_*ptt')), "%s.ptt" % root_output)      # Preparing the Protein Table (ppt table)
    with open(("%s.gff" % root_output), "w") as outgff, open("%s.gbk" % root_output, "r") as ingbk:
        GFF.write(SeqIO.parse(ingbk, "genbank"), outgff)                # Printing the GFF output
    create_genbank_submission_files(args, root_output)
        
    # Cleaning the intermediate files
    eprint("Cleaning all intermediate files")
    remove_files()
    eprint("Done.")
    
    # Final statement
    eprint("\nGenome annotation done!")
    eprint("The GenBank file is %s.gbk" % root_output)
    eprint("The GFF3 file is %s.gff" % root_output)
    eprint("The table file for GenBank submission is %s.tbl" % root_output)
    eprint("The FASTA file for GenBank submission is %s.fasta" % root_output)
    eprint("The table file with all protein information is %s.csv" % root_output)
