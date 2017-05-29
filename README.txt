NAME: VirAnnot (virannot.py)

AUTHOR: Enrique Gonzalez-Tortuero <enrique.gonzaleztortuero@ucc.ie>

DESCRIPTION: VirAnnot is a script written in Python 2.7 that annotates viral 
genomes automatically (using a de novo algorithm) and predict the function of
their proteins using BLAST and HHSEARCH/HHPRED.

REQUIREMENTS:

Before using this script, the following Python modules and programs should be installed:

* Python modules:
	- BCBio (https://github.com/chapmanb/bcbio-nextgen)
	- Biopython (Bio module; Cock et al. 2009)

* Programs:
	- GNU Parallel (Tange 2011): it is used to parallelize BLAST and HHSUITE. The program is publicly available at https://www.gnu.org/software/parallel/ under the GPLv3 licence.
	- LASTZ (Harris 2007): it is used to predict the circularity of the contigs. The program is publicly available at https://github.com/lastz/lastz under the MIT licence.
	- Prodigal (Hyatt et al. 2010): it is used to predict the ORFs. When the contig is smaller than 20,000 bp, MetaProdigal (Hyatt et al. 2012) is automatically activated instead of normal Prodigal. This program is publicly available at https://github.com/hyattpd/prodigal/releases/ under the GPLv3 licence.
	- BLAST+ (Camacho et al. 2008): it is used to predict the function of the predicted proteins according to homology. This suite is publicly available at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ under the GPLv2 licence. Databases are available at ftp://ftp.ncbi.nlm.nih.gov/blast/db/
	- HHSUITE (Söding 2005): it is used to predict the function of the predicted proteins according to Hidden Markov Model-Hidden Markov Model (HMM-HMM) comparisons. First, the sequences are aligned against reference databases using HHblits (Remmert et al. 2011) and, then, the resulting multiple sequence alignment is converted into a HMM and compared with known HMM databases using HHsearch/HHpred (Hildebrand et al. 2009). This suite is publicly available at https://github.com/soedinglab/hh-suite under the GPLv3 licence. Databases are available at http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/
	- ARAGORN (Laslett and Canback 2004): it is used to predict tRNA sequences in the contig. This program is publicly available at http://mbio-serv2.mbioekol.lu.se/ARAGORN/ under the GPLv2 licence.
	- Tandem Repeats Finder (TRF; Benson 1999): it is used to predict the tandem repeats in your contig. This program is freely available at https://tandem.bu.edu/trf/trf.html under a custom licence.
	- Inverted Repeats Finder (IRF; Warburton et al. 2004): it is used to predict the inverted repeats in your contig. This program is freely available at https://tandem.bu.edu/irf/irf.download.html under a custom licence.

Although you can install the programs manually, we strongly recommend the use of the Docker image to create an environment for virannot. The link to the Docker image is:
	https://github.com/vimalkvn/sysadminbio/tree/master/docker-images/virannot

However, you will need to download the databases for BLAST and HHSUITE:

BLAST: https://ftp.ncbi.nlm.nih.gov/blast/db/
HHSUITE: http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/

Note that this bioinformatic pipeline only takes protein databases (i.e. "nr", "swissprot"...)!

PARAMETERS:

The program has the following two kind of arguments:

	Mandatory parameters
	--------------------

--input FASTAFILE					Input file as a nucleotidic FASTA file. It can contains multiple sequences (e.g. metagenomic contigs)
--blastdb BLASTDB					BLAST database that will be used for the protein function prediction. The database MUST be for amino acids.
--hhblitsdb HHBLITSDB					HHblits database that will be used for the first step of the protein function prediction. In this case, HHBLITSDB should be in the format "/full/path/to/db1/db1 (without the extension _a3m_db)"
--hhsearchdb HHSEARCHDB [HHSEARCHDB ...]		HHsearch/HHpred database/s that will be used for the second step of the protein function prediction. You can use more than a single database for the analysis. In that case, it will be in the format "/full/path/to/db1/db1_hhm_db /full/path/to/db2/db2_hhm_db"
--modifiers TEXTFILE					Input file as a plain text file with the modifiers per every FASTA header according to SeqIn (https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html). All modifiers must be written in a single line and are separated by a single space character. No space should be placed besides the = sign. For example: [organism=Serratia marcescens subsp. marcescens] [sub-species=marcescens] [strain=AH0650_Sm1] [moltype=DNA] [tech=wgs] [gcode=11] [country=Australia] [isolation-source=sputum]. This line will be copied and printed along with the record name as the definition line of every contig sequence.

	Advanced parameters
	-------------------

--readlength INT					Read length for the circularity prediction (default: 101 bp)
--out OUTPUTNAME					Name of the outputs files without extensions, as the program will add them automatically. By default, the program will use the input name as the output.
--locus STRING						Name of the contigs. If the input is a multiFASTA file, please put a general name as the program will add the number of the contig at the end of the name. By default, the name of the contigs will be "LOC".
--threads INT						Number of threads/CPUs. By default, the program will use 1 CPU.
--noparallel						Using multithreading BLAST and HHPRED instead of using parallel BLAST and HHPRED. Only recommendable when GNU Parallel is not installed. By default, this option is disabled.
--gff 							Printing the output as a General Feature Format (GFF) version 3. It is a flat table file with contains 9 columns of data (see http://www.ensembl.org/info/website/upload/gff3.html for more information). By default, the program will not print the GFF3 file (--gff False).
--blastevalue FLOAT					BLAST e-value threshold. By default, the threshold will be 1e-05.
--hhsuiteevalue FLOAT					HHSUITE e-value threshold. By default, the threshold is 1e-03.
--typedata CON|VRL|PHG					GenBank Division: One of the following codes:
								VRL - Eukaryotic/Archaea virus)
								PHG - Phages
								CON - Contig
							By default, the program will consider every sequence as a contig (CON)
--gcode NUMBER						Number of GenBank translation table. At this moment, the available options are:
								1  - Standard genetic code [Eukaryotic]
								11 - Bacteria/Archaea/Phages)
								4  - Mycoplasma/Spiroplasma
								25 - Gracilibacteria & Candidate division SR1
								6  - Protozoa [nuclear]
							By default, the program will use the translation table no. 11
--mincontigsize INT					Minimum contig length to be considered in the final files. By default, the program only consider from 200 bp.
--idthr FLOAT						Identity threshold to consider that a protein belong to a specific hit. By default, the threshold is 50.0 %
--coverthr FLOAT					Coverage threshold to consider that a protein belong to a specific hit. By default, the threshold is 50.0 %
--maxfilt INT (>100)					Max number of hits allowed to pass 2nd prefilter in HHblits. By default, the max number of hits is 20000
--neffmax FLOAT [1.00-20.00]				Skip further search iterations in HHblits when diversity of query MSA becomes larger than this value. By default, this value is 10
--diffid FLOAT (>0.01)					Max allowed difference between the ID percentages of BLAST and HHsearch/HHpred. By default, the allowed difference is 5.00 % and we do not recommended to change such value.
--blastexh						Use of exhaustive BLAST to predict the proteins by homology according to Fozo et al. (2010). In this case, the search will be done using a word size of 2, a gap open penalty of 8, a gap extension penalty of 2, the PAM70 matrix instead of the BLOSUM62 and no compositional based statistics. This method is more accurate to predict the functions of the proteins but it is slower than BLAST default parameters. By default, exhaustive BLAST is disabled.
--hhsearchexh						Use of exhaustive HHsearch/HHpred to tackle proteins of unknown function according to Fidler et al. (2016). In this case, the MSA is created with HHblits using 8 iterations and with a dynamic progressive e-value acceptance from 1e-03 to 1e-02. After that, HHsearch/HHpred use the Maximum ACcuracy (MAC) alignment algorithm to re-align the MSA and an e-value threshold of 1e-02. This method allows to predict a putative function when there are a lot of unknown function proteins but it is slower than HHsearch/HHpred default parameters. By default, exhaustive HHsearch/HHpred is disabled.

REFERENCES:

	- Benson G (2008) "Tandem repeats finder: a program to analyze DNA sequences." Nucleic Acids Research 27, 573–80.
	- Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10: 421.
	- Fidler DR, Murphy SE, Courtis K, Antonodiou P, El-Tohamy R, Ient J, Levine TP (2016) Using HHsearch to tackle proteins of unknown function: a pilot study with PH domains. Traffic 17: 1214–26.
	- Fozo EM, Makarova KS, Shabalina SA, Yutin N, Koonin EV, Storz G (2010) Abundance of type I toxin-antitoxin systems in bacteria: searches for new candidates and discovery of novel families. Nucleic Acids Research 38: 3743-59.
	- Harris RS (2007) Improved pairwise alignment of genomic DNA. Ph.D. Thesis, The Pennsylvania State University. 
	- Hildebrand A, Remmert A, Biegert A, Söding J (2009) Fast and accurate automatic structure prediction with HHpred. Proteins 77: 128-32.
	- Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11: 119.
	- Hyatt D, Locascio PF, Hauser LJ, Uberbacher EC (2012) Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics 28: 2223-30.
	- Laslett D, Canback B (2004) ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences. Nucleic Acids Research 32, 11–16.
	- Remmert M, Biegert A, Hauser A, Söding J (2011) HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment. Nature Methods 9: 173-5. 
	- Söding J (2005) Protein homology detection by HMM-HMM comparison. Bioinformatics 21: 951-60.
	- Tange O (2011) GNU Parallel - The Command-Line Power Tool. ;login: The USENIX Magazine 36:42-7.
	- Warburton PE, Giordano J, Cheung F, Gelfand Y, Benson G (2004) Inverted repeat structure of the human genome: The X-chromosome contains a preponderance of large, highly homologous inverted repeats that contain testes genes. Genome Research 14: 1861-9.

HISTORY: 

v 0.1.0 - Current version of the program.
