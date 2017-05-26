# VirAnnot
De novo automatic Viral Annotation

VirAnnot is a script written in Python 2.7 that annotates viral genomes automatically (using a de novo algorithm) and predict the function of their proteins using BLAST and HHSEARCH/HHPRED.

REQUIREMENTS:

Before using this script, the following Python modules and programs should be installed:

* Python modules:
	- BCBio (https://github.com/chapmanb/bcbio-nextgen)
	- Biopython (Bio module; Cock et al. 2009)

* Programs:
	- LASTZ (Harris 2007): it is used to predict the circularity of the contigs. The program is publicly available at https://github.com/lastz/lastz under the MIT licence.
	- Prodigal (Hyatt et al. 2010): it is used to predict the ORFs. When the contig is smaller than 20,000 bp, MetaProdigal (Hyatt et al. 2012) is automatically activated instead of normal Prodigal. This program is publicly available at https://github.com/hyattpd/prodigal/releases/ under the GPLv3 licence.
	- BLAST+ (Camacho et al. 2008): it is used to predict the function of the predicted proteins according to homology. This suite is publicly available at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ under the GPLv2 licence. Databases are available at ftp://ftp.ncbi.nlm.nih.gov/blast/db/
	- HHSUITE (Söding 2005): it is used to predict the function of the predicted proteins according to Hidden Markov Model-Hidden Markov Model (HMM-HMM) comparisons. First, the sequences are aligned against reference databases using HHblits (Remmert et al. 2011) and, then, the resulting multiple sequence alignment is converted into a HMM and compared with known HMM databases using HHsearch/HHpred (Hildebrand et al. 2009). This suite is publicly available at https://github.com/soedinglab/hh-suite under the GPLv3 licence. Databases are available at http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/
	- ARAGORN (Laslett and Canback 2004): it is used to predict tRNA sequences in the contig. This program is publicly available at http://mbio-serv2.mbioekol.lu.se/ARAGORN/ under the GPLv2 licence.
	- Tandem Repeats Finder (TRF; Benson 1999): it is used to predict the tandem repeats in your contig. This program is freely available at https://tandem.bu.edu/trf/trf.html under a custom licence.
	- Inverted Repeats Finder (IRF; Warburton et al. 2004): it is used to predict the inverted repeats in your contig. This program is freely available at https://tandem.bu.edu/irf/irf.download.html under a custom licence.

Although you can install the programs manually, we strongly recommend the use of the Docker image to create an environment for virannot. The link to the Docker image is https://github.com/vimalkvn/sysadminbio/tree/master/docker-images/virannot

However, you will need to download the databases for BLAST and HHSUITE:
* BLAST: https://ftp.ncbi.nlm.nih.gov/blast/db/
* HHSUITE: http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/

Note that this bioinformatic pipeline only takes protein databases (i.e. "nr", "swissprot"...)!

When using this program, you must to cite their use:

HISTORY OF THE SOURCE CODE:

v 0.1.0 - Current version of the program.
