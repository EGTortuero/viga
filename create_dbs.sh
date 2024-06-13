#! /bin/bash

## Installing databases
echo "VIGA Databases installation and update script"
echo "------------------------"
echo ""
echo "This script will install and/or update all databases needed to run VIGA before its execution."
echo "All programs will be installed in:"
echo $PWD/databases
#echo "In case that you are running non-Debian-based UNIX distros, please replace all"
#echo "lines with 'apt' by 'dnf' in this script!"

# Checking if all required programs are installed
if ! foobar_loc="$(type -p "esl-reformat")" || [[ -z $foobar_loc ]]; then
	echo "You need to run first the install.sh script to download all databases. You need to install EASEL (HMMER 3)!"
	exit
fi
if ! foobar_loc="$(type -p "cmpress")" || [[ -z $foobar_loc ]]; then
	echo "You need to run first the install.sh script to download all databases. You need to install INFERNAL!"
	exit
fi
if ! foobar_loc="$(type -p "diamond")" || [[ -z $foobar_loc ]]; then
	echo "You need to run first the install.sh script to download all databases. You need to install DIAMOND!"
	exit
fi
if ! foobar_loc="$(type -p "makeblastdb")" || [[ -z $foobar_loc ]]; then
	echo "You need to run first the install.sh script to download all databases. You need to install BLAST+!"
	exit
fi
if ! foobar_loc="$(type -p "hmmpress")" || [[ -z $foobar_loc ]]; then
	echo "You need to run first the install.sh script to download all databases. You need to install HMMER 3!"
	exit
fi

# Creating databases folder
mkdir databases
cd databases

# Downloading RFAM db and formatting
echo "Downloading RFAM and formatting for its use in INFERNAL"
mkdir rfam
cd rfam
curl -O ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm
cd ..
echo "Done"
echo ""

# Downloading all RefSeq_Viral_Proteins db
echo "Downloading RefSeq Viral Proteins"
curl -O https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
gunzip viral.1.protein.faa.gz
mv viral.1.protein.faa refseq_viral_proteins.faa
echo "Done"
echo ""

# Formatting RefSeq_Viral_Proteins for its use in DIAMOND
echo "Formatting the database to be used in DIAMOND"
mkdir RefSeq_Viral_DIAMOND
cd RefSeq_Viral_DIAMOND
cp ../refseq_viral_proteins.faa .
diamond makedb --in refseq_viral_proteins.faa -d refseq_viral_proteins
rm refseq_viral_proteins.faa
cd ..
echo "Done"
echo ""

# Formatting RefSeq_Viral_Proteins for its use in BLAST
echo "Formatting the database to be used in BLAST"
mkdir RefSeq_Viral_BLAST
cd RefSeq_Viral_BLAST
cp ../refseq_viral_proteins.faa .
makeblastdb -in refseq_viral_proteins.faa -dbtype prot -out refseq_viral_proteins
rm refseq_viral_proteins.faa
cd ..
rm refseq_viral_proteins.faa
echo "Done"
echo ""

# Downloading RVDB
echo "Downloading RVDB and formatting for its use in HMMer"
mkdir pvogs_rvdb
cd pvogs_rvdb
curl -O https://rvdb-prot.pasteur.fr/files/U-RVDBv28.0-prot.hmm.xz
curl -O https://rvdb-prot.pasteur.fr/files/U-RVDBv28.0-prot-hmm-txt.tar.xz
unxz U-RVDBv28.0-prot.hmm.xz
unxz U-RVDBv28.0-prot-hmm-txt.tar.xz
tar xvf U-RVDBv28.0-prot-hmm-txt.tar
#{ echo annot/*.txt | xargs cat; } > U-RVDBv26.0.txt # Pending an script to transform these files into a table with annotation ID and 
mv U-RVDBv28.0-prot.hmm RVDB_28.0_only.hmm

# Downloading VOGs and VFAM
echo "Downloading VOGs and formatting for its use in HMMer"
curl -O https://fileshare.lisc.univie.ac.at/vog/vog224/vog.hmm.tar.gz
curl -O https://fileshare.lisc.univie.ac.at/vog/vog224/vog.annotations.tsv.gz
tar zxvf vog.hmm.tar.gz &> /dev/null
gunzip vog.annotations.tsv.gz
{ echo hmm/*.hmm | xargs cat; } > vog_only.hmm
rm -rf hmm/

echo "Downloading VFAM and formatting for its use in HMMer"
curl -O https://fileshare.lisc.univie.ac.at/vog/vog224/vfam.hmm.tar.gz
curl -O https://fileshare.lisc.univie.ac.at/vog/vog224/vfam.annotations.tsv.gz
tar zxvf vfam.hmm.tar.gz &> /dev/null
gunzip vfam.annotations.tsv.gz
{ echo hmm/*.hmm | xargs cat; } > vfam_only.hmm
rm -rf hmm/

# Downloading PHROGs
echo "Downloading PHROGs and formatting for its use in HMMer"
curl -O https://phrogs.lmge.uca.fr/downloads_from_website/MSA_phrogs.tar.gz
curl -O https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv
tar zxvf MSA_phrogs.tar.gz &> /dev/null
cd MSA_Phrogs_M50_FASTA/
for X in $(ls *.fma | sed "s/\.fma//"); 
	do 
	esl-reformat -o $X.sto stockholm $X.fma;
	hmmbuild --cpu 8 $X.hmm $X.sto;
	done
cd ..
{ echo MSA_Phrogs_M50_FASTA/*.hmm | xargs cat; } > phrogs_only.hmm
rm -rf MSA_Phrogs_M50_FASTA

# Formatting the database
#cat pvogs_only_mod.hmm RVDB_24.1_only.hmm vog_only.hmm > pvogs_vogs_RVDB.hmm
#cat vog_only.hmm RVDB_26.0_only.hmm phrogs_only.hmm > vogs_RVDB_phrogs.hmm
#hmmpress -f vogs_RVDB_phrogs.hmm
#cd ../..
echo "Done"
