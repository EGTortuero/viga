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

# Creating databases folder
mkdir databases
cd databases

# Downloading RFAM db and formatting
echo "Downloading RFAM and formatting for its use in INFERNAL"
mkdir rfam
cd rfam
curl -O ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm &> /dev/null
cd ..
echo "Done"
echo ""

# Downloading all RefSeq_Viral_Proteins db
echo "Downloading RefSeq Viral Proteins"
curl -O https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
curl -O https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz
curl -O https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.protein.faa.gz
curl -O https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.protein.faa.gz
gunzip viral.1.protein.faa.gz
gunzip viral.2.protein.faa.gz
gunzip viral.3.protein.faa.gz
gunzip viral.4.protein.faa.gz
cat viral.1.protein.faa viral.2.protein.faa viral.3.protein.faa viral.4.protein.faa > refseq_viral_proteins.faa
rm viral.1.protein.faa viral.2.protein.faa viral.3.protein.faa viral.4.protein.faa
echo "Done"
echo ""

# Formatting RefSeq_Viral_Proteins for its use in DIAMOND
echo "Formatting the database to be used in DIAMOND"
mkdir RefSeq_Viral_DIAMOND
cd RefSeq_Viral_DIAMOND
cp ../refseq_viral_proteins.faa .
diamond makedb --in refseq_viral_proteins.faa -d refseq_viral_proteins &> /dev/null
rm refseq_viral_proteins.faa
cd ..
echo "Done"
echo ""

# Formatting RefSeq_Viral_Proteins for its use in BLAST
echo "Formatting the database to be used in BLAST"
mkdir RefSeq_Viral_BLAST
cd RefSeq_Viral_BLAST
cp ../refseq_viral_proteins.faa .
makeblastdb -in refseq_viral_proteins.faa -dbtype prot -out refseq_viral_proteins &> /dev/null
rm refseq_viral_proteins.faa
cd ..
rm refseq_viral_proteins.faa
echo "Done"
echo ""

# Downloading PVOGs
echo "Downloading PVOGs, VOGs and RVDB and formatting for its use in HMMer"
mkdir pvogs_rvdb
cd pvogs_rvdb
curl -O https://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz
#curl -O http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz
tar zxvf AllvogHMMprofiles.tar.gz &> /dev/null
{ echo AllvogHMMprofiles/*.hmm | xargs cat; } > pvogs_only.hmm
sed 's/NAME  VOG/NAME  PVOG/g' < pvogs_only.hmm > pvogs_only_mod.hmm
rm -rf AllvogHMMprofiles pvogs_only.hmm

# Downloading RVDB 
curl -O https://rvdb-prot.pasteur.fr/files/U-RVDBv23.0-prot.hmm.xz
unxz U-RVDBv23.0-prot.hmm.xz
mv U-RVDBv23.0-prot.hmm RVDB_23.0_only.hmm

# Downloading VOGs
curl -O http://fileshare.csb.univie.ac.at/vog/latest/vog.hmm.tar.gz
tar zxvf vog.hmm.tar.gz &> /dev/null
mkdir AllVOGHMMprofiles
mv VOG*hmm AllVOGHMMprofiles
{ echo AllVOGHMMprofiles/*.hmm | xargs cat; } > vog_only.hmm
rm -rf AllVOGHMMprofiles

## Downloading PHROGs
#curl -O https://phrogs.lmge.uca.fr/downloads_from_website/HMM_phrog.tar.gz
#tar zxvf HMM_phrog.tar.gz &> /dev/null
#{ echo HMM_phrog/*.hmm | xargs cat; } > phrogs_only.hmm
#

# Formatting the database
cat pvogs_only_mod.hmm RVDB_23.0_only.hmm vog_only.hmm > pvogs_vogs_RVDB.hmm
hmmpress -f pvogs_vogs_RVDB.hmm &> /dev/null
cd ../..
echo "Done"
