#! /bin/bash

## Installing databases

# Checking if all required programs are installed
if ! foobar_loc="$(type -p "cmpress")" || [[ -z $foobar_loc ]]; then
	echo "You need to run first the installer.sh script to download all databases. You need to install INFERNAL!"
	exit
fi
if ! foobar_loc="$(type -p "diamond")" || [[ -z $foobar_loc ]]; then
	echo "You need to run first the installer.sh script to download all databases. You need to install DIAMOND!"
	exit
fi
if ! foobar_loc="$(type -p "makeblastdb")" || [[ -z $foobar_loc ]]; then
	echo "You need to run first the installer.sh script to download all databases. You need to install BLAST+!"
	exit
fi

# Creating databases folder
mkdir databases
cd databases

# Downloading RFAM db and formatting
echo "Downloading RFAM and formatting for its use in INFERNAL"
mkdir rfam
cd rfam
curl -O ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz &> /dev/null
gunzip Rfam.cm.gz
cmpress Rfam.cm &> /dev/null
cd ..
echo "Done"
echo ""

# Downloading all RefSeq_Viral_Proteins db
echo "Downloading RefSeq Viral Proteins"
curl -O ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz &> /dev/null
curl -O ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz &> /dev/null
curl -O ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.protein.faa.gz &> /dev/null
gunzip viral.1.protein.faa.gz
gunzip viral.2.protein.faa.gz
gunzip viral.3.protein.faa.gz
cat viral.1.protein.faa viral.2.protein.faa viral.3.protein.faa > refseq_viral_proteins.faa
rm viral.1.protein.faa viral.2.protein.faa viral.3.protein.faa
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

# Downloading PVOGs and formatting 
echo "Downloading PVOGs and formatting for its use in HMMer"
mkdir pvogs
cd pvogs
curl -O ftp://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz &> /dev/null
tar zxvf AllvogHMMprofiles.tar.gz &> /dev/null
{ echo AllvogHMMprofiles/*.hmm | xargs cat; } > pvogs.hmm
rm -rf AllvogHMMprofiles
cd ../..
echo "Done"