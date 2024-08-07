#! /bin/bash

## Installing databases
echo "VIGA installation script"
echo "------------------------"
echo ""
echo "This script will install all required software to run VIGA before its execution."
echo "All programs will be installed in:"
echo $PWD/programs/bin
#echo "In case that you are running non-Debian-based UNIX distros, please replace all"
#echo "lines with 'apt' by 'dnf' in this script!"

# Creation of the folder "programs"
mkdir programs
cd programs/
mkdir bin

# Checking if all required programs are installed
echo ""
echo "Installing LASTZ..."
wget https://github.com/lastz/lastz/archive/refs/tags/1.04.22.tar.gz &>/dev/null
tar zxvf 1.04.22.tar.gz	&>/dev/null
cd lastz-1.04.22/src
sed 's/definedForAll = -Wall -Wextra -Werror -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE/definedForAll = -Wall -Wextra -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE/' Makefile > Makefile2 # It will remove all warnings as errors(!)
mv Makefile2 Makefile
cd ..
make &>/dev/null
make test &>/dev/null
mv src/lastz ../bin/
mv src/lastz_D ../bin/
echo "Done."
cd ..

echo ""
echo "Installing Aragorn..."
mkdir aragorn
cd aragorn
wget http://www.trna.se/ARAGORN/Downloads/aragorn1.2.41.c &>/dev/null
gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.41.c &>/dev/null
cp aragorn ../bin/
echo "Done."
cd ..

echo ""
echo "Installing Infernal..."
wget http://eddylab.org/infernal/infernal-1.1.5.tar.gz &>/dev/null
tar xvfz infernal-1.1.5.tar.gz &>/dev/null
cd infernal-1.1.5
./configure --prefix $PWD/../ &>/dev/null
make &>/dev/null
make check &>/dev/null
make install &>/dev/null
echo "Done."
cd ..

echo ""
echo "Installing PILER-CR..."
wget https://www.drive5.com/pilercr/pilercr1.06.tar.gz &>/dev/null
tar xvfz pilercr1.06.tar.gz &>/dev/null
cd pilercr1.06
rm -rf pilercr pilercr.exe
make &>/dev/null
cp pilercr ../bin
echo "Done."
cd ..

echo ""
echo "Installing Prodigal..."
git clone https://github.com/hyattpd/Prodigal.git &>/dev/null
cd Prodigal
make &>/dev/null
cp prodigal ../bin
echo "Done."
cd ..

echo ""
echo "Installing Prodigal-GV..."
git clone https://github.com/apcamargo/prodigal-gv.git &>/dev/null
cd prodigal-gv
make &>/dev/null
cp prodigal-gv ../bin
echo "Done."
cd ..

echo ""
echo "Installing Diamond..."
wget https://github.com/bbuchfink/diamond/archive/refs/tags/v2.1.8.tar.gz &>/dev/null
tar xzfv v2.1.8.tar.gz &>/dev/null
cd diamond-2.1.8
mkdir bin
cd bin
cmake -DCMAKE_INSTALL_PREFIX=${PWD}/../.. .. &>/dev/null
make -j4 &>/dev/null
make install &>/dev/null
echo "Done."
cd ../..

echo ""
echo "Installing BLAST (It will take time)..."
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz.md5 &>/dev/null
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz &>/dev/null
md5sum --quiet --check ncbi-blast-2.15.0+-x64-linux.tar.gz.md5
tar xvfz ncbi-blast-2.15.0+-x64-linux.tar.gz &>/dev/null
cp -r ncbi-blast-2.15.0+-x64-linux/bin/* bin/
echo "Done."

echo ""
echo "Installing HMMER 3..."
wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz &>/dev/null
tar zxf hmmer-3.4.tar.gz &>/dev/null
cd hmmer-3.4
./configure --prefix $PWD/../ &>/dev/null
make &>/dev/null
make check &>/dev/null
make install &>/dev/null
cd easel
make install &>/dev/null
echo "Done."
cd ..

#echo ""
#echo "Installing HHSuite..."
#git clone https://github.com/soedinglab/hh-suite.git
#mkdir -p hh-suite/build
#cd hh-suite/build
#cmake -DCMAKE_INSTALL_PREFIX=$PWD/../../ .. &>/dev/null
#make -j 4 &>/dev/null
#make install &>/dev/null
#echo "Done."
#cd ..

# Cleaning the folder to harbour only the binary files and needed libraries
echo ""
echo "Cleaning all temporary folders and downloaded tarballs"
rm *.tar.gz
rm ncbi-blast-2.14.1+-src.tar.gz.md5
rm -rf Prodigal/ aragorn/ diamond-2.1.8/ infernal-1.1.5/ lastz-1.04.22/ pilercr1.06/ 
rm -rf ncbi-blast-2.14.1+-src/ hmmer-3.4/ 
#rm hh-suite/

# Creating the databases
#echo ""
#echo "Creating the corresponding databases:"
#sh create_dbs.sh
