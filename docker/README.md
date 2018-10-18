A Docker image for VIGA
=======================
This is a [Docker](https://www.docker.com) image for
running [VIGA](https://github.com/EGTortuero/viga) in a
container. If you have installed Docker, you can pull this
image using:

	docker pull vimalkvn/viga

In addition to the **VIGA** script from the Github repository,
this image has all dependencies pre-installed — Aragorn, HMMER,
Prodigal, BioPython, LASTZ, Tandem and Inverted repeats finder

## Note regarding databases

If you have all the necessary databases downloaded and formatted
correctly, you can start using this container immediately.
The [create_dbs.sh](https://github.com/EGTortuero/viga) script
included in this repository automates this process of
downloading and formatting databases.

The methods below assume these databases (`--blastdb`,
`--rfamdb` and `--hmmerdb`) are located in
`/data/databases`. If this is not the
case, please modifythe wrapper script or the command
accordingly. A simple approach would be to save all
databases in one location and export the top level directory
using the `-v` option of the `docker` command.

## Using this Docker image
Once databases are setup, select from one of the following
methods to run a container — using the wrapper
script or running a container directly.

In both cases, input files are assumed to be present in
the directory where the command is run.

Other VIGA program options can be specified and are
sent as is.

### Option 1. Use the wrapper script (recommended)

**Note:** This script requires **sudo** access as it calls the
`docker` command. Administrators can restrict sudo access to
just this script in `/etc/sudoers` if required.

Download the [run-viga](https://raw.githubusercontent.com/vimalkvn/viga/master/docker/run-viga)
wrapper script available from this repository, save it in a
location accessible in the **PATH**, for example, `/usr/local/bin`,
make it executable and then run it like this:

	sudo run-viga --input rubella.fasta \
	--diamonddb /data/databases/refseq_viral_proteins \
	--blastdb /data/databases/refseq_viral_proteins \
	--hmmerdb /data/databases/pvogs.hmm \
	--rfamdb /data/databases/Rfam.cm \
	--modifiers modifiers.txt

### Option 2. Run a container directly using docker

**Note:** The user(s) running this command should be in
the **docker** group.

	docker run --rm \
	-e LOGNAME=$(logname) \
	-e USER=$(logname) \
	-u ${UID}:${UID} \
	-v /data/databases:/data/databases:ro \
	-v $(pwd):/wdir \
	-w /wdir \
	vimalkvn/viga \
	python /program/VIGA.py \
	--input rubella.fasta \
	--diamonddb /data/databases/refseq_viral_proteins \
	--blastdb /data/databases/refseq_viral_proteins \
	--hmmerdb /data/databases/pvogs.hmm \
	--rfamdb /data/databases/Rfam.cm \
	--modifiers modifiers.txt

The `-e` and `-u` options are used so that the container
has permissions to read input files and the generated output
files will be owned by the user running the container and not
root.
