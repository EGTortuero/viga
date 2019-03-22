FROM ubuntu:trusty
LABEL maintainer="Vimalkumar Velayudhan <vimalkumarvelayudhan@gmail.com>"
ENV VIGA_VERSION="0.11.0"
ENV VIGA_COMMIT_ID="072DC21"
ENV LASTZ_VERSION="1.04.00"
ENV DIAMOND_VERSION="0.9.22"

# Install packages available in repositories
RUN sed -i 's/^# \(deb http.*trusty-backports main restricted\)/\1 universe/' \
/etc/apt/sources.list && \
apt-get update && \
apt-get install -y --no-install-recommends \
aragorn \
build-essential \
cmake \
hmmer \
infernal \
libc6-i386 \
libpthread-stubs0-dev \
ncbi-blast+ \
prodigal \
python-biopython \
python-pip \
python-scipy \
python-six \
wget \
zlib1g-dev && \
pip install bcbio-gff

# Download, compile and install lastz
RUN wget https://github.com/lastz/lastz/archive/${LASTZ_VERSION}.tar.gz && \
tar -xf ${LASTZ_VERSION}.tar.gz && \
cd lastz-${LASTZ_VERSION} && \
# remove -Werror from Makefile to fix compile errors
sed -i 's/-Werror //' src/Makefile && \
make && \
install -m 0755 src/lastz /usr/local/bin/ && \
install -m 0755 src/lastz_D /usr/local/bin/ && \
cd .. && rm -rf lastz-*

# Download, compile and install pilercr
RUN wget http://www.drive5.com/pilercr/pilercr1.06.tar.gz && \
tar -xf pilercr1.06.tar.gz && \
cd pilercr1.06 && \
make && \
install -m 0755 pilercr /usr/local/bin/ && \
cd .. && rm -rf pilercr1.06*

# Download Inverted Repeats Finder 3.07 and Tandem Repeats Finder 4.09
RUN wget -O irf https://tandem.bu.edu/irf/downloads/irf307.linux.exe && \
wget -O trf https://tandem.bu.edu/trf/downloads/trf409.linux64 && \
install -m 0755 irf /usr/local/bin/ && \
install -m 0755 trf /usr/local/bin/ && \
rm irf trf

# Download and compile diamond
RUN wget https://github.com/bbuchfink/diamond/archive/v${DIAMOND_VERSION}.tar.gz && \
tar -xf v${DIAMOND_VERSION}.tar.gz && \
cd diamond-${DIAMOND_VERSION} && \
mkdir bin && cd bin && \
cmake .. && \
make && make install

RUN apt-get -y purge --auto-remove python-pip \
	build-essential cmake libpthread-stubs0-dev \
	zlib1g-dev && rm -rf /var/lib/apt/lists/*

# Download viga script
RUN mkdir /program && \
	wget -P /program https://raw.githubusercontent.com/EGTortuero/viga/${VIGA_COMMIT_ID}/VIGA.py && \
	chmod +x /program/VIGA.py && apt-get -y purge wget

# create and set user for container
RUN useradd -r -s /sbin/nologin viga
USER viga

CMD printf "%s\n" \
	"This image runs VIGA version ${VIGA_VERSION} from Github commit ${VIGA_COMMIT_ID}" \
	" " \
	"Please run the container like this:" \
	"  docker run --rm vimalkvn/viga python /program/VIGA.py" \
	" " \
	"or use the run-viga wrapper script"
