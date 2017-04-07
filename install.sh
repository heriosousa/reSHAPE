#!/bin/bash

BASEDIR=$(pwd)

#Installing a5-miseq
echo 'export PATH=$PATH:'$BASEDIR'/a5-miseq/bin' >> ~/.bashrc

echo "--------------------------------------------------------------------------------"
echo "Verifying a5-miseq setup"
./a5-miseq/test.a5.sh

#Installing GARM

echo 'export GARMBIN='$BASEDIR'/garm/bin' >> ~/.bashrc
echo 'export GARMLIB='$BASEDIR'/garm/lib' >> ~/.bashrc
echo 'export MUMBIN='$BASEDIR'/garm/MUMmer3.22' >> ~/.bashrc
echo 'export AMOSBIN='$BASEDIR'/garm/amos-3.0.0/bin' >> ~/.bashrc
echo 'export AMOSLIB='$BASEDIR'/garm/amos-3.0.0/lib' >> ~/.bashrc 
echo 'export PATH=$PATH:'$BASEDIR'/garm' >> ~/.bashrc	# add GARM.pl to PATH

source ~/.bashrc

echo "--------------------------------------------------------------------------------"
echo "Verifying GARM setup"
./garm/config.pl

#Installing Canu
cd canu/src
make -j >/dev/null

echo "--------------------------------------------------------------------------------"
echo "Verifying Canu setup"
cd ../Linux-amd64/
curl -L -o p6.25x.fastq http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq
./bin/canu -p ecoli -d ecoli-auto genomeSize=4.8m -pacbio-raw p6.25x.fastq
cd ..

echo "Check every stage log to verify setup errors!"
