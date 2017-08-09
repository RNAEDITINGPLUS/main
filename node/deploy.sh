#!/usr/bin/env bash
HOMEDIR=$(cd `dirname $0`; pwd)

echo "Making a dir to save reference files.\n"
mkdir ref
cd ref
echo "Initiating download.\n"
echo "============================================"
echo "Downloading reference files (hg38).\n"
echo "============================================"
wget https://www.rnaeditplus.org/Public/ref/reference.tar.gz
echo "============================================"
echo "Decomposing\n"
echo "============================================"
tar -zxvf reference.tar.gz
sed -i "s?{REFPATH}?$HOMEDIR/ref?g" $HOMEDIR/node/config.conf
cd ..
echo "============================================"
echo "Installing system packages.\n"
echo "============================================"
sudo apt-get install libmysqld-dev
sudo apt-get install python-dev
echo "============================================"
echo "Installing python packages.\n"
echo "============================================"
pip install -r prerequisite.txt
pip install fisher
tar -zxvf REDItools-1.0.4.tar.gz
cd reditools-1.0.4
sudo python setup.py install
wget http://mirmap.ezlab.org/downloads/miRmap-1.1.tar.gz
tar -zxvf miRmap-1.1.tar.gz
cd miRmap-1.1
