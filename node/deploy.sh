#!/usr/bin/env bash
HOMEDIR=$(cd `dirname $0`; pwd)

echo "============================================\n"
echo "Enter database host:\n"
echo "============================================\n"
read dbh
sed -i "s?{DBH}?$dbh?g" $HOMEDIR/node/config.conf
echo "============================================\n"
echo "Enter database user:\n"
echo "============================================\n"
read dbu
sed -i "s?{DBU}?$dbu?g" $HOMEDIR/node/config.conf
echo "============================================\n"
echo "Enter database password:\n"
echo "============================================\n"
read dbpw
sed -i "s?{DBPW}?$dbpw?g" $HOMEDIR/node/config.conf
echo "============================================\n"
echo "Enter database name:\n"
echo "============================================\n"
read dbn
sed -i "s?{DBN}?$dbn?g" $HOMEDIR/node/config.conf
echo "============================================\n"
echo "Enter database port:\n"
echo "============================================\n"
read dbpt
sed -i "s?{DBPT}?$dbpt?g" $HOMEDIR/node/config.conf
echo "============================================\n"
echo "Downloading reference files (hg38).\n"
echo "============================================\n"
wget https://www.rnaeditplus.org/Public/ref/reference.tar.gz
echo "============================================\n"
echo "Decomposing\n"
echo "============================================\n"
tar -zxvf reference.tar.gz
sed -i "s?{REFPATH}?$HOMEDIR/ref?g" $HOMEDIR/node/config.conf
sed -i "s?{RUNPATH}?$HOMEDIR?g" $HOMEDIR/node/config.conf
echo "============================================\n"
echo "Installing system packages.\n"
echo "============================================\n"
sudo apt-get install libmysqld-dev
sudo apt-get install python-dev
echo "============================================\n"
echo "Installing python packages.\n"
echo "============================================\n"
pip install -r prerequisite.txt
pip install fisher
cd reditools-1.0.4
sudo python setup.py install
echo "============================================\n"
echo "Compiling libsvm packages.\n"
echo "============================================\n"
cd libsvm-3.22/python
make
cd ..
mv libsvm.so.2 ..
echo "============================================\n"
echo "Successfully configured the computing node\n"
echo "============================================\n"