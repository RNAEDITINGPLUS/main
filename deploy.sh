#!/usr/bin/env bash
HOMEDIR=$(cd `dirname $0`; pwd)

echo "============================================\n"
echo "Installing Apache2\n"
echo "============================================\n"
sudo apt-get install apache2
echo "============================================\n"
echo "Installing MySQL\n"
echo "============================================\n"
sudo apt install mysql-server mysql-client
echo "============================================\n"
echo "Installing PHP\n"
echo "============================================\n"
sudo apt-get install libapache2-mod-php5 php5 php5-gd php5-mysql
sudo apt-get install php-mbstring
sudo apt-get install php-gettext
echo "============================================\n"
echo "Enter database host:\n"
echo "============================================\n"
read dbh
sed -i "s?{DBH}?$dbh?g" $HOMEDIR/webapp/ThinkPHP/Conf/convention.php
echo "============================================\n"
echo "Enter database user:\n"
echo "============================================\n"
read dbu
sed -i "s?{DBU}?$dbu?g" $HOMEDIR/webapp/ThinkPHP/Conf/convention.php
echo "============================================\n"
echo "Enter database password:\n"
echo "============================================\n"
read dbpw
sed -i "s?{DBPW}?$dbpw?g" $HOMEDIR/webapp/ThinkPHP/Conf/convention.php
echo "============================================\n"
echo "Enter database name:\n"
echo "============================================\n"
read dbn
sed -i "s?{DBN}?$dbn?g" $HOMEDIR/webapp/ThinkPHP/Conf/convention.php
echo "============================================\n"
echo "Enter database port:\n"
echo "============================================\n"
read dbpt
sed -i "s?{DBPT}?$dbpt?g" $HOMEDIR/webapp/ThinkPHP/Conf/convention.php
echo "============================================\n"
echo "Creating tables\n"
echo "============================================\n"
mysql -h localhost -u $DBU -p$DBPW < database.sql
salt=$RANDOM
sed -i "s?{SALT}?$salt?g" $HOMEDIR/webapp/listener.php
sed -i "s?{SALT}?$salt?g" $HOMEDIR/node/sender.php
sed -i "s?{RUNPATH}?$HOMEDIR/node?g" $HOMEDIR/webapp/listener.php
echo "============================================\n"
echo "Successfully configured the frontend node\n"
echo "============================================\n"