#!/usr/bin/python2.7
# Copyright (C) 2016 Li Yao <yaoli95@outlook.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

#import mysql.connector
import MySQLdb
from repConfig import getConfig


def module_exists(moduleName):
    try:
        __import__(moduleName)
    except ImportError:
        return False
    else:
        return True


def mysqlConnection(buffer=1):
    try:
        connection = MySQLdb.connect(host=getConfig("db", "host"), user=getConfig("db", "user"), passwd=getConfig("db", "password"), db=getConfig("db", "db_name"), port=int(getConfig("db", "port")), local_infile=1)
        #connection = MySQLdb.connect(user=getConfig("db", "user"), password=getConfig("db", "password"), database=getConfig("db", "db_name"), host=getConfig("db", "host"), port=int(getConfig("db", "port")))
    except Exception, e:
        print e
    if buffer == 1:
        try:
            cursor = connection.cursor()
        except Exception, e:
            print e
            return 1
        #cursor = connection.cursor(buffered=True)
    else:
        cursor = connection.cursor()
    return connection, cursor