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

import sqlite3
import getopt
import sys
from repConfig import getConfig
import os

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ht:s:a:b:f:e:", ["help", "type=","seq=", "url1=", "url2=", "fallback=", "expr="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit()
    flag = 0
    seq = None; url = None; fb = None; url1 = None; url2 = None; expr = None;
    if len(opts)==0:
        help()
        sys.exit()
    for o, a in opts:
        if o in ("-h","--help"):
            help()
            sys.exit()
        elif o in ("-t", "--type"):
            flag = int(a)
        elif o in ("-s", "--seq"):
            seq = a
        elif o in ("-a", "--url1"):
            url1 = a
        elif o in ("-b", "--url2"):
            url2 = a
        elif o in ("-f", "--fallback"):
            fb = a
        elif o in ("-e", "--expr"):
            expr = a
    if flag == None or (flag != 2 and fb == None) or (flag != 2 and url2 == None):
        sys.exit()
        
    if int(flag) == 1:
        if int(seq) == 2:
            url = url1+';'+url2
        else:
            url = url2
        res = put(seq, url, fb)
        print res
    elif int(flag) == 2:
        run()
    elif int(flag) == 3:
        if expr:
            res = put(flag, url2, fb, expr)
        else:
            res = put(flag, url2, fb)
        if res >= 1:
            print res
        else:
            print res        

def help():
    print "Queue implement"
    print "Options:"
    print "\t-h, --help:\n\t\t show this help message and exit"
    print "Mandatory:"
    print "\t-t, --type:\n\t\t 1 means add new task, 2 means a cron call"
    print "\t-s, --seq:\n\t\t Experiment type. If the data you provide is produced by paired-end RNA-seq, type 1; if the data is singled-end RNA-seq, type 2."
    print "\t-u, --url:\n\t\t The url to access your fastq file. If you have paired-end RNA-seq file, type RNA-seq F;RNA-seq R"
    print "\t-f, --fallback:\n\t\t If you would like to use hg19 as the reference genome, type 1, elsewise, type 0"

def conSQLite(db):
    if (db == 'memory'):
        con = sqlite3.connect(":memory:")
    else:
        try:
            con = sqlite3.connect(db)
        except Exception,e:
            print e
            return (None, None)
    cursor = con.cursor()
    return con, cursor


conn, cursor = conSQLite(getConfig("queue", "dbFile"))

if conn == None or cursor == None:
    import sys
    sys.exit(0)
    
#conn, cursor = conSQLite('queue.db')
def put(ty, job, fallback, ex=None):
    if ex:
        query = """INSERT INTO `queue` (`t`, `job`, `expr`, `status`, `fb`) VALUES (%s, '%s', '%s', 0, %s);""" % (ty, job, ex, fallback)
    else:
        query = """INSERT INTO `queue` (`t`, `job`, `status`, `fb`) VALUES (%s, '%s', 0, %s);""" % (ty, job, fallback)

    try:
        cursor.execute(query)
        id = cursor.lastrowid
        conn.commit()
    except Exception,e:
        print e
        return 0
    return id

def update(jid, field, value):
    query = """UPDATE `queue` SET `%s` = '%s' WHERE `id` = %s;""" % (field, value, jid)
    
    try:
        cursor.execute(query)
        conn.commit()
    except Exception, e:
        print e
        return 0
    return 1

def get():
    query = """SELECT `id`, `t`, `job`, `expr`, `fb` FROM `queue` WHERE `status` = 0 ORDER BY `id` LIMIT 1;"""
    cursor.execute(query)
    res = cursor.fetchone()
    if res != None:
        id, t, job, ex, fb = res
        update(id, 'status', 1)
        return id, t, job, ex, fb
    else:
        return 0


def check():
    #query = """SELECT COUNT(`job`) FROM `queue` WHERE `status` = 1;"""
    #cursor.execute(query)
    #runningTask = cursor.fetchone()
    #if runningTask != None:
    #    return runningTask[0]
    #else:
    query = """SELECT COUNT('job') FROM `queue` WHERE `status` > 0;"""
    cursor.execute(query)
    cTask = cursor.fetchone()
    if cTask != None and int(cTask[0]) <= 8:
        return 0
    else:
        return cTask[0]


def run():
    # if os.geteuid() != 0:
    # exit("You need to have root privileges to run this script.\nPlease try again, this time using 'sudo'. Exiting.")
    import dispatch
    loaded = check()
    jobStore = 0
    if loaded == 0:
        tmp = get()
        if tmp:
            jid, jt, job, ex, jfb = tmp
            jobStore = dispatch.repDispatch(jt, job, jfb, jid, ex)
            if jobStore != 0:
                if jt != 3:
                    update(jid, 'result', jobStore)
                update(jid, 'status', -1)
            else:
                update(jid, 'status', -1)

if __name__ == '__main__':
    main()
