#!/usr/bin/env python
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

import subprocess
from repConfig import getConfig
from os import path
from sigEngine import randSig
from multiprocessing import cpu_count
from sys import argv
from queue import update
from database import mysqlConnection
import cache


def ftpSize(hostName, path):
    from ftplib import FTP
    try:
        ftp = FTP(hostName)
    except:
        print 'FTP error!'
        return -1
    try:
        ftp.login()
        ftp.voidcmd('TYPE I')
        size = ftp.size(path)
        ftp.quit()
    except:
        print 'can not login anonymously or connect error or function exec error!'
        return 0
    return size


def httpSize(hostName, path):
    import httplib
    try:
        conn = httplib.HTTPConnection(hostName)
        conn.request("GET", path)
        resp = conn.getresponse()
    except:
        print 'connect error!'
        return 0
    return int(resp.getheader("content-length"))


def getRemoteSize(head, hostName, path):
    size = -1
    if head == 'ftp':
        size = ftpSize(hostName, path)
    if head == 'http' or head == 'https':
        size = httpSize(hostName, path)
    if head == 'loclf':
        if os.path.exists(path):
            size = os.path.getsize(path)
        else:
            size = 0
    return size


def getPath(url):
    import re
    regSite = re.compile(r"^(?:(?:http|https|ftp):\/\/)?([\w\-_]+(?:\.[\w\-_]+)+)", re.IGNORECASE)
    head = str(url).split(":")[0]
    try:
        port = str(url).split(':')[2].split('/')[0]
    except:
        port = ""
    matched = regSite.match(url)
    
    if regSite.match(url):
        hostName = regSite.match(url).group(1)
        if len(port)>1:
            hostName = hostName + ':' + port
        path = str(url).split(hostName)[1]
    else:
        head = "loclf"
        hostName = "local"
        path = url
    return (head, hostName, path)


def getRemoteSizeFactory(url):
    url = url.strip()
    if len(url) <= 4:
        exit(-1)
    urlArr = url.split(" ")
    totalSize = 0
    for url in urlArr:
        head, hostName, path = getPath(url)
        if len(hostName) <= 4 or len(path) < 1:
            exit(-1)
        totalSize += getRemoteSize(head, hostName, path)
    return totalSize


def wget(url, path):
    upperLimit = 5 * 1024 * 1024 * 1024
    if getRemoteSizeFactory(url) <= upperLimit:
        proc = subprocess.Popen(["wget", "-P", path, url, "-q"], shell=False, stdout=None,
                                stderr=subprocess.STDOUT).wait()
        return proc
    else:
        return 1


def notify(jobId, err):
    import urllib2
    url = 'http://localhost/dispatch/autorun/sender.php?yl=' + str(jobId) + '&er=' + str(err)
    try:
        response = urllib2.urlopen(url)
        html = response.read()
    except:
        print 'Fail to open url'
    return 0


def updateADAR(adar1, adar2, jid):
    con, cur = mysqlConnection()
    query = """UPDATE `jobs` SET `adar1` = '%s', `adar2` = '%s' WHERE `id` = %s;""" % (adar1, adar2, jid)
    try:
        cur.execute(query)
        con.commit()
        con.close()
        return 1
    except Exception, e:
        print e
        return 0


def repDispatch(t, remoteFile, fallback=0, jobId=0, expr=None):
    folderName = randSig()
    if int(fallback) == 1:
        refType = "fallback"
    else:
        refType = "dispatch"

    storeFolder = path.join(getConfig(refType, "mappingStore"), folderName)

    if int(t) == 1:  # singled-end RNA-seq
        if remoteFile.find(";") == -1:
            return 0
        else:
            tmp = remoteFile.split(";")
        dnaProc = wget(tmp[0], getConfig(refType, "remoteFile"))
        rnaProc = wget(tmp[1], getConfig(refType, "remoteFile"))
        update(jobId, 'status', 2)
        cmpProc = subprocess.Popen(
            [getConfig("program", "DnaRna"), "-i", path.join(getConfig(refType, "remoteFile"), path.basename(tmp[0])),
             "-j", path.join(getConfig(refType, "remoteFile"), path.basename(tmp[1])), "-f",
             getConfig(refType, "refgenome"), "-o", storeFolder, "-c", "10,1", "-Q", "33,64", "-q", "25,25", "-m",
             "20,20", "-s", "2", "-g", "1", "-u", "-a", "6-0", "-v", "2", "-n0.0", "-N0.0", "-V"], shell=False,
            stdout=None, stderr=subprocess.STDOUT).wait()
        update(jobId, 'status', 3)
        spProc = subprocess.Popen(
            [getConfig("program", "selectPos"), "-i", path.join(storeFolder, "outTable"), "-d", "12", "-c", "2", "-C",
             "10", "-v", "2", "-V", "0", "-f", "0.1", "-F", "1.0", "-e", "-u", "-o",
             path.join(storeFolder, "candidates.rep")], shell=False, stdout=None, stderr=subprocess.STDOUT).wait()
        annProc = subprocess.Popen(["AnnotateTable.py", "-a", getConfig(refType, "aluDatabase"), "-i",
                                    path.join(storeFolder, "candidates.rep"), "-u", "-c", "1,2,3", "-n", "RepMask",
                                    "-o", path.join(storeFolder, "candidates.rmsk.rep")], shell=False, stdout=None,
                                   stderr=subprocess.STDOUT).wait()

        # remove positions not annotated in SINE
        filterProc = subprocess.Popen(
            [getConfig("program", "filterTable"), "-i", path.join(storeFolder, "candidates.rmsk.rep"), "-f",
             getConfig(refType, "aluDatabase"), "-F", "SINE", "-E", "-o",
             path.join(storeFolder, "candidates.rmsk.alu.rep"), "-p"], shell=False, stdout=None,
            stderr=subprocess.STDOUT).wait()

        annProc = subprocess.Popen([getConfig("program", "annot"), "-a", getConfig(refType, "refGeneDatabase"), "-i",
                                    path.join(storeFolder, "candidates.rmsk.alu.rep"), "-u", "-c", "1,2", "-n",
                                    "RefSeq", "-o", path.join(storeFolder, "candidates.rmsk.alu.ann.rep")], shell=False,
                                   stdout=None, stderr=subprocess.STDOUT).wait()
        update(jobId, 'status', 4)
        strokes(remoteFile, jobId, refType)

        cufProc = subprocess.Popen(
            [getConfig("program", "cufflinks"), "-o", path.join(storeFolder, "expression.rep"), "-p",
             str(cpu_count() - 1), "-G", getConfig(refType, "refGeneDatabase"),
             path.join(getConfig(refType, "remoteFile"), path.basename(tmp[0]))], shell=False, stdout=None,
            stderr=subprocess.STDOUT).wait()

        if annProc == 0:
            return folderName
        else:
            return 0
    elif int(t) == 2:  # paired-end RNA-seq
        try:
            if remoteFile.find(";") == -1:
                notify(jobId, 3)
                return 0
            else:
                tmp = remoteFile.split(";")
            downloadProc1 = wget(tmp[0], storeFolder)
            if downloadProc1:
                notify(jobId, 3)
                return 0

            downloadProc2 = wget(tmp[1], storeFolder)
            if downloadProc2:
                notify(jobId, 3)
                return 0

            readsPath1 = path.join(storeFolder, path.basename(tmp[0]))
            readsPath2 = path.join(storeFolder, path.basename(tmp[1]))

            if downloadProc1 or downloadProc2:
                notify(jobId, 3)
                return 0

            if remoteFile.count(".gz") == 2:
                unzipProc1 = subprocess.Popen([getConfig("program", "gzip"), readsPath1], shell=False, stdout=None,
                                              stderr=subprocess.STDOUT).wait()
                readsPath1 = readsPath1.replace(".gz", "")
                unzipProc2 = subprocess.Popen([getConfig("program", "gzip"), readsPath2], shell=False, stdout=None,
                                              stderr=subprocess.STDOUT).wait()
                readsPath2 = readsPath2.replace(".gz", "")
            else:
                notify(jobId, 3)
                return 0

            if unzipProc1 or unzipProc2:
                notify(jobId, 3)
                return 0
        except Exception, e:
            print e
            notify(jobId, 3)
            return 0

        try:
            update(jobId, 'status', 2)
            # rnaSeq1 = path.join(getConfig(refType, "remoteFile"), path.basename(tmp[0]))
            # rnaSeq2 = path.join(getConfig(refType, "remoteFile"), path.basename(tmp[1]))
            mapProc = subprocess.Popen(
                [getConfig("program", "map"), "-q", "-x", getConfig("datasets", "hisat_index"), "-1", readsPath1, "-2",
                 readsPath2, "-S", path.join(storeFolder, getConfig("map", "outputName")), "-p", str(cpu_count() - 1),
                 "--dta-cufflinks"], shell=False, stdout=None, stderr=subprocess.STDOUT).wait()
            if mapProc:
                notify(jobId, 4)
                return 0
            update(jobId, 'status', 3)
            bamProc = subprocess.Popen([getConfig("program", "sam"), "view", "-@", str(cpu_count() - 1), "-bS",
                                        path.join(storeFolder, getConfig("map", "outputName")), "-o",
                                        path.join(storeFolder, getConfig("map", "bamName"))], shell=False, stdout=None,
                                       stderr=subprocess.STDOUT).wait()
            if bamProc:
                notify(jobId, 4)
                return 0
            update(jobId, 'status', 4)
            sortProc = subprocess.Popen(
                [getConfig("program", "sam"), "sort", path.join(storeFolder, getConfig("map", "bamName")),
                 path.join(storeFolder, getConfig("map", "sortedName"))], shell=False, stdout=None,
                stderr=subprocess.STDOUT).wait()
            if sortProc:
                notify(jobId, 4)
                return 0
            update(jobId, 'status', 5)
            expProc = subprocess.Popen(
                [getConfig("program", "expr"), "-p", "7", "-q", "-G", getConfig("map", "refGene"), "-o",
                 path.join(storeFolder, getConfig("map", "expName")),
                 path.join(storeFolder, getConfig("map", "sortedName")) + '.bam'], shell=False, stdout=None,
                stderr=subprocess.STDOUT).wait()
            expFile = open(path.join(storeFolder, getConfig("map", "expFile")), 'r')
            ea1 = ea2 = ''
            for line in expFile:
                items = line.split('\t')
                if items[4] == 'ADAR':
                    ea1 = items[9]
                elif items[4] == 'ADARB1':
                    ea2 = items[9]
            expFile.close()
            updateADAR(ea1, ea2, jobId)
            # if expProc:
            #    return 0
            # update(jobId, 'status', 6)
            cmpProc = subprocess.Popen(
                [getConfig("program", "Denovo"), "-i", path.join(storeFolder, getConfig("map", "sortedName")) + '.bam',
                 "-f", getConfig(refType, "refgenome"), "-o", storeFolder, "-t", str(cpu_count() - 1)], shell=False,
                stdout=None, stderr=subprocess.STDOUT).wait()
            update(jobId, 'status', 6)
            spProc = subprocess.Popen(
                [getConfig("program", "selectPos"), "-i", path.join(storeFolder, "outTable"), "-d", "-1", "-c", "2",
                 "-C", "0", "-v", "3", "-f", "0.1", "-e", "-o", path.join(storeFolder, "candidates.rep")], shell=False,
                stdout=None, stderr=subprocess.STDOUT).wait()
            update(jobId, 'status', 7)
            annProc = subprocess.Popen([getConfig("program", "annot"), "-a", getConfig(refType, "aluDatabase"), "-i",
                                        path.join(storeFolder, "candidates.rep"), "-u", "-c", "1,2,3", "-n", "RepMask",
                                        "-o", path.join(storeFolder, "candidates.rmsk.rep")], shell=False, stdout=None,
                                       stderr=subprocess.STDOUT).wait()
            update(jobId, 'status', 8)

            # remove positions not annotated in SINE
            # filterProc = subprocess.Popen(["FilterTable.py", "-i", path.join(storeFolder, "candidates.rmsk.rep"), "-f", getConfig(refType, "aluDatabase"), "-F", "SINE", "-E", "-o", path.join(storeFolder, "candidates.rmsk.alu.rep"), "-p"], shell=False, stdout = None, stderr = subprocess.STDOUT).wait()
            # remove snp
            # filterProc = subprocess.Popen([getConfig("program", "filterTable"), "-i", path.join(storeFolder, "candidates.rmsk.rep"), "-s", getConfig(refType, "dbsnp"), "-S", "snp", "-E", "-o", path.join(storeFolder, "candidates.rmsk.snp.rep"), "-p"], shell=False, stdout = None, stderr = subprocess.STDOUT).wait()
            # update(jobId, 'status', 9)
            annProc = subprocess.Popen(
                [getConfig("program", "annot"), "-a", getConfig(refType, "refGeneDatabase"), "-i",
                 path.join(storeFolder, "candidates.rmsk.rep"), "-u", "-c", "1,2", "-n", "RefSeq", "-o",
                 path.join(storeFolder, "candidates.rmsk.ann.rep")], shell=False, stdout=None,
                stderr=subprocess.STDOUT).wait()
            update(jobId, 'status', 9)
            # cufProc = subprocess.Popen([getConfig("program", "cufflinks"), "-o", path.join(storeFolder, "expression.rep"), "-p", str(cpu_count()-1), "-G", getConfig(refType, "refGeneDatabase"), path.join(getConfig(refType, "remoteFile"), path.basename(tmp[0]))], shell=False, stdout = None, stderr = subprocess.STDOUT).wait()
            filterProc = subprocess.Popen(
                ["python", getConfig("program", "filtertable"), "-s", getConfig("datasets", "snpbed"), "-i",
                 path.join(storeFolder, "candidates.rmsk.ann.rep"), "-o",
                 path.join(storeFolder, "candidates.rmsk.ann.rep.desnp")], shell=False, stdout=None,
                stderr=subprocess.STDOUT).wait()
            from main import save_candidates
            save_candidates(path.join(storeFolder, "candidates.rmsk.ann.rep.desnp"), jobId)
            cache.main(path.join(storeFolder, "candidates.rmsk.ann.rep.desnp"), jobId)
            update(jobId, 'status', 10)
            strokes(path.join(storeFolder, "candidates.rmsk.ann.rep.desnp.nd"), jobId, refType, expr)
            #strokes(path.join(storeFolder, "candidates.rmsk.ann.rep"), jobId, refType, expr)
            update(jobId, 'status', -1)
            annProc = 0
            if annProc == 0:
                return folderName
            else:
                return 0
        except Exception, e:
            print e
            notify(jobId, -1)
            return 0
    elif int(t) == 3:
        try:
            downProc = wget(str(remoteFile), getConfig(refType, "remoteFile"))
            localFile = path.join(getConfig(refType, "remoteFile"), path.basename(remoteFile))
        except Exception, e:
            print e
            notify(jobId, 3)
            return 0
        try:
            update(jobId, 'status', 4)
            from main import save_candidates
            save_candidates(localFile, jobId)
            ef_mir, ef_utr, ef_mis, ef_alt = cache.main(localFile, jobId)
            toAppendData = {
                'mir': ef_mir,
                'utr': ef_utr,
                'mis': ef_mis,
                'alt': ef_alt,
            }
            update(jobId, 'status', 10)
            strokes(localFile+'.nd', jobId, refType, expr, toAppendData)
                #strokes(remoteFile, jobId, refType, expr)
            update(jobId, 'status', -1)
        except Exception, e:
            print e
            notify(jobId, -1)
            return 0
    else:
        notify(jobId, 4)
        return 0


def strokes(job, jobId, refType, expr='no', cas=None):
    from main import analysis
    if path.exists(job):
        analysis(job, jobId, cas, refType, expr)
    else:
        downProc = wget(str(job), getConfig(refType, "remoteFile"))
        analysis(path.join(getConfig(refType, "remoteFile"), path.basename(job)), jobId, cas, refType, expr)

# if __name__ == "__main__":
#    if len(argv) == 4:
#        p = repDispatch(argv[1], argv[2], int(argv[3]))
#        print "Done!@hpc"
#    else:
#        print "Dispather for RNA Editing Plus"
#        print "Powered by Yao Li"
