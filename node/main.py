#! /usr/bin/env python
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

from pandas import DataFrame
import pandas as pd
import subprocess
from repConfig import getConfig
from database import mysqlConnection

cnx, cursor = mysqlConnection()

def readTable(filePath):
    table = pd.read_table(filePath)
    return table

def readFPKM(filePath):
    expression = pd.read_table(filePath)
    return expression

def identifyUTR3(table, jobId, hg19=0, fpkm='no'):
    from utrOntology import utrOntologymiRNA
    counter = 0
    #editingOn3UTR = table[table['RefSeq_feat'] == '3UTR']
    editingOn3UTR = table[table['RefSeq_feat'].str.contains('3UTR')]
    result = {}
 
    for x in editingOn3UTR.index:
        event = editingOn3UTR.ix[x]
        chr = event['Region'] if event['Region'].find('chr') != -1 else 'chr'+event['Region']
        alt3UTR = utrOntologymiRNA(event['RefSeq_gid'], chr, event['Position'], 'A>I', jobId, fallback=0, expr=fpkm)
        
        if alt3UTR:
            counter += 1
    
    return counter

def identifymiRNA(table, jobId, hg19=0):
    import miRNAOntology
    symbols = list(table['RefSeq_gid'])

    num = len(symbols)
    mirs = []
    for i in range(num):
        if symbols[i].lower().find('mir') != -1:
            mirs.append(i)
    editingOnmiRNA = table.ix[mirs]
    
    return miRNAOntology.factory(editingOnmiRNA['Region'], editingOnmiRNA['Position'], jobId)

def identifymisSense(table, jobId, hg19=0):
    import missenseOntology
    counter = 0
    #editingOn3UTR = table[table['RefSeq_feat'] == '3UTR']
    condition = table[table['RefSeq_feat'].str.contains('CDS')]
    cdsi = []
    for x in condition.index:
        cdsi.append(x)
    editingOnCDS = table.ix[cdsi]
    
    return missenseOntology.factory(editingOnCDS['Region'], editingOnCDS['Position'], jobId)
    
def identifySplicing(table, jobId, hg19=0):

    import splicingOntology
    feats = list(table['RefSeq_feat'])

    num = len(feats)
    introns = []
    '''
    editingOnIntron = table[table['RefSeq_feat'] == 'intron']
    '''
    for i in range(num):
        if feats[i].lower().find('intron') == 0:
            introns.append(i)
    editingOnIntron = table.ix[introns]
    cols = list(table.columns)
    if 'Ref' in cols and 'Edited' in cols:
        return splicingOntology.factory(editingOnIntron['RefSeq_gid'], editingOnIntron['Region'], editingOnIntron['Position'], jobId, 1, editingOnIntron['Ref'], editingOnIntron['Edited'])
    else:
        return splicingOntology.factory(editingOnIntron['RefSeq_gid'], editingOnIntron['Region'], editingOnIntron['Position'], jobId)

def checkIntegrity(table, rawFile):
    cols = list(table.columns)
    if cols[0] == 'Region' or cols[1] == 'Position' or cols[2] == 'Strand':
        if 'RefSeq_gid' in cols and 'RefSeq_feat' in cols:
            return 0
        else:
            an1proc = subprocess.Popen([getConfig("program", "annot"), "-a", getConfig("dispatch", "refGeneDatabase"), "-i", rawFile, '-u', '-c', '1,2', '-n', 'RefSeq', '-o', rawFile+'.1ann'], shell=False, stdout = None, stderr = subprocess.STDOUT).wait()
            an2proc = subprocess.Popen([getConfig("program", "annot"), "-a", getConfig("dispatch", "aluDatabase"), "-i", rawFile+'.1ann', '-u', '-c', '1,2,3', '-n', 'RepMask', '-o', rawFile+'.ann'], shell=False, stdout = None, stderr = subprocess.STDOUT).wait()
            return 1 #needs to be annotated
    return 2 #error format

#def reformTable(table):

def updateJob(jid, utr, mir, spli, mis):
    from database import mysqlConnection
    n, c = mysqlConnection()
    query = """UPDATE `jobs` SET `utr`=%s, `mir`=%s, `splicing`=%s , `mis`=%s, `status`=5 WHERE `trace`=%s;""" % (utr, mir, spli, mis, jid)
    try:
        c.execute(query)
        n.commit()
    except Exception, e:
        print e
        return 0
    return 1

def reformTable(table, jobId):
    conservedCols = ('job', 'region', 'position', 'strand', 'frequency', 'repmask_gid', 'refseq_feat', 'refseq_gid')
    toDrop = []
    for i in table.columns:
        if i.lower() not in conservedCols:
            #del table[i]
            toDrop.append(i)
    reformedTable = table.drop(toDrop, axis=1)
    table['Job'] = jobId
    reformedTable['Job'] = jobId
    reformedTable.sort_index(axis=1)
    reformedTable.to_sql(getConfig("datasets", "candidate"), cnx, flavor='mysql', if_exists='append', index=False)

def notify(jobId):
    import urllib2
    url = 'http://localhost/dispatch/autorun/sender.php?yl='+str(jobId)
    try:
        response = urllib2.urlopen(url)
        html = response.read()
    except:
        print 'Fail to open url'
    return 0

def checkExpTable(exprTable):
    cols = list(exprTable.columns)
    if 'gene_short_name' in cols and ('FPKM' in cols or 'TPM' in cols):
        return 0
    else:
        return 1

def updateADAR(jobId, a1, a2):
    from database import mysqlConnection
    n, c = mysqlConnection()
    query = """UPDATE `jobs` SET `adar1`='%s', `adar2`='%s' WHERE  `trace`=%s;""" % (a1, a2, jobId)
    try:
        c.execute(query)
        n.commit()
        n.close()
    except Exception, e:
        print e
        return 0
    return 1

def getADAR(exprTable):
    adar1Pool = exprTable[exprTable['gene_short_name'] == 'ADAR']['FPKM']
    for t in adar1Pool:
        adar1 = str(t)
    adar2Pool = exprTable[exprTable['gene_short_name'] == 'ADARB1']['FPKM']
    for t in adar1Pool:
        adar2 = str(t)
    return adar1, adar2

def save_candidates(input_file, job_id):
    table = readTable(input_file)
    if checkIntegrity(table, input_file) != 0:
        table = readTable(input_file+".ann")
        if checkIntegrity(table, input_file+".ann") != 0:
            return 0
    reformTable(table, job_id)

def analysis(input, jobId, cachedNum=None, fallback=0, fpkm='no'):
    table = readTable(input)
    if checkIntegrity(table, input) != 0:
        table = readTable(input+".ann")
        if checkIntegrity(table, input+".ann") != 0:
            return 0
    # reformTable(table, jobId)
    
    if fpkm != 'no':
        expr = readTable(fpkm)
        if checkExpTable(expr) == 0:
                
            a1, a2 = getADAR(expr)
            updateADAR(jobId, a1, a2)

            u = identifyUTR3(table, jobId, fallback, expr)
            m = identifymiRNA(table, jobId, fallback)
        else:
            u = identifyUTR3(table, jobId, fallback)
            m = identifymiRNA(table, jobId, fallback)
    else:
        u = identifyUTR3(table, jobId, fallback)
        m = identifymiRNA(table, jobId, fallback)
    e = identifymisSense(table, jobId, fallback)
    s = identifySplicing(table, jobId, fallback)
    if cachedNum is not None:
        updateJob(jobId, u+cachedNum['utr'], m+cachedNum['mir'], s+cachedNum['alt'], e+cachedNum['mis'])
    else:
        updateJob(jobId, u, m, s, e)
    notify(jobId)
    #print "Job %s finished. %s %s %s" % (jobId, u, m, s)
    return 0

if __name__ == '__main__':
    analysis(r'cdst.txt', 1, fallback=0, fpkm='no')
