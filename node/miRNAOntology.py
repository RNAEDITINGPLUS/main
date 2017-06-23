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

from database import mysqlConnection
from repConfig import getConfig
import os, re, time, sigEngine, repthermo, pysam
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import mirbase
import subprocess, pickle
from os import path
from pandas import DataFrame
import pandas as pd
from classifier import *
import mirmap
import mirmap.library_link

cnx, cursor = mysqlConnection()
scoreMatrix = {'T':1, 'G':2, 'R':4, 'C':8}
utrLibrary = {}
mirLibrary = {}
tgMapLibrary = {}
regulated = []
released = []
mirTable = pysam.Tabixfile(getConfig("datasets", "mirg"))

def isInSeedRegion(miRNAStart, editedPosition, strand):
    if strand == '+':
        if (editedPosition >= miRNAStart + 1 and editedPosition <= miRNAStart + 7):
            return True
        else:
            return False
    else:
        if (editedPosition <= miRNAStart - 1 and editedPosition >= miRNAStart - 7):
            return True
        else:
            return False

def convertRefGenomeSeedR():
    """Convert reference genome"""
    query = "SELECT * FROM miediting WHERE mirna LIKE '%p';"
    cursor.execute(query)
    
    for (rec_id, chromosome, position, gene, strand, mirna, seed_pos, annot1, annot2, alu, raw_seq, a2g, a2t, a2c) in cursor:
        queryCursor = cnx.cursor()
        query = "SELECT start, end FROM mirgenome WHERE chromosome='%s' AND strand='%s' AND start<=%d AND start>=%d" % (chromosome, strand, int(position) - 1, int(position) - 7)
        queryCursor.execute(query)
        try:
            data = queryCursor.fetchall()
            if (len(data)):
                query = "UPDATE miediting SET seed_pos=%d WHERE rec_id=%d" % (int(position) - int(data[0][0]) , rec_id)
                queryCursor.execute(query)
                cnx.commit()
        except  :
            print "Error: when executing "+query

def buildEditedSequenceOnMature(mirna, startPos, editPos, strand, matureSeq):
    """Build Edited mature miRNA(A>I)"""
    seqs = {}
    
    seqs['raw'] = matureSeq
    if strand == '+':
        relEditPos = editPos - startPos + 1
    else:
        relEditPos = startPos - editPos + 1
    
    if relEditPos == 0:
        seqs['a2g'] = "G" + seqs['raw'][1:]
        seqs['a2t'] = "T" + seqs['raw'][1:]
        seqs['a2c'] = "C" + seqs['raw'][1:]
        edited = "I" + seqs['raw'][1:]
    else:
        seqs['a2g'] = seqs['raw'][:relEditPos-1]+"G"+seqs['raw'][relEditPos:]
        seqs['a2t'] = seqs['raw'][:relEditPos-1]+"T"+seqs['raw'][relEditPos:]
        seqs['a2c'] = seqs['raw'][:relEditPos-1]+"C"+seqs['raw'][relEditPos:]
        edited = seqs['raw'][:relEditPos-1]+"I"+seqs['raw'][relEditPos:]

    records = []
    records.append(SeqRecord(Seq(seqs['a2g']), id=mirna+"{trick}I2G", description=""))
    records.append(SeqRecord(Seq(seqs['a2t']), id=mirna+"{trick}I2T", description=""))
    records.append(SeqRecord(Seq(seqs['raw']), id=mirna+"{trick}RAW", description=""))
    #print editPos, seqs
    return (records, edited, seqs)

def buildEditedSequenceOnMatureForMirMap(mirna, startPos, editPos, strand, matureSeq, ref='A', to='I'):
    """Build Edited mature miRNA(A>I)"""
    seqs = {}
    seqs['raw'] = matureSeq
    if strand == '+':
        relEditPos = editPos - startPos + 1
    else:
        relEditPos = startPos - editPos + 1
    if ref=='A' and to=='I':
        if relEditPos == 0:
            seqs['a2g'] = "G" + seqs['raw'][1:]
            seqs['a2t'] = "T" + seqs['raw'][1:]
            seqs['a2c'] = "C" + seqs['raw'][1:]
            edited = "I" + seqs['raw'][1:]
        else:
            seqs['a2g'] = seqs['raw'][:relEditPos-1]+"G"+seqs['raw'][relEditPos:]
            seqs['a2t'] = seqs['raw'][:relEditPos-1]+"T"+seqs['raw'][relEditPos:]
            seqs['a2c'] = seqs['raw'][:relEditPos-1]+"C"+seqs['raw'][relEditPos:]
            edited = seqs['raw'][:relEditPos-1]+"I"+seqs['raw'][relEditPos:]
    else:
        if relEditPos == 0:
            seqs['mut'] = to + seqs['raw'][1:]
        else:
            seqs['mut'] = seqs['raw'][:relEditPos-1]+to+seqs['raw'][relEditPos:]
    return seqs

def callMirMap(mir, utr):
    ret = mirMapAnalysis(mir, utr)
    if ret == -1.0:
        return 0
    else:
        record = ret.pop()
        return ret, record

def editingAnalysis(mir, utr):
    score = {}
    tags = {}
    reports = {}
    score['rs'] = mirMapAnalysis(mir['raw'], utr)
    if score['rs'] == -1.0:
        tags['rs'] = [-1.0]
    else:
        reports['rs'] = score['rs'].pop()
        feats, labels = formatInputForLibSVM(score['rs'])
        tags['rs'] = predictBySVM(getConfig('datasets', 'mirml'), feats, labels)
    score['gs'] = mirMapAnalysis(mir['a2g'], utr)
    if score['gs'] == -1.0:
        tags['gs'] = [-1.0]
    else:
        reports['gs'] = score['gs'].pop()
        feats, labels = formatInputForLibSVM(score['gs'])
        tags['gs'] = predictBySVM(getConfig('datasets', 'mirml'), feats, labels)
    score['ts'] = mirMapAnalysis(mir['a2t'], utr)
    if score['ts'] == -1.0:
        tags['ts'] = [-1.0]
    else:
        reports['ts'] = score['ts'].pop()
        feats, labels = formatInputForLibSVM(score['ts'])
        tags['ts'] = predictBySVM(getConfig('datasets', 'mirml'), feats, labels)
    score['cs'] = mirMapAnalysis(mir['a2c'], utr)
    if score['cs'] == -1.0:
        tags['cs'] = [-1.0]
    else:
        reports['cs'] = score['cs'].pop()
        feats, labels = formatInputForLibSVM(score['cs'])
        tags['cs'] = predictBySVM(getConfig('datasets', 'mirml'), feats, labels)
    if tags['rs'][0] == 1.0:
        if tags['gs'][0] == -1.0 and tags['ts'][0] == -1.0 and tags['cs'] == -1.0:
            return 1, score['rs'] #inactive
        else:
            return 0
    elif tags['gs'][0] == 1.0 or tags['ts'][0] == 1.0 or tags['cs'] == 1.0:
        if tags['gs'][0] == 1.0 and tags['ts'][0] == 1.0:
            return 2, score['gs']
        elif tags['gs'][0] == 1.0:
            return 3, score['gs']
        elif tags['ts'][0] == 1.0:
            return 4, score['ts']
        elif tags['cs'][0] == 1.0:
            return 5, score['cs']
    else:
        return 0

def mirMapAnalysis(mir, utr):
    mim = mirmap.mm(utr, mir)
    mim.libs = mirmap.library_link.LibraryLink(getConfig('program', 'mirMapLib'))
    mim.find_potential_targets_with_seed(allowed_lengths = [6,7], allowed_gu_wobbles={6:1, 7:1}, allowed_mismatches={6:2, 7:2}, take_best=False)
    
    try:
        scores = [mim.dg_duplex, mim.dg_binding, mim.dg_duplex_seed, mim.dg_binding_seed, mim.report()]
    except Exception, e:
        return -1
    return scores

def buildEditedSequenceOnSeed(mirna, startPos, editPos, strand, matureSeq):
    """Build Edited mature miRNA(A>I)"""
    if strand == '+':
        relEditPos = editPos - startPos + 1
    else:
        relEditPos = startPos - editPos + 1
    if (relEditPos < 2 or relEditPos > 8):
        return 0, 0
    seqs = {}
    seqs['raw'] = matureSeq
    if relEditPos == 1:
        seqs['a2g'] = "G"+seqs['raw'][2:9]
        seqs['a2t'] = "T"+seqs['raw'][2:9]
        seqs['a2c'] = "C"+seqs['raw'][2:9]
        edited = "I"+seqs['raw'][2:9]
        seqs['raw'] = seqs['raw'][1:9]
    else:
        seqs['a2g'] = seqs['raw'][0:relEditPos-1]+"G"+seqs['raw'][relEditPos:9]
        seqs['a2t'] = seqs['raw'][0:relEditPos-1]+"T"+seqs['raw'][relEditPos:9]
        seqs['a2c'] = seqs['raw'][0:relEditPos-1]+"C"+seqs['raw'][relEditPos:9]
        edited = seqs['raw'][0:relEditPos-1]+"I"+seqs['raw'][relEditPos:9]
        seqs['raw'] = seqs['raw'][0:9]
    
    records = []
    records.append(mirna + '-g\t' + seqs['a2g'] + '\t9606')
    records.append(mirna + '-t\t' + seqs['a2t'] + '\t9606')
    records.append(mirna + '-c\t' + seqs['a2c'] + '\t9606')
    records.append(mirna + '-r\t' + seqs['raw'] + '\t9606')
    return (records, edited)

def getUTR(gene):
    global utrLibrary
    if utrLibrary.has_key(gene):
        return utrLibrary[gene]
    else:
        print 'No record for %s' % gene
        return 0

def getmiR(mir):
    global mirLibrary
    if mirLibrary.has_key(mir):
        return mirLibrary[mir]
    else:
        print 'No record for %s' % mir
        return 0

def factory(chr, pos, jobId):
    global utrLibrary
    global tgMapLibrary
    global regulated
    global released
    
    if len(utrLibrary) == 0:
        loadUTRSeq()
    if len(tgMapLibrary) == 0:
        loadTGMap()

    counter = 0
    for x in chr.index:
        mirRes = whichMiRNA(chr.ix[x], pos.ix[x])

        if mirRes != 0:
            name, start, seed, strand, mature = mirRes
        else:
            print 'Cannot find the miRNA'
            continue
        fastas, a2i, rawList = buildEditedSequenceOnMature(name, start, pos.ix[x], strand, mature)
	
        if seed == 1:
            new_t, old_t, com_t, dif_t, sig = basicSupportForSeed(name, start, pos.ix[x], chr.ix[x], strand, mature, rawList, jobId)
        else:
            new_t, old_t, com_t, dif_t, sig = targetsPredictBymiRandaA2I(fastas, name, pos.ix[x], chr.ix[x], jobId)
	
        if new_t != 0 or old_t != 0 or com_t != 0 or dif_t != 0:
            counter += 1
            recordMiRNAEdit(name, pos.ix[x], chr.ix[x], seed, a2i, old_t, new_t, com_t, sig, jobId)
	
    reg = list(set(regulated).difference(set(released)))
    regS = ', '.join(reg)
    rel = list(set(released).difference(set(regulated)))
    relS = ', '.join(rel)
    recordRegulation(jobId, regS, relS)
    return counter

def buildMature4TS(name, tax, seqs):
    matures = []
    matures.append(name+'\t'+tax+'\t'+name+'-I2G\t'+seqs['a2g']+'\n')
    matures.append(name+'\t'+tax+'\t'+name+'-I2C\t'+seqs['a2c']+'\n')
    matures.append(name+'\t'+tax+'\t'+name+'-I2T\t'+seqs['a2t']+'\n')
    matures.append(name+'\t'+tax+'\t'+name+'-RAW\t'+seqs['raw']+'\n')
    return matures

def basicSupportForSeed(mir, start, pos, chr, strand, mature, rawList, jobId):
    global utrLibrary
    global tgMapLibrary
    sig = sigEngine.generatemiRNASig(mir, pos, chr, jobId)
    outFile = path.join(getConfig("dispatch", "mpool"), "out_"+sig)
    matureList = buildMature4TS(mir, '9606', rawList)
    miRSeeds, sA2I = buildEditedSequenceOnSeed(mir, int(start), pos, strand, mature)
    if not miRSeeds:
        return 0, 0, 0, 0, 'err'
    fasta = buildEditedSequenceOnMatureForMirMap(mir, start, pos, strand, mature)

    al = targetsPredictByTargetScan(miRSeeds, matureList, jobId)
    white = getDiffFromTargetScanOutput(al)
    step2Run = white[2]
    newList = white[0]
    dieList = white[1]
    tiList = white[3]
    newTPair = white[4]
    dieTPair = white[5]
    stepTPair = white[6]
    tIssuePair = white[7]

    dyUTR = [utrLibrary[white] for white in newList if utrLibrary.has_key(white)]
    dyUTR.extend([utrLibrary[black] for black in dieList if utrLibrary.has_key(black)])

    result = []

    for key, utrSym in enumerate(newList):  #high confidence new targets
        utr = getUTR(utrSym)
        if utr != 0:
            try:
                score = editingAnalysis(fasta, utr[newTPair[key][0]-25:newTPair[key][1]+20].strip())
                if score != 0:
                    status = score[0]
                    if status == 1:
                        item = 'raw'
                    elif status == 2:
                        item = 'a2g & a2t'
                    elif status == 3:
                        item = 'a2g'
                    elif status == 4:
                        item = 'a2t'
                    else:
                        return -1
                    if status != 0 and status != 1:
                        score = score[1]
                        regulated.append(tgMapLibrary[utrSym])
                        line = mir+'\t'+sig+'\t'+str(jobId)+'\tNEW\t'+item+'\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+utr[newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                    else:
                        regulated.append(tgMapLibrary[utrSym])
                        line = mir+'\t'+sig+'\t'+str(jobId)+'\tNEW\traw\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+utr[newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                else:
                    regulated.append(tgMapLibrary[utrSym])
                    line = mir+'\t'+sig+'\t'+str(jobId)+'\tNEW\traw\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+utr[newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                    result.append(line)
            except Exception, e:
                print e
                regulated.append(tgMapLibrary[utrSym])
                line = mir+'\t'+sig+'\t'+str(jobId)+'\tNEW\traw\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+utr[newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                result.append(line)
    for key, utrSym in enumerate(dieList):  #high confidence die targets       
        utr = getUTR(utrSym)
        if utr != 0:
            try:
                score = editingAnalysis(fasta, utr[dieTPair[key][0]-25:dieTPair[key][1]+20].strip())
                if score != 0:
                    status = score[0]
                    if status == 1:
                        item = 'raw'
                    elif status == 2:
                        item = 'a2g & a2t'
                    elif status == 3:
                        item = 'a2g'
                    elif status == 4:
                        item = 'a2t'
                    else:
                        return -1
                    if status == 1:
                        score = score[1]
                        released.append(tgMapLibrary[utrSym])
                        line = mir+'\t'+sig+'\t'+str(jobId)+'\tDIE\t'+item+'\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+utr[dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                    else:
                        released.append(tgMapLibrary[utrSym])
                        line = mir+'\t'+sig+'\t'+str(jobId)+'\tDIE\traw\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t \t \t \t \t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+utr[dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                else:
                    released.append(tgMapLibrary[utrSym])
                    line = mir+'\t'+sig+'\t'+str(jobId)+'\tDIE\traw\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t \t \t \t \t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+utr[dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                    result.append(line)
            except:
                released.append(tgMapLibrary[utrSym])
                line = mir+'\t'+sig+'\t'+str(jobId)+'\tDIE\traw\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t \t \t \t \t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+utr[dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                result.append(line)
    for key, utrSym in enumerate(step2Run):
        utr = getUTR(utrSym)
        if utr != 0:
            try:
                score = editingAnalysis(fasta, utr[stepTPair[key][0]-10:stepTPair[key][1]+25].strip())
                if score != 0:
                    status = score[0]
                    if status == 1:
                        item = 'raw'
                    elif status == 2:
                        item = 'a2g & a2t'
                    elif status == 3:
                        item = 'a2g'
                    elif status == 4:
                        item = 'a2t'
                    else:
                        return -1                    
                    if status != 1.0 and status != 0:
                        score = score[1]
                        regulated.append(tgMapLibrary[utrSym])
                        line = mir+'\t'+sig+'\t'+str(jobId)+'\tNEW\t'+item+'-so\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+utr[newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                    elif status == 1.0:
                        score = score[1]
                        released.append(tgMapLibrary[utrSym])
                        line = mir+'\t'+sig+'\t'+str(jobId)+'\tDIE\t'+item+'-so\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+utr[dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                    else:
                        line = mir+'\t'+sig+'\t'+str(jobId)+'\tRAW\t'+item+'\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t \t \t \t \t'+str(stepTPair[key][0])+'\t'+str(stepTPair[key][1])+'\t'+utr[stepTPair[key][0]-10:stepTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                else:
                    line = mir+'\t'+sig+'\t'+str(jobId)+'\tRAW\traw\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t \t \t \t \t'+str(stepTPair[key][0])+'\t'+str(stepTPair[key][1])+'\t'+utr[stepTPair[key][0]-10:stepTPair[key][1]+25].strip()+'\n'
                    result.append(line)
            except:
                line = mir+'\t'+sig+'\t'+str(jobId)+'\tRAW\traw\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t \t \t \t \t'+str(stepTPair[key][0])+'\t'+str(stepTPair[key][1])+'\t'+utr[stepTPair[key][0]-10:stepTPair[key][1]+25].strip()+'\n'
                result.append(line)
    for key, utrSym in enumerate(tiList):
        utr = getUTR(utrSym)
        if utr != 0:
            try:
                score = editingAnalysis(fasta, utr[tIssuePair[key][0]-10:tIssuePair[key][1]+25].strip())
                if score != 0:
                    status = score[0]
                    if status == 1:
                        item = 'raw'
                    elif status == 2:
                        item = 'a2g & a2t'
                    elif status == 3:
                        item = 'a2g'
                    elif status == 4:
                        item = 'a2t'
                    else:
                        return -1
                    if status != 1.0 and status != 0:
                        score = score[1]
                        regulated.append(tgMapLibrary[utrSym])
                        line = mir+'\t'+sig+'\t'+str(jobId)+'\tNEW\t'+item+'-so\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(tIssuePair[key][0])+'\t'+str(tIssuePair[key][1])+'\t'+utr[tIssuePair[key][0]-10:tIssuePair[key][1]+25].strip()+'\n'
                        result.append(line)
                    elif status == 1.0:
                        score = score[1]
                        released.append(tgMapLibrary[utrSym])
                        line = mir+'\t'+sig+'\t'+str(jobId)+'\tDIE\t'+item+'-so\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+utr[dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                    else:
                        released.append(tgMapLibrary[utrSym])
                        line = mir+'\t'+sig+'\t'+str(jobId)+'\tDIE\ta2t-so\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t-\t-\t-\t-\t'+str(tIssuePair[key][0])+'\t'+str(tIssuePair[key][1])+'\t'+utr[tIssuePair[key][0]-10:tIssuePair[key][1]+25].strip()+'\n'
                        result.append(line)
                else:
                    released.append(tgMapLibrary[utrSym])
                    line = mir+'\t'+sig+'\t'+str(jobId)+'\tDIE\ta2t-so\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t-\t-\t-\t-\t'+str(tIssuePair[key][0])+'\t'+str(tIssuePair[key][1])+'\t'+utr[tIssuePair[key][0]-10:tIssuePair[key][1]+25].strip()+'\n'
                    result.append(line)
            except:
                released.append(tgMapLibrary[utrSym])
                line = mir+'\t'+sig+'\t'+str(jobId)+'\tDIE\ta2t\t'+tgMapLibrary[utrSym]+'\t'+utrSym+'\t-\t-\t-\t-\t'+str(tIssuePair[key][0])+'\t'+str(tIssuePair[key][1])+'\t'+utr[tIssuePair[key][0]-10:tIssuePair[key][1]+25].strip()+'\n'
                result.append(line)
    fileA = open(outFile, 'w')
    fileA.writelines(result)
    fileA.close()
    
    query = """LOAD DATA LOCAL INFILE '%s' INTO TABLE %s (`mirna`, `sig`, `job`, `tag`, `way`, `gene_symbol`, `transcript_id`, `dg_duplex`, `dg_binding`, `dg_duplex_seed`, `dg_binding_seed`, `utr_start`, `utr_end`, `utr3`);""" % (outFile, getConfig("datasets", "mirs"))
    cursor.execute(query)
    cnx.commit()
    #new_t, old_t, com_t, dif_t
    
    return len(newList), len(step2Run)+len(dieList), len(step2Run), len(dieList), sig

def loadUTRSeq():
    global utrLibrary
    tmp = file(getConfig('datasets', 'utrLibrary'), 'rb')
    utrLibrary = pickle.load(tmp)
    tmp.close()

def loadTGMap():
    global tgMapLibrary
    tmp = file(getConfig('datasets', 'tgLibrary'), 'rb')
    tgMapLibrary = pickle.load(tmp)
    tmp.close()    

def convertRefGenomeMatureR():
    query = "SELECT * FROM miediting WHERE mirna LIKE '%p';"
    cursor.execute(query)
    
    for (rec_id, chromosome, position, gene, strand, mirna, seed_pos, annot1, annot2, alu, raw_seq, a2g, a2t, a2c) in cursor:
        queryCursor = cnx.cursor()
        query = "SELECT start, end FROM mirgenome WHERE chromosome='%s' AND strand='%s' AND start<=%d AND start>=%d" % (chromosome, strand, int(position) - 1, int(position) - 7)
        queryCursor.execute(query)
        try:
            data = queryCursor.fetchall()
            if (len(data)):
                query = "UPDATE miediting SET seed_pos=%d WHERE rec_id=%d" % (int(position) - int(data[0][0]) , rec_id)
                queryCursor.execute(query)
                cnx.commit()
        except :
            print query

def parsemiRanda(file, sig, jobId):
    with open(file) as f:
        all = f.read()
        n = re.sub("{trick}", "\t"+sig+"\t"+str(jobId)+"\t", all)
    
    fw = open(file, "w")
    fw.write(n)
    fw.close()
    
def whichMiRNA(chr, pos):
    """Find which miRNA it is"""
    global mirTable
    miRNA = None
    for hit in mirTable.fetch(reference=chr, start=int(pos)-1, end=int(pos), parser=pysam.asGTF()):
        chr, source, type, start, end, score, strand, phase, attributes = hit
        start = int(start)
        end = int(end)
        attributes = attributes.split(';')
        miRNA = attributes[2].split('=')[1]
        seq = attributes[4].split('=')[1]
        if miRNA != None:
            if strand == '+':
                seed = 1 if isInSeedRegion(start, pos, strand) else 0
                return (miRNA, start, seed, strand, seq)
            else:
                seed = 1 if isInSeedRegion(end, pos, strand) else 0
                return (miRNA, end, seed, strand, seq)
    return 0

def miRanda2gene(file, allList=[]):
    test= []
    with open(file) as f:
        rows = f.readlines()
        i = 0
        for row in rows:
            ele = row.split()
            query = "SELECT * FROM utr3 WHERE infasta = '%s';" % ele[4]
            cursor.execute(query)
            geneSymbol = cursor.fetchone()
            if len(allList) > 0:
                if geneSymbol:
                    if geneSymbol[6] not in allList:
                        continue
            if (geneSymbol == None and ele[4] != "NA"):
                geneSymbol = "Unknown("+ele[4]+")"
                newRows.append(row.replace(ele[4], ele[4].replace("hg38_refGene_","")+"\t"+geneSymbol))
            elif ele[4] != "NA" :
                newRows.append(row.replace(ele[4], ele[4].replace("hg38_refGene_","")+"\t"+str(geneSymbol[6])))
            i = i + 1
    fw = open(file, "w")
    fw.writelines(newRows)
    fw.close()

def isUnexpressed(expressionTable, gene_symbol):
    cols = list(expressionTable.columns)
    record = expressionTable[expressionTable['gene_short_name'] == gene_symbol]
    if 'FPKM' in cols:
        value = list(record['FPKM'])
    elif 'TPM' in cols:
        value = list(record['TPM'])
    else:
        return 1
    if len(value) == 0:
        return 2
    else:
        if float(value[0]) == 0:
            return 1
        else:
            return 0

def diffSetsBymiRanda(file, sig):
    mut = []
    wt = []
    rows = []
    taggedWithWay = []
    file_handler = open(file, 'r')
    rows = file_handler.readlines()
    '''
    with open(file) as f:
        rows = f.readlines()
        for row in rows:
            row = row.split()
            if row[3] == "RAW":
                wt.append(row[4])
            else:
                mut.append(row[4])
    '''
    for row in rows:
        row = row.split()
        if row[3] == "RAW":
            wt.append(row[4])
        else:
            mut.append(row[4])
    newT = set(mut).difference(set(wt))
    dieT = set(wt).difference(set(mut))
    comT = set(wt).intersection(set(mut))
    
    for row in rows:
        items = row.split('\t')
        if items[4] in newT:
            regulated.append(items[4])
            row = row.replace(items[4], 'NEW\t'+items[4])
        elif items[4] in dieT:
            released.append(items[4])
            row = row.replace(items[4], 'DIE\t'+items[4])
        else:
            row = row.replace(items[4], 'RAW\t'+items[4])
        taggedWithWay.append(row)
    fp = open(file, mode='w')
    fp.writelines(taggedWithWay)
    fp.close()
    #mut = list(set(mut))
    #wt = list(set(wt))
    #a.difference(b) --> in a not in b
    #diff = list(set(mut).difference(set(wt)))
    #union = list(set(mut).intersection(set(wt)))

    #for di in diff:
    #    query = """INSERT INTO `miranda_diff` (`sig`, `gene_symbol`) VALUES ("%s", "%s");""" % (sig, di)
    #    cursor.execute(query)
    #    cnx.commit()

    return (len(mut), len(wt), len(comT), len(dieT))

def uniqueTopmiRanda(miRandaOP, newT=[], dieT=[]):
    file = open(miRandaOP, 'r')
    rawSet = {}
    a2gSet = {}
    a2tSet = {}
    rawIndex = {}
    a2gIndex = {}
    a2tIndex = {}
    uniqueFile = file.readlines()
    for k, line in enumerate(uniqueFile):
        items = line.split('\t')
        if items[3] == 'RAW':
            if rawSet.has_key(items[4]):
                if rawSet[items[4]] > float(items[5]):
                    continue
            rawSet[items[4]] = float(items[5])
            rawIndex[items[4]] = k
        elif items[3] == 'I2G':
            if a2gSet.has_key(items[4]):
                if a2gSet[items[4]] > float(items[5]):
                    continue
            a2gSet[items[4]] = float(items[5])
            a2gIndex[items[4]] = k
        elif items[3] == 'I2T':
            if a2tSet.has_key(items[4]):
                if a2tSet[items[4]] > float(items[5]):
                    continue
            a2tSet[items[4]] = float(items[5])
            a2tIndex[items[4]] = k
    tmp = rawIndex.values()
    tmp.extend(a2gIndex.values())
    tmp.extend(a2tIndex.values())
    tmp.sort()
    newFile = open(miRandaOP, 'w')
    for lineNum in tmp:
        items = uniqueFile[lineNum].split('\t')
        if len(dieT) > 0:
            if items[4] in dieT and items[3] == 'RAW':
                continue
        newFile.write(uniqueFile[lineNum])
    newFile.close()


def targetsPredictBymiRandaA2I4ML(seqs, utrs, mirnaName, position, chromosome, jobId, inactive=[], newt=[]):
    #Set output files
    tag = str(time.time())
    miFile = path.join(getConfig("dispatch", "mpool"), "mir_"+tag+".fasta")
    utrFile = path.join(getConfig("dispatch", "mpool"), "ulib_"+tag+".fasta")
    outFile = path.join(getConfig("dispatch", "mpool"), "out_"+tag)

    SeqIO.write(seqs, miFile, "fasta")
    SeqIO.write(utrs, utrFile, "fasta")
    #Run miRanda
    os.system(os.path.join(getConfig("program", "miranda")+" "+miFile+" "+utrFile+" -out "+outFile+" -sc 60"))
    #Generate a signature for this editing event
    sig = sigEngine.generatemiRNASig(mirnaName, position, chromosome, jobId)
    parsemiRanda(outFile, sig, jobId)
    
    res = []
    
    uniqueTopmiRanda(outFile, newt, inactive)
    
    new_t, old_t, com_t, dif_t = diffSetsBymiRanda(outFile, sig)
    query = """LOAD DATA LOCAL INFILE '%s' INTO TABLE %s (`mirna`, `sig`, `job`, `tag`, `gene_symbol`, `score`, `energy`, `mi_start`, `mi_end`, `utr_start`, `utr_end`, `match_len`, `identity`, `similarity`, `mir_seq`, `lines`, `utr_seq`);""" % (outFile, getConfig("datasets", "mirtargets"))
    cursor.execute(query)
    cnx.commit()
    #delete(outFile)
    return new_t, old_t, com_t, dif_t, sig

def targetsPredictBymiRandaA2I(seqs, mirnaName, position, chromosome, jobId, whitelist=[]):
    #Set output files
    tag = str(time.time())
    miFile = path.join(getConfig("dispatch", "mpool"), "mir_"+tag+".fasta")
    outFile = path.join(getConfig("dispatch", "mpool"), "out_"+tag)

    SeqIO.write(seqs, miFile, "fasta")

    #Run miRanda
    os.system(os.path.join(getConfig("program", "miranda")+" "+miFile+" "+getConfig("datasets", "utr3")+" -out "+outFile))
    #Generate a signature for this editing event
    sig = sigEngine.generatemiRNASig(mirnaName, position, chromosome, jobId)
    parsemiRanda(outFile, sig, jobId)
    uniqueTopmiRanda(outFile)
    
    try:
        (new_t, old_t, com_t, dif_t) = diffSetsBymiRanda(outFile, sig)
    except :
        return 0,0,0,0
    
    query = """LOAD DATA LOCAL INFILE '%s' INTO TABLE %s (`mirna`, `sig`, `job`, `tag`, `way`, `gene_symbol`, `score`, `energy`, `mi_start`, `mi_end`, `utr_start`, `utr_end`, `match_len`, `identity`, `similarity`, `mir_seq`, `lines`, `utr_seq`);""" % (outFile, getConfig("datasets", "mirtargets"))
    cursor.execute(query)
    cnx.commit()

    #clean
    delete(miFile)
    delete(outFile)
    
    return (new_t, old_t, com_t, dif_t, sig)

def delete(filepath):
    try:
        from os import remove, path
    except:
        return 1
    if path.isfile(filepath):
        remove(filepath)
    
    return 0

def recordMiRNAEdit(mirName, pos, chr, seed, seq, ot, nt, ct, sig, jobId):
    try:
        query = """INSERT INTO `mirediting` (`mirna`, `accession`, `edit_pos_raw`, `edit_pos_chr`, `event`, `role`, `sequence`, `old_t`, `new_t`, `com_t`, `sig`, `job`) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', %s);""" % (mirName, mirbase.getmiRnaAcc(mirName), pos, chr, "A>I", seed, seq, ot, nt, ct, sig, jobId)
        cursor.execute(query)
        cnx.commit()
    except Exception, e:
        #print e
        return 0

def recordRegulation(jobId, reg, rel):
    try:
        query = """INSERT INTO `mir_overview` (`job_id`, `reg`, `rel`) VALUES ('%s', '%s', '%s');""" % (jobId, reg, rel)
        cursor.execute(query)
        cnx.commit()
    except Exception, e:
        #print e
        return 0

def getDiffFromTargetScanOutput(fileName):
    global tgMapLibrary
    '''
    Score distribution:
    1/2/3/8/9/10/11   --> new target
    4                 --> inactive target
    5/6/7/12/13/14/15 --> no change
    ===============================
    2016/06/13 Update
    Score distribution:
    1/2/3/8/9/10/11   --> new target
    4/5               --> inactive target
    6/7/12/13/14/15   --> no change
    '''
    tsTable = pd.read_table(fileName)
    genes = tsTable.groupby(tsTable['a_Gene_ID'])
    genemes = []
    newT = []
    dieT = []
    step2 = []
    tIssue = []
    pairInfoRaw = {}
    pairInfoA2G = {}
    pairInfoA2T = {}
    pairInfoA2C = {}
    pairInfoTIssue = {}
    newTPair = []
    dieTPair = []
    comTPair = []
    tIssuePair = []
    
    for name, group in genes:
        if tgMapLibrary.has_key(name):
            geneSymbol = tgMapLibrary[name]
        else:
            continue
        #genemes.append(name)
        tmp = tsTable[tsTable['a_Gene_ID'] == name]
        score = 0
        mirID = tmp['miRNA_family_ID']
        for mir in set(mirID):
            item = tmp[tmp['miRNA_family_ID'] == mir]
            if mir.find("-g") != -1:
                score += scoreMatrix['G']
                pairInfoA2G[name] = (int(list(item['UTR_start'])[0]), int(list(item['UTR_end'])[0]))
            elif mir.find("-a") != -1:
                #score += scoreMatrix['A']
                pairInfoA2A[name] = (item['UTR_start'], item['UTR_end'])
            elif mir.find("-t") != -1:
                score += scoreMatrix['T']
                pairInfoA2T[name] = (int(list(item['UTR_start'])[0]), int(list(item['UTR_end'])[0]))
            elif mir.find("-r") != -1:
                score += scoreMatrix['R']
                pairInfoRaw[name] = (int(list(item['UTR_start'])[0]), int(list(item['UTR_end'])[0]))
            elif mir.find("-c") != -1:
                score += scoreMatrix['C']
                pairInfoA2C[name] = (int(list(item['UTR_start'])[0]), int(list(item['UTR_end'])[0]))
        if score in (1,2,3,8,9,10,11):
            if geneSymbol not in newT:
                if pairInfoA2G.has_key(name):
                    newTPair.append(pairInfoA2G[name])
                    newT.append(name)
                elif pairInfoA2T.has_key(name):
                    newTPair.append(pairInfoA2T[name])
                    newT.append(name)
                elif pairInfoA2C.has_key(name):
                    newTPair.append(pairInfoA2C[name])
                    newT.append(name)
        elif score == 4:
            if geneSymbol not in dieT:
                if pairInfoRaw.has_key(name):
                    dieTPair.append(pairInfoRaw[name])
                    dieT.append(name)
        elif score == 5:
            if geneSymbol not in tIssue:
                if pairInfoA2T.has_key(name):
                    tIssuePair.append(pairInfoA2T[name])
                    tIssue.append(name)
        else:
            if geneSymbol not in step2:
                if pairInfoRaw.has_key(name):
                    comTPair.append(pairInfoRaw[name])
                    step2.append(name)
    return newT, dieT, step2, tIssue, newTPair, dieTPair, comTPair, tIssuePair


def targetsPredictByTargetScan(seqs, matures, jobId, taxonomy=9606):
    tag = str(time.time())
    miFile = path.join(getConfig("dispatch", "tspool"), "mir_"+tag+".fasta")
    mimFile = path.join(getConfig("dispatch", "tspool"), "mirm_"+tag+".fasta")
    outFile = path.join(getConfig("dispatch", "tspool"), "out_"+tag)
    f = open(miFile, "w")
    lines = [l+'\n' for l in seqs]
    f.writelines(lines)
    f.close()
    del(lines)
    f = open(mimFile, "w")
    f.writelines(matures)
    f.close()
    os.system(getConfig("program", "targetScan") + " " +miFile+ " " + getConfig("datasets", "ts_utr") + " " + outFile)
    return outFile
    
def perdiect(table = 'miediting'):
    query = "SELECT mirna FROM `%s` WHERE seed_pos != 0;" % table
    cursor.execute(query)    
    miRNAInfo = cursor.fetchall()
    mirs = []
    for mirna in miRNAInfo:
        targetsPredictBymiRandaA2I(str(mirna[0]))

    cursor.close()
    cnx.close()

def mirna2fasta(mirna, out = 'mirna_all.fasta'):
    result = []
    if mirna == 'all':
        query = "SELECT mature1_id, mature1_seq, mature1_acc FROM mirna;"
        cursor.execute(query)
        tmp = cursor.fetchall()
        query = "SELECT mature2_id, mature2_seq, mature2_acc FROM mirna;"
        cursor.execute(query)
        tmp.append (cursor.fetchall())
        for x in tmp:
            result.append(SeqRecord(Seq(str(x[1])), 
                                 id=str(x[0]),
                                 description=str(x[2])))

    SeqIO.write(result, out, "fasta")

if __name__ == '__main__':
    factory('chr14', 101040790, 123)