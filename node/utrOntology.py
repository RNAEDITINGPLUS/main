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
from dblite import conSQLite
from repConfig import getConfig
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os, pickle, re, pysam, mirmap, subprocess
from os import path
from sigEngine import generateUTR3Sig, anySig
from pandas import DataFrame
import pandas as pd
from mirbase import getmiRnaSymbol
import mirmap.library_link
from classifier import *

con, cursor = mysqlConnection()
scoreMatrix = {'T':1, 'G':2, 'R':4, 'C':8}
mirLibrary = {}
pf = open(getConfig("datasets", "utrComp"), mode='rb')
utr3Dict = pickle.load(pf)
mirTable = pysam.Tabixfile(getConfig("datasets", "utrbed"))

createTableSQL = """
CREATE TABLE [utr3_prediction] (
  [rec_id] INTEGER PRIMARY KEY AUTOINCREMENT, 
  [mirna] CHAR(255), 
  [fastatag] CHAR(255), 
  [repeat_chk] CHAR(255), 
  [sig] CHAR(255), 
  [tag] CHAR(255), 
  [gene_symbol] CHAR(255), 
  [score] DECIMAL(5, 2), 
  [energy] DECIMAL(5, 2), 
  [mi_start] INTEGER, 
  [mi_end] INTEGER, 
  [utr_start] INTEGER, 
  [utr_end] INTEGER, 
  [match_len] INTEGER, 
  [identity] DECIMAL(5, 2), 
  [similarity] INTEGER(5, 2), 
  [mir_seq] CHAR(255), 
  [lines] CHAR(255), 
  [utr_seq] CHAR(255));
"""
createIndex = """CREATE INDEX [rc] ON [utr3_prediction] ([repeat_chk]);"""

def buildEditedUTR3(chr, pos, event, out):
    '''build editied 3utr file'''
    
    query = """SELECT sequence, refseq_id, gene_symbol, start, end, strand FROM %s WHERE chromosome LIKE '%s%%' AND start <= %s AND end >= %s;""" % (getConfig("datasets", "refutr3"), chr, pos, pos)
    cursor.execute(query)
    raw = cursor.fetchall()
    utrs = []

    for utr in raw:
        if (utr[5] == '-'):
            editPos = utr[3] - int(pos)
            if event == "A>I":
                #utrs.append(SeqRecord(Seq(str(utr[0][:editPos])+"C"+str(utr[0][editPos+1:])), id=str(utr[1]+"{trick}I2C{tab}"+str(utr[2])), description=""))
                utrs.append(SeqRecord(Seq(str(utr[0][:editPos])+"G"+str(utr[0][editPos+1:])), id=str(utr[1]+"{trick}I2G{tab}"+str(utr[2])), description=""))
                utrs.append(SeqRecord(Seq(str(utr[0][:editPos])+"T"+str(utr[0][editPos+1:])), id=str(utr[1]+"{trick}I2T{tab}"+str(utr[2])), description=""))
                utrs.append(SeqRecord(Seq(str(utr[0])), id=str(utr[1]+"{trick}RAW{tab}"+str(utr[2])), description=str(utr[2])))
        else:
            editPos = int(pos) - utr[3]
            if event == "A>I":
                #utrs.append(SeqRecord(Seq(str(utr[0][:editPos])+"C"+str(utr[0][editPos+1:])), id=str(utr[1]+"{trick}I2C{tab}"+str(utr[2])), description=""))
                utrs.append(SeqRecord(Seq(str(utr[0][:editPos])+"G"+str(utr[0][editPos+1:])), id=str(utr[1]+"{trick}I2G{tab}"+str(utr[2])), description=""))
                utrs.append(SeqRecord(Seq(str(utr[0][:editPos])+"T"+str(utr[0][editPos+1:])), id=str(utr[1]+"{trick}I2T{tab}"+str(utr[2])), description=""))
                utrs.append(SeqRecord(Seq(str(utr[0])), id=str(utr[1]+"{trick}RAW{tab}"+str(utr[2])), description=""))
    
    SeqIO.write(utrs, out, "fasta")
    return out

def takeLongest(utrArr):
    maxValue = 0
    maxIndex = 0
    thisLength = 0

    for k, v in enumerate(utrArr):
        thisLength = v[4] - v[3]
        if thisLength > maxValue:
            maxValue = thisLength
            maxIndex = k
    return k

def buildEditedUTR3TS(chr, pos, event, out, gs):
    global mirTable
    utrs = []
    returnSeq = {'a2g': '', 'a2t': '', 'raw': ''}
    editPos = 0
    for hit in mirTable.fetch(reference=chr, start=int(pos), end=int(pos)+1, parser=pysam.asBed()):
        #chr, start, end, name, score, strand = hit
        chr, start, end, name, score, strand, seq = hit
        if strand == '-':
            editPos = int(end) - int(pos)
        else:
            editPos = int(pos) - int(start) - 1
        editedSeqG = seq[:editPos]+'G'+seq[editPos+1:]
        editedSeqT = seq[:editPos]+'T'+seq[editPos+1:]
        tag = name
        seq = seq.replace('\n', '')
        transcript = name.split('_')[0]
        allSeq = ''
        for transPart in utr3Dict[transcript]:
            tmp = transPart.split('-')
            for shit in mirTable.fetch(reference=chr, start=int(tmp[1]), end=int(tmp[2]), parser=pysam.asBed()):
                schr, sstart, send, sname, sscore, sstrand, sseq = shit
                if sstart == tmp[1] and send == tmp[2] and sname.find(transcript) != -1:
                    if sname == tag:
                        returnSeq['a2g'] += editedSeqG
                        returnSeq['a2t'] += editedSeqT
                        returnSeq['raw'] += sseq
                    else:
                        returnSeq['a2g'] += sseq
                        returnSeq['a2t'] += sseq
                        returnSeq['raw'] += sseq
        utrs.append(gs+'-g\t9606\t'+returnSeq['a2g']+'\n')
        utrs.append(gs+'-t\t9606\t'+returnSeq['a2t']+'\n')
        utrs.append(gs+'-r\t9606\t'+returnSeq['raw']+'\n')
    outFile = open(out, 'w')
    outFile.writelines(utrs)
    outFile.close()
    if len(utrs):
    	return out, returnSeq, editPos
    else:
    	return 0, 0, 0

def buildEditedUTR3TSOld(chr, pos, event, out, gs):
    '''build editied 3utr file'''
    
    query = """SELECT sequence, refseq_id, gene_symbol, start, end, strand FROM %s WHERE chromosome LIKE '%s%%' AND start <= %s AND end >= %s;""" % (getConfig("datasets", "refutr3"), chr, pos, pos)
    cursor.execute(query)
    raw = cursor.fetchall()
    utrs = []
    returnSeq = {}
    index = takeLongest(raw)
    utr = raw[index]
    #for utr in raw:
    utrS = str(utr[0]).strip().replace('T', 'U').replace('t', 'u')
    if (utr[5] == '-'):
        editPos = utr[3] - int(pos)
        if event == "A>I":
            utrs.append(str(utr[2])+'-g\t9606\t'+utrS[:editPos]+"G"+utrS[editPos+1:]+'\n')
            utrs.append(str(utr[2])+'-t\t9606\t'+utrS[:editPos]+"T"+utrS[editPos+1:]+'\n')
            utrs.append(str(utr[2])+'-r\t9606\t'+utrS+'\n')
            returnSeq['a2g'] = utrS[:editPos]+"G"+utrS[editPos+1:]
            returnSeq['a2t'] = utrS[:editPos]+"T"+utrS[editPos+1:]
            returnSeq['raw'] = utrS
    else:
        editPos = int(pos) - utr[3]
        if event == "A>I":
            utrs.append(str(utr[2])+'-g\t9606\t'+str(utrS[:editPos])+"G"+str(utrS[editPos+1:])+'\n')
            utrs.append(str(utr[2])+'-t\t9606\t'+str(utrS[:editPos])+"T"+str(utrS[editPos+1:])+'\n')
            utrs.append(str(utr[2])+'-r\t9606\t'+str(utrS)+'\n')
            returnSeq['a2g'] = str(utrS[:editPos])+"G"+str(utrS[editPos+1:])
            returnSeq['a2t'] = str(utrS[:editPos])+"T"+str(utrS[editPos+1:])
            returnSeq['raw'] = str(utrS)

    outFile = open(out, 'w')
    outFile.writelines(utrs)
    outFile.close()
    #SeqIO.write(utrs, out, "fasta")
    return out, returnSeq, editPos

def loadMir():
    global mirLibrary
    tmp = file(getConfig('datasets', 'mir_map'), 'rb')
    mirLibrary = pickle.load(tmp)
    tmp.close()

def miRNATargetsPredictionByTargetScan(utrFile, outFile):
    #print  os.system(os.path.join(getConfig("program", "miranda")+" "+getConfig("analyse", "miRNAFasta")+" "+utrFile+" -out "+outFile))
    #print getConfig("program", "miranda"), getConfig("analyse", "miRNAFasta"), utrFile, "-out", outFile
    os.system(getConfig("program", "targetScan") + " " +getConfig("datasets", "ts_mir")+ " " + utrFile + " " + outFile)
    return outFile

def miRNATargetsPredictionBymiRanda(utrFile, outFile):
    #print  os.system(os.path.join(getConfig("program", "miranda")+" "+getConfig("analyse", "miRNAFasta")+" "+utrFile+" -out "+outFile))
    #print getConfig("program", "miranda"), getConfig("analyse", "miRNAFasta"), utrFile, "-out", outFile 
    proc = subprocess.Popen([getConfig("program", "miranda"), getConfig("analyse", "miRNAFasta"), utrFile, "-out", outFile], shell=False, stdout = subprocess.PIPE, stderr = subprocess.STDOUT).wait()
    return outFile

def fileParse(file, sig, pos, chr, event = "A>I"):
    with open(file) as f:
        all = f.read()
        n = re.sub("{trick}", "\t"+sig+"\t", all)
        n = re.sub("{tab}", "\t", n)
    return n

def addRepeatChk(content):
    i = 0
    rows = content.split("\n")
    rows.remove('')
    for row in rows:
        ele = row.split("\t")
        repeatMarker = anySig(ele[0]+","+ele[1]+","+ele[2]+","+ele[5]+","+ele[6]+","+ele[7]+","+ele[8]+","+ele[9]+","+ele[10]+","+ele[11]+","+ele[12]+","+ele[13])
        
        try:
            rows[i] = row.replace(ele[1], ele[1].replace("hg38_refGene_","")+"\t"+repeatMarker)
        except:
            print ele
        i = i + 1
    return rows

def diffSetsBymiRanda(rows, sig):
    mut = []
    wt = []
    for row in rows:
        row = row.split()
        if row[4] == "RAW":
            wt.append(row[0])
        else:
            mut.append(row[0])
    mut = list(set(mut))
    wt = list(set(wt))
    diff = list(set(mut).difference(set(wt)))
    union = list(set(mut).intersection(set(wt)))

    for di in diff:
        query = """INSERT INTO `miranda_diff` (`sig`, `gene_symbol`) VALUES ("%s", "%s");""" % (sig, di)
        cursor.execute(query)
        con.commit()

    return (len(mut), len(wt), len(union), len(diff))

def findBestMatch(reports):
    maxValue = 0
    maxIndex = 0
    
    mirMapReport = reports.split('\t');

    if len(mirMapReport) == 0:
        return ''
    else:
        indexs = len(mirMapReport)/5
    
    for i in range(1, indexs+1):
        index = 3 + 5 * (i-1)
        thisV = mirMapReport[index].count('|')
        if maxValue < thisV:
            maxValue = thisV
            maxIndex = index
    indexPrime = index - 1
    return mirMapReport[indexPrime]

def removeUnexpressed(expressionTable, gene_symbol):
    cols = list(expressionTable.columns)
    #gs = str(gene_symbol).replace("-3p","").replace("-5p", "").replace("hsa-miR-", "MIR").replace("hsa-let-", "LET")
    gs = getmiRnaSymbol(con, cursor, (gene_symbol).replace("-3p","").replace("-5p", ""))
    record = expressionTable[expressionTable['gene_short_name'] == gs]
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

def editingAnalysis(mir, utr, start, end):
    score = {}
    tags = {}
    reports = {}
    score['rs'] = mirMapAnalysis(mir, utr['raw'][start:end].strip())
    u = ''
    
    if score['rs'] == -1.0:
        tags['rs'] = -1.0
    else:
        reports['rs'] = score['rs'].pop()
        u = findBestMatch(reports['rs'])
        feats, labels = formatInputForLibSVM(score['rs'])
        tags['rs'] = predictBySVM(getConfig('datasets', 'utrml'), feats, labels)
    score['gs'] = mirMapAnalysis(mir, utr['a2g'][start:end].strip())
    if score['gs'] == -1.0:
        tags['gs'] = -1.0
    else:
        reports['gs'] = score['gs'].pop()
        u = findBestMatch(reports['gs'])
        feats, labels = formatInputForLibSVM(score['gs'])
        tags['gs'] = predictBySVM(getConfig('datasets', 'utrml'), feats, labels)

    score['ts'] = mirMapAnalysis(mir, utr['a2t'][start:end].strip())
    if score['ts'] == -1.0:
        tags['ts'] = -1.0
    else:
        reports['ts'] = score['ts'].pop()
        u = findBestMatch(reports['ts'])
        feats, labels = formatInputForLibSVM(score['ts'])
        tags['ts'] = predictBySVM(getConfig('datasets', 'utrml'), feats, labels)

    if tags['rs'][0] == 1.0:
        if tags['gs'] == -1.0 and tags['ts'] == -1.0:
            return 1, score['rs'], u #inactive
        else:
            return 0

    elif tags['gs'][0] == 1.0 or tags['ts'][0] == 1.0:
        if tags['gs'][0] == 1.0 and tags['ts'][0] == 1.0:
            return 2, score['gs'], u
        elif tags['gs'][0] == 1.0:
            return 3, score['gs'], u
        elif tags['ts'][0] == 1.0:
            return 4, score['ts'], u
    else:
        return 0

def reverseUTR(utr):
    return utr[::-1] 

def mirMapAnalysis(mir, utr):
    mim = mirmap.mm(utr, mir)
    mim.libs = mirmap.library_link.LibraryLink(getConfig('program', 'mirMapLib'))
    tmp = mim.find_potential_targets_with_seed(allowed_lengths = [6,7], allowed_gu_wobbles={6:1, 7:1}, allowed_mismatches={6:2, 7:2}, take_best=False)
    try:
        scores = [mim.dg_duplex, mim.dg_binding, mim.dg_duplex_seed, mim.dg_binding_seed, mim.report()]
    except Exception, e:
        print e
        return -1
    return scores
    
def getDiffFromTargetScanOutput(fileName, startP, endP):
    global tgMapLibrary
    '''
    Score distribution:
    1/2/3/8/9/10/11 --> new target
    4     --> inactive target
    5/6/7/12/13/14/15 --> no change
    '''
    tsTable = pd.read_table(fileName)
    mirs = tsTable.groupby(tsTable['miRNA_family_ID'])
    genemes = []
    newT = []
    dieT = []
    step2 = []
    tIssue = []
    pairInfoRaw = {}
    pairInfoA2G = {}
    pairInfoA2T = {}
    newTPair = []
    dieTPair = []
    comTPair = []
    a2gDict = []
    a2tDict = []
    tIssuePair = []
    
    for name, group in mirs:
        tmp = tsTable[tsTable['miRNA_family_ID'] == name]
        score = 0
        geneID = tmp['a_Gene_ID']
        for gene in set(geneID):
            item = tmp[tmp['a_Gene_ID'] == gene]
            #print item
            if gene.find("-g") != -1:
                score += scoreMatrix['G']
                tmpStartArr = list(item['UTR_start'])
                tmpEndArr = list(item['UTR_end'])
                #if len(tmpStartArr) > 1:
                #    print 'aaa'
                flag = 0
                for k, v in enumerate(tmpStartArr):
                    if v >= startP:
                        if tmpEndArr[k] <= endP:
                            flag = 1
                            pairInfoA2G[name] = (int(list(item['UTR_start'])[k]), int(list(item['UTR_end'])[k]))
                            a2gDict.append(name)
                if flag == 0:
                    pairInfoA2G[name] = (int(list(item['UTR_start'])[0]), int(list(item['UTR_end'])[0]))
                    a2gDict.append(name)
            elif gene.find("-a") != -1:
                #score += scoreMatrix['A']
                pairInfoA2A[name] = (item['UTR_start'], item['UTR_end'])
            elif gene.find("-t") != -1:
                score += scoreMatrix['T']
                #pairInfoA2T[name] = (int(list(item['UTR_start'])[0]), int(list(item['UTR_end'])[0]))
                #a2tDict.append(name)
                tmpStartArr = list(item['UTR_start'])
                tmpEndArr = list(item['UTR_end'])
                flag = 0
                for k, v in enumerate(tmpStartArr):
                    if v >= startP:
                        if tmpEndArr[k] <= endP:
                            flag = 1
                            pairInfoA2T[name] = (int(list(item['UTR_start'])[k]), int(list(item['UTR_end'])[k]))
                            a2tDict.append(name)
                if flag == 0:
                    pairInfoA2T[name] = (int(list(item['UTR_start'])[0]), int(list(item['UTR_end'])[0]))
                    a2tDict.append(name)
            elif gene.find("-r") != -1:
                score += scoreMatrix['R']
                pairInfoRaw[name] = (int(list(item['UTR_start'])[0]), int(list(item['UTR_end'])[0]))
        geneSymbol = name.replace('-r', '').replace('-g', '').replace('-a', '').replace('-t', '')
        if score in (1,2,3,8,9,10,11):
            if geneSymbol not in newT:
                if pairInfoA2G.has_key(geneSymbol):
                    newTPair.append(pairInfoA2G[geneSymbol])
                    newT.append(geneSymbol)
                elif pairInfoA2T.has_key(geneSymbol):
                    newTPair.append(pairInfoA2T[geneSymbol])
                    newT.append(geneSymbol)
        elif score == 4:
            if geneSymbol not in dieT:
                if pairInfoRaw.has_key(geneSymbol):
                    dieTPair.append(pairInfoRaw[geneSymbol])
                    dieT.append(geneSymbol)
        elif score == 5:
            if geneSymbol not in tIssue:
                if pairInfoA2T.has_key(name):
                    tIssuePair.append(pairInfoA2T[name])
                    tIssue.append(name)        
        else:
            if geneSymbol not in step2:
                if pairInfoRaw.has_key(geneSymbol):
                    comTPair.append(pairInfoRaw[geneSymbol])
                    step2.append(geneSymbol)
    return newT, dieT, step2, tIssue, newTPair, dieTPair, comTPair, tIssuePair, a2gDict, a2tDict

def utrOntologymiRNA(gene, chr, pos, event, jid, fallback=0, expr='no'):
    global mirLibrary
    loadMir()
    sig = generateUTR3Sig(gene, pos, chr, jid, event)
    #utrf = buildEditedUTR3(chr, pos, event, path.join(getConfig("dispatch", "mpool"), "ue_"+sig+".fasta"))
    #pref = miRNATargetsPredictionBymiRanda(utrf, path.join(getConfig("dispatch", "mpool"), "up_"+sig+".fasta"))
    utrf, seqArr, relP = buildEditedUTR3TS(chr, pos, event, path.join(getConfig("dispatch", "tspool"), "ue_"+sig+".fasta"), gene)
    if utrf == 0:
    	print gene, chr, pos
    	return 0
    pref = miRNATargetsPredictionByTargetScan(utrf, path.join(getConfig("dispatch", "tspool"), "up_"+sig+".out"))
    #parsedFile = fileParse(pref, sig, pos, chr, event)
    startPos = relP - 25
    endPos = relP + 25
    white = getDiffFromTargetScanOutput(pref, startPos, endPos)
    step2Run = white[2]
    newList = white[0]
    dieList = white[1]
    tiList = white[3]
    newTPair = white[4]
    dieTPair = white[5]
    stepTPair = white[6]
    tIssuePair = white[7]
    a2gDicti = white[8]
    a2tDicti = white[9]
    new = 0; die = 0; raw = 0;
    outFile = path.join(getConfig("dispatch", "mpool"), "out_"+sig)
    
    result = []
    for key, mirSym in enumerate(newList):  #high confidence new targets
        try:
            if mirLibrary.has_key(mirSym):
                score = editingAnalysis(mirLibrary[mirSym], seqArr, newTPair[key][0]-18, newTPair[key][1]+18)
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
                    utr = seqArr[item]
                    if status != 0 and status != 1:
                        utrS = score[2]
                        score = score[1]
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tNEW\t'+item+'\t'+mirSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+utrS.strip()+'\n'
                        result.append(line)
                        new += 1
                    else:#SVM negative
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tNEW\t'+item+'\t'+mirSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+utr[newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                        new += 1
                else:#not supported by SVM
                    if mirSym in a2gDicti:
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tNEW\ta2g\t'+mirSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+seqArr['a2g'][newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                        new += 1
                    elif mirSym in a2tDicti:
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tNEW\ta2t\t'+mirSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+seqArr['a2t'][newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                        new += 1
            else:
                if mirSym in a2gDicti:
                    line = gene+'\t'+sig+'\t'+str(jid)+'\tNEW\ta2g\t'+mirSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+seqArr['a2g'][newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                    result.append(line)
                    new += 1
                elif mirSym in a2tDicti:
                    line = gene+'\t'+sig+'\t'+str(jid)+'\tNEW\ta2t\t'+mirSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+seqArr['a2t'][newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                    result.append(line)
                    new += 1
        except Exception, e:
            if mirSym in a2gDicti:
                line = gene+'\t'+sig+'\t'+str(jid)+'\tNEW\ta2g\t'+mirSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+seqArr['a2g'][newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                result.append(line)
                new += 1
            elif mirSym in a2tDicti:
                line = gene+'\t'+sig+'\t'+str(jid)+'\tNEW\ta2t\t'+mirSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+seqArr['a2t'][newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                result.append(line)
                new += 1
                
    for key, mirSym in enumerate(dieList):  #high confidence die targets 
        try:
            if mirLibrary.has_key(mirSym):
                score = editingAnalysis(mirLibrary[mirSym], seqArr, dieTPair[key][0]-18, dieTPair[key][1]+18)
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
                    utr = seqArr[item]
                    if status != 0 and status != 1:
                        utrS = score[2]
                        score = score[1]
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\t'+item+'\t'+mirSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+utrS.strip()+'\n'
                        result.append(line)
                        die += 1
                    else:#SVM negative
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\t'+item+'\t'+mirSym+'\t \t \t \t \t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+utr[dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                        die += 1
                else:#SVM negative
                    if mirSym in a2gDicti:
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\ta2g\t'+mirSym+'\t \t \t \t \t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+seqArr['a2g'][dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                        die += 1
                    elif mirSym in a2tDicti:
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\ta2t\t'+mirSym+'\t \t \t \t \t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+seqArr['a2t'][dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                        die += 1
                    else:
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\traw\t'+mirSym+'\t \t \t \t \t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+seqArr['a2t'][dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                        die += 1
            else:#SVM negative
                if mirSym in a2gDicti:
                    line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\ta2g\t'+mirSym+'\t \t \t \t \t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+seqArr['a2g'][dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                    result.append(line)
                    die += 1
                elif mirSym in a2tDicti:
                    line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\ta2t\t'+mirSym+'\t \t \t \t \t'+str(dieTPair[key][0])+'\t'+str(dieTPair[key][1])+'\t'+seqArr['a2t'][dieTPair[key][0]-10:dieTPair[key][1]+25].strip()+'\n'
                    result.append(line)
                    die += 1
        except Exception, e:
            if mirSym in a2gDicti:
                line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\ta2g\t'+mirSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+seqArr['a2g'][newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                result.append(line)
                die += 1
            elif mirSym in a2tDicti:
                line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\ta2t\t'+mirSym+'\t \t \t \t \t'+str(newTPair[key][0])+'\t'+str(newTPair[key][1])+'\t'+seqArr['a2t'][newTPair[key][0]-10:newTPair[key][1]+25].strip()+'\n'
                result.append(line)
                die += 1
    for key, mirSym in enumerate(step2Run):
        try:
            if mirLibrary.has_key(mirSym):
                score = editingAnalysis(mirLibrary[mirSym], seqArr, stepTPair[key][0]-18, stepTPair[key][1]+18)
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
                    utr = seqArr[item]
                    if status != 1.0 and status != 0:
                        utrS = score[2]
                        score = score[1]
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tNEW\t'+item+'-so\t'+mirSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(stepTPair[key][0])+'\t'+str(stepTPair[key][1])+'\t'+utrS.strip()+'\n'
                        result.append(line)
                        new += 1
                    elif status == 1.0:
                        score = score[1]
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\t'+item+'-so\t'+mirSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(stepTPair[key][0])+'\t'+str(stepTPair[key][1])+'\t'+utr[stepTPair[key][0]-10:stepTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                        die += 1
                    else:
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tRAW\t'+item+'\t'+mirSym+'\t \t \t \t \t'+str(stepTPair[key][0])+'\t'+str(stepTPair[key][1])+'\t'+utr[stepTPair[key][0]-10:stepTPair[key][1]+25].strip()+'\n'
                        result.append(line)
                        raw += 1
                else:
                    line = gene+'\t'+sig+'\t'+str(jid)+'\tRAW\traw\t'+mirSym+'\t \t \t \t \t'+str(stepTPair[key][0])+'\t'+str(stepTPair[key][1])+'\t'+seqArr['a2g'][stepTPair[key][0]-10:stepTPair[key][1]+25].strip()+'\n'
                    result.append(line)
                    raw += 1
            else:
                line = gene+'\t'+sig+'\t'+str(jid)+'\tRAW\traw\t'+mirSym+'\t \t \t \t \t'+str(stepTPair[key][0])+'\t'+str(stepTPair[key][1])+'\t'+seqArr['a2g'][stepTPair[key][0]-10:stepTPair[key][1]+25].strip()+'\n'
                result.append(line)
                raw += 1
        except Exception, e:
            line = gene+'\t'+sig+'\t'+str(jid)+'\tRAW\traw\t'+mirSym+'\t \t \t \t \t'+str(stepTPair[key][0])+'\t'+str(stepTPair[key][1])+'\t'+seqArr['a2g'][stepTPair[key][0]-10:stepTPair[key][1]+25].strip()+'\n'
            result.append(line)
            raw += 1
    for key, mirSym in enumerate(tiList):
        try:
            if mirLibrary.has_key(mirSym):
                score = editingAnalysis(mirLibrary[mirSym], seqArr, tIssuePair[key][0]-18, tIssuePair[key][1]+18)
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
                    utr = seqArr[item]
                    if status != 1.0 and status != 0:
                        utrS = score[2]
                        score = score[1]
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tNEW\t'+item+'-so\t'+mirSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(tIssuePair[key][0])+'\t'+str(tIssuePair[key][1])+'\t'+utrS.strip()+'\n'
                        result.append(line)
                        new += 1
                    elif status == 1.0:
                        score = score[1]
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\t'+item+'-so\t'+mirSym+'\t'+str(score[0])+'\t'+str(score[1])+'\t'+str(score[2])+'\t'+str(score[3])+'\t'+str(tIssuePair[key][0])+'\t'+str(tIssuePair[key][1])+'\t'+utr[tIssuePair[key][0]-10:tIssuePair[key][1]+25].strip()+'\n'
                        result.append(line)
                        die += 1
                    else:
                        line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\t'+item+'-so\t'+mirSym+'\t-\t-\t-\t-\t'+str(tIssuePair[key][0])+'\t'+str(tIssuePair[key][1])+'\t'+utr[tIssuePair[key][0]-10:tIssuePair[key][1]+25].strip()+'\n'
                        result.append(line)
                        raw += 1
                else:
                    line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\ta2t-so\t'+mirSym+'\t-\t-\t-\t-\t'+str(tIssuePair[key][0])+'\t'+str(tIssuePair[key][1])+'\t'+seqArr['a2g'][tIssuePair[key][0]-10:tIssuePair[key][1]+25].strip()+'\n'
                    result.append(line)
                    raw += 1
            else:
                line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\ta2t-so\t'+mirSym+'\t-\t-\t-\t-\t'+str(tIssuePair[key][0])+'\t'+str(tIssuePair[key][1])+'\t'+seqArr['a2g'][tIssuePair[key][0]-10:tIssuePair[key][1]+25].strip()+'\n'
                result.append(line)
                raw += 1
        except Exception, e:
            line = gene+'\t'+sig+'\t'+str(jid)+'\tDIE\ta2t-so\t'+mirSym+'\t-\t-\t-\t-\t'+str(tIssuePair[key][0])+'\t'+str(tIssuePair[key][1])+'\t'+seqArr['a2g'][tIssuePair[key][0]-10:tIssuePair[key][1]+25].strip()+'\n'
            result.append(line)
            raw += 1
    fileA = open(outFile, 'w')
    fileA.writelines(result)
    fileA.close()
    
    query = """LOAD DATA LOCAL INFILE '%s' INTO TABLE %s (`gene`, `sig`, `job`, `tag`, `way`, `mir`, `dg_duplex`, `dg_binding`, `dg_duplex_seed`, `dg_binding_seed`, `utr_start`, `utr_end`, `utr3`);""" % (outFile, getConfig("analyse", "utr3Targets"))
    cursor.execute(query)
    con.commit()
    #new_t, old_t, com_t, dif_t
    recordUTR3Edit(gene, pos, chr, die, new, raw, sig, jid)
    return new+die
    
def recordUTR3Edit(gene, pos, chr, ot, nt, ct, sig, jobId):
    try:
        query = """INSERT INTO `%s` (`gene`, `edit_pos_chr`, `edit_pos_raw`, `event`, `old_t`, `new_t`, `com_t`, `sig`, `job`) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');""" % (getConfig("datasets", "utrPri"), gene, chr, pos, "A>I", ot, nt, ct, sig, jobId)
        cursor.execute(query)
        con.commit()
    except Exception, e:
        print e
        return 0