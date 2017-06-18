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

from sigEngine import randSig, generateSSSig
import os, string, re, subprocess, pysam, time
from repConfig import getConfig
from random import Random
from database import mysqlConnection
from os import path

con, cursor = mysqlConnection()
filterThreshold = string.atof(getConfig("splicing", "mesSCT"))
mesVT = string.atof(getConfig("splicing", "mesVT"))
bpThreshold = string.atof(getConfig("splicing", "bpThreshold"))
bpVarThreshold = string.atof(getConfig("splicing", "bpVarThreshold"))
'''
Citation:
Eng L, Coutinho G, Nahas S, Yeo G, Tanouye R, Babaei M, Dork T, Burge C, Gatti RA. Nonclassical splicing mutations in the coding and noncoding regions of the ATM Gene: maximum entropy estimates of splice junction strengths. Hum Mutat. 2004 Jan; 23(1):67-76
FO Desmet, Hamroun D, Lalande M, Collod-Beroud G, Claustres M, Beroud C. Human Splicing Finder: an online bioinformatics tool to predict splicing signals. Nucleic Acid Research, 2009
'''
#splicingSite5Motif = re.compile("[CA]AGGT[AG]AG[ACGT]", re.IGNORECASE|re.DOTALL)
splicingSite5Motif = re.compile("[CA][ACGT][ACGT]GT[ACGT][ACGT][ACGT][ACGT]", re.IGNORECASE|re.DOTALL)
splicingSite3Motif = re.compile("[ATCG]{18}AG[AG][ATCG][ATCG]", re.IGNORECASE|re.DOTALL)
branchMotif = re.compile("[CGT][ACGT][CT][CT][CG]A[CGT]", re.IGNORECASE|re.DOTALL)

branchSiteMatris = {
    0:{'A': 0.0,    'C':4.25,   'G':2.62,   'T':2.72},
    1:{'A': 0.0,    'C':6.87,   'G':2.29,   'T':3.5},
    2:{'A': 0.0,    'C':25.72,  'G':0.0,    'T':1.88},
    3:{'A': 0.0,    'C':6.05,   'G':0.0,    'T':15.09},
    4:{'A': 0.0,    'C':11.82,  'G':6.89,   'T':0.0},
    5:{'A': 29.63,  'C':0.0,    'G':0.0,    'T':0.0},
    6:{'A': 0.0,    'C':6.62,   'G':3.89,   'T':2.72},
}
STATUSCODE = {
    'INACTIVE5SS': 1,
    'ENHENCED5SS': 2,
    'WEAKEND5SS':  3,
    'INACTIVE3SS':  4,
    'ENHENCED3SS':  5,
    'WEAKEND3SS':  6,
    'NEWBS':  7,
    'INACTIVEBS':  8,
    'WEAKENDBS':  9,
}
intronInfo = pysam.Tabixfile(getConfig("datasets", "intronInfo"))

def callMaxEntScan3ss(sequence, sig):
    #if len(sequence) != 23:
    #    print "The 3ss sequence must be 23 bases long [20 bases in the intron][3 bases in the exon]"
    #    return 0
    fileName = path.join(getConfig("dispatch", "etpool"), "MaxEntScan_" + sig)
    with open(fileName, 'w') as file:
        file.writelines(sequence)
    maxEntScan3 = subprocess.Popen(["perl", getConfig("program", "ss3"), fileName], shell=False, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    score = maxEntScan3.stdout.read()
    try:
        os.remove(fileName)
    except Exception, e:
        print e
    return score.split()

def callMaxEntScan5ss(sequence, sig):
    #if len(sequence) != 9:
    #    print "The 5ss sequence must be 9 bases long [3 bases in the exon][6 bases in the intron]"
    #    return 0
    fileName = path.join(getConfig("dispatch", "etpool"), "MaxEntScan_" + sig)
    with open(fileName, 'w') as file:
        file.writelines(sequence)
    maxEntScan5 = subprocess.Popen(["perl", getConfig("program", "ss5"), fileName], shell=False, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    score = maxEntScan5.stdout.read()
    os.remove(fileName)
    return score.split()

def isChangedSS(raw, edited, type):
    '''This function will return a tuple, when the first element equals to 1, the edited sequence can be seen as a 3' splicing site.
    If the second element equals to -1, it means the editing breaks the splice site. In the other case, if it returns 1, we consider that the editing creates a new splicing site'''
    runWell = 0
    changed = 0

    if type == 3:
        rawScore = callMaxEntScan3ss(raw)
        editedScore = callMaxEntScan3ss(edited)
    else:
        rawScore = callMaxEntScan5ss(raw)
        editedScore = callMaxEntScan5ss(edited)
    
    rs = string.atof(rawScore)
    es = string.atof(editedScore)

    if es > string.atof(getConfig("splicing", "mesSCT")):
        runWell = 1

    if ((es - rs) / rs) < -1 * float(getConfig("splicing", "mesVT")):
        changed = -1
    elif ((es - rs) / rs) > float(getConfig("splicing", "mesVT")):
        changed = 1

    return (runWell, changed, rs, es)

def getIntronInfo(geneSymbol, chr, position):
    global intronInfo
    """
    GET Intron Info
    Desc: get intron info by providing gene official symbol and mutation's position
    Return: If there is a hited record, a list will be returned.
    [rec_name, chromosome, start, end, sequence] 
    """
    result = []
    for iii in intronInfo.fetch(reference=chr, start=int(position), end=int(position)+1, parser=pysam.asBed()):
        chromosome, start, end, recName, geneSymbol, strand, sequence = iii
        result.append([recName, chromosome, int(start), int(end), strand, sequence.replace('\n', '')])
    '''
    query = """SELECT `rec_name`, `chromosome`, `start`, `end`, `strand` FROM %s WHERE `gene_symbol`='%s' AND `start`<=%s AND `END`>=%s;""" % (getConfig("datasets", "intronInfo"), geneSymbol, position, position)
    cursor.execute(query)
    iii = cursor.fetchall()
    if len(iii) == 0:
        return 0
    result = []
    #print list(ii)[0]
    for ii in iii:
        if ii != None:
            ii = list(ii)
        else:
            continue
        if len(ii) == 5:
            query = """SELECT `seq` FROM %s WHERE `rec_name` = '%s';""" % (getConfig("datasets", "intronSeq"), ii[0])
            cursor.execute(query)
            seq = cursor.fetchone()
            if seq != None:
                seq = seq[0]
            else:
                continue
            ii.append(seq)
        else:
            continue
        if ii[2] > ii[3]:   #negative strand
            tmp = ii[2]
            ii[2] = ii[3]
            ii[3] = tmp
        result.append(ii)
    '''    
    return result

def isChangeSS(raw, edited, sig, type=3):
    rawDict = []
    editDict = []
    filteredRawDict = {}
    filteredEditDict = {}
    filterThreshold = string.atof(getConfig("splicing", "mesSCT"))

    if type == 5:
        '''
        for ss5 in splicingSite5Motif.finditer(raw):
            rawDict.append((ss5.start(), ss5.end(), ss5.group(0)))
        for ss5 in splicingSite5Motif.finditer(edited):
            editDict.append((ss5.start(), ss5.end(), ss5.group(0)))
        '''
        rawDict.append((0, 9, raw))
        editDict.append((0, 9, edited))
        rawScore = callMaxEntScan5ss([ss5[2]+"\n" for ss5 in rawDict], sig)
        editScore = callMaxEntScan5ss([ss5[2]+"\n" for ss5 in editDict], sig)
    else:
        '''
        for ss3 in splicingSite3Motif.finditer(raw):
            rawDict.append((ss3.start(), ss3.end(), ss3.group(0)))
        for ss3 in splicingSite3Motif.finditer(edited):
            editDict.append((ss3.start(), ss3.end(), ss3.group(0)))
        '''
        rawDict.append((0, 9, raw))
        editDict.append((0, 9, edited))
        rawScore = callMaxEntScan3ss([ss3[2]+"\n" for ss3 in rawDict], sig)
        editScore = callMaxEntScan3ss([ss3[2]+"\n" for ss3 in editDict], sig)
    if (rawScore==editScore): 
        return 0, 0, 0, 0
    
    rawRecordsNum = len(rawDict)
    editRecordsNum = len(editDict)
    for i in range(rawRecordsNum):
        score = string.atof(rawScore[i])
        #if score >= filterThreshold:
        filteredRawDict[str(rawDict[i][0])+","+str(rawDict[i][1])] = (score, rawDict[i][2])
    for i in range(editRecordsNum):
        score = string.atof(editScore[i]) 
        if score >= filterThreshold:
            filteredEditDict[str(editDict[i][0])+","+str(editDict[i][1])] = (score, editDict[i][2])
    #del(rawDict);del(editDict);del(rawScore);del(editScore)
    inactiveSplicingSite = []
    newPotentialSplicingSite = []
    for (k,v) in filteredRawDict.items():
        if k not in filteredEditDict.keys():
            return (1*type, v[0], 0, 0)  #inactive
        elif (v[0] != filteredEditDict[k][0]):
            variation = (filteredEditDict[k][0] - v[0]) / v[0]
            if (variation > float(getConfig("splicing", "mesVT"))):
                return (2*type, v[0], filteredEditDict[k][0], variation)  #enhanced splicing
            elif (variation < -1*float(getConfig("splicing", "mesVT"))):
                if (variation < -1*float(getConfig("splicing", "mesVT"))):
                    return (1*type, v[0], filteredEditDict[k][0], variation)  #weakened splicing
                else:
                    return (3*type, v[0], filteredEditDict[k][0], variation)  #weakened splicing

    for (k, v) in filteredEditDict.items():
        if k not in filteredRawDict.keys():
            return (4*type, 0, v[0], 0)  #new splicing site
    return 0, 0, 0, 0

def defineAGEZ(seq, pos3SS):
    #window = seq[]
    pos = pos3SS - 12
    firstTag = 1
    while pos > 0:
        if (seq[pos] == 'G' and seq[pos-1] == 'A'):
            if firstTag == 1:
                firstTag = 0
            else:
                return pos
        pos -= 1
    return 0

def bpPositionWeightMatrix(seq, start, end):
    score = 0.0
    maxScore = 0.0
    for i in range(start, end - 7):
        for j in range(7):
            score += branchSiteMatris[j][seq[i+j].upper()]
        if score >= bpThreshold:
            return score
        score = 0.0
    return 0

def isChangeBranchSite(raw, edited, start, end):
    rawScore = bpPositionWeightMatrix(raw, start, end)
    editedScore = bpPositionWeightMatrix(edited, start, end)

    if rawScore == 0 and editedScore != 0:
        return(16, rawScore, editedScore, editedScore)
    elif rawScore != 0 and editedScore == 0:
        return(17, rawScore, editedScore, editedScore)
    elif rawScore == 0 and editedScore == 0:
        return 0, 0, 0, 0

    variation = (editedScore - rawScore) / rawScore
    if variation < 0 and ( -1*variation > bpVarThreshold):
        return (18, rawScore, editedScore, variation)
    elif variation > 0 and ( variation > bpVarThreshold):
        return (19, rawScore, editedScore, variation)
    return 0, 0, 0, 0

def isMutChangeSplicing(geneName, chr, mutPos, raw='A', muted='G', jid=0):
    flag = 0
    ssRangeFlag = 0
    intronSet = getIntronInfo(geneName, chr, mutPos)
    uniSig = generateSSSig(geneName, mutPos, jid, raw+'>'+muted)
    if intronSet == 0: 
        return 0, 0, 0, 0, 0
    for intron in intronSet:
        orderOfIntron = intron[0].split('_')[2]
        transcript = intron[0].split('_')[0].split('.')[0]
        rangeOf5SS = range(10, 16)
        rangeOf3SS = range(intron[3] - intron[2] - 10, intron[3] - intron[2]+10)
        ter5Of3SS = intron[3] - intron[2] + 6
    
        if intron[4] == '-':
            relativeMutPos = intron[3] - mutPos + 10
        else:
            relativeMutPos = mutPos - intron[2] + 9
        mutSeq = intron[5][:relativeMutPos]+muted+intron[5][relativeMutPos+1:]
        #print mutPos, intron[5], mutSeq

        if relativeMutPos in rangeOf5SS:
            ssRangeFlag = 1
            if (intron[5][relativeMutPos].lower() != raw.lower()):
                print "Please check your sequence and mutation position(%s, %s, %s)" % (mutPos, intron[5][relativeMutPos].lower(), raw.lower())
                continue
                #return 0, 0, 0, 0, 0
            ss5ret = list(isChangeSS(intron[5][7:16], mutSeq[7:16], uniSig, 5))
            if ss5ret != 0:
                ss5ret.append(orderOfIntron)
                #return ss5ret
                if ss5ret[0] != 0:
                    flag = 1
                    recordSplicing(geneName, chr, mutPos, ss5ret[0], ss5ret[1], ss5ret[2], ss5ret[3], ss5ret[4], jid, transcript)
        elif relativeMutPos in rangeOf3SS:
            ssRangeFlag = 1
            if (intron[5][relativeMutPos].lower() != raw.lower()):
                print "Please check your sequence and mutation position(%s, %s, %s)" % (mutPos, intron[5][relativeMutPos].lower(), raw.lower())
                #return 0, 0, 0, 0, 0
                continue
            ss3ret = list(isChangeSS(intron[5][-30:-7], mutSeq[-30:-7], uniSig, 3))
            if ss3ret != 0:
                ss3ret.append(orderOfIntron)
                #return ss3ret
                if ss3ret[0] != 0:
                    flag = 1
                    recordSplicing(geneName, chr, mutPos, ss3ret[0], ss3ret[1], ss3ret[2], ss3ret[3], ss3ret[4], jid, transcript)
        else:
            bs = defineAGEZ(intron[5], ter5Of3SS)
            if relativeMutPos in range(bs, ter5Of3SS):
                ssRangeFlag = 1
                bsRet = list(isChangeBranchSite(intron[5], mutSeq, bs, ter5Of3SS))
                bsRet.append(orderOfIntron)
                #return bsRet
                if bsRet[0] != 0:
                    flag = 1
                    recordSplicing(geneName, chr, mutPos, bsRet[0], bsRet[1], bsRet[2], bsRet[3], bsRet[4], jid, transcript)
        #return 0, 0, 0, 0, 0
    if ssRangeFlag == 0:
        print geneName, mutPos, 'not in ss region'
    return flag

def factory(gene, chr, pos, jobId, extend=0, ref=None, edited=None):
    counter = 0
    for x in gene.index:
        if extend == 1:
            t = isMutChangeSplicing(gene.ix[x], chr.ix[x], pos.ix[x], ref.ix[x], edited.ix[x], jid=jobId)
        else:
            t = isMutChangeSplicing(gene.ix[x], chr.ix[x], pos.ix[x], jid=jobId)
        if t != 0:
            counter += 1
    return counter

def recordSplicing(gene, chr, pos, t, r, n, v, order, jobId, transcript):
    try:
        query = """INSERT INTO `splicing_event` (`gene`, `chromosome`, `pos`, `type`, `raw_score`, `new_score`, `variation`, `order`, `job`, `transcript`) VALUES ('%s', '%s', %s, %s, %s, '%s', %s, %s, %s, '%s');""" % (gene, chr, pos, t, r, n, v, order, jobId, transcript)
        cursor.execute(query)
        con.commit()
    except Exception, e:
        print "Splicing Error:"+str(e)
        return 0
