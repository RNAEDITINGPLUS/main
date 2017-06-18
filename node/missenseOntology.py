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

import pysam, pickle, repConfig
from Bio.Seq import transcribe, translate
from Bio.Alphabet import IUPAC
from database import mysqlConnection
aaTable = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid', 'C': 'Cysteine', 'Q': 'Glutamine', 
    'E': 'Glutamic acid', 'G': 'Glycine', 'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'M': 'Methionine', 
    'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine', 'T': 'Threonine','W': 'Tryptophan', 'Y': 'Tyrosine', 
    'V': 'Valine', 'K': 'Lysine', '*': 'Stop Codon'
}
sdf = open(repConfig.getConfig("datasets", "cdslibrary"))
seqDict = pickle.load(sdf)
ceTable = pysam.Tabixfile(repConfig.getConfig("datasets", "cdsbed"))
cnx, cursor = mysqlConnection()
def isMissenseMut(chr, pos, job):
    global seqDict, ceTable
    for hit in ceTable.fetch(reference=chr, start=int(pos), end=int(pos)+1, parser=pysam.asBed()):
        chromosome, start, end, name, gene, strand, relStart = hit
        tmp = name.split('_')
        transcriptID = tmp[0]
        rawSeq = seqDict[transcriptID]
        if strand == '+':
            relPos = int(pos) - int(start) + int(relStart)-1
        else:
            relPos = int(end) - int(pos) + int(relStart)
        editedSeq = rawSeq[:relPos]+'G'+rawSeq[relPos+1:]
        rawSeqr = transcribe(rawSeq)
        rawPr = translate(rawSeqr)
        editedSeqr = transcribe(editedSeq)
        editedPr = translate(editedSeqr)
        tag, rp, fa, ta = seqCompare(rawPr, editedPr)
        if tag:
            recordMissense(chr, pos, job, gene, transcriptID, rp, aaTable[fa], aaTable[ta])
            return 1
        else:
            return 0

def seqCompare(raw, edited):
    tag = 0
    pos = 0
    for index, aa in enumerate(raw):
        if aa != edited[index]:
            tag = 1
            pos = index - 1
            return tag, pos, aa, edited[index]
    return 0, 0, 0, 0

def recordMissense(chr, pos, job, gene, transcript, relp, fromA, toA):
    try:
        sql = """INSERT INTO `%s` (`job`, `chromosome`, `position`, `gene`, `transcript`, `relpos`, `fr`, `to`) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');""" % (repConfig.getConfig("datasets", "mist"), job, chr, pos, gene, transcript, relp, fromA, toA)
        cursor.execute(sql)
        cnx.commit()
    except Exception, e:
        print e
        return 0
    
def factory(chr, pos, jobId):
    counter = 0
    for x in chr.index:
        if isMissenseMut(chr.ix[x], pos.ix[x], jobId):
            counter += 1
    return counter