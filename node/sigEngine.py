import hashlib
from repConfig import getConfig

def generatemiRNASig(mirna, position, chromosome, jid, event="A>I"):
    sig = hashlib.md5()
    sig.update(mirna+str(chromosome)+str(position)+str(jid)+event+getConfig("safety", "sault"))
    return sig.hexdigest()

def generateUTR3Sig(utr, position, chromosome, jid, event="A>I"):
    sig = hashlib.md5()
    sig.update(utr+str(chromosome)+str(position)+str(jid)+event+getConfig("safety", "sault"))
    return sig.hexdigest()

def anySig(combined):
    sig = hashlib.md5()
    sig.update(combined)
    return sig.hexdigest()

def generateSSSig(gene, position, jid, event="A>I"):
    sig = hashlib.md5()
    sig.update(gene+str(position)+str(jid)+event+getConfig("safety", "sault"))
    return sig.hexdigest()

def randSig():
    import time
    sig = hashlib.md5()
    sig.update(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    return sig.hexdigest()