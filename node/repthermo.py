from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

def miScore(seq):
    """miScore
    Kume, H., Hino, K., Galipon, J., & Ui-Tei, K. (2014).
    A-to-I editing in the miRNA seed region regulates target mRNA selection and silencing efficiency.
    Nucleic Acids Research, 42(15), 10050�C10060. http://doi.org/10.1093/nar/gku662
    """
    try:
        score = mt.Tm_NN(seq[1:9]) - 0.5 * mt.Tm_NN(seq[0:5])
    except:
        score = -1995

    return score