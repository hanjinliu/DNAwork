comp = {"A": "T", "T": "A", "G": "C", "C": "G",
        "a": "t", "t": "a", "g": "c", "c": "g"
        }

comp_ex = {"A": "T", "T": "A", "G": "C", "C": "G",
           "a": "t", "t": "a", "g": "c", "c": "g",
           "N": "N", "n": "n", "W": "W", "w": "w",
           "S": "S", "s": "s", "R": "Y", "r": "y",
           "Y": "R", "y": "r"
           }

codon = {
        'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
        'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
        'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : '*',
        'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

        'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R',
        'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R',
        'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R',
        'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R',

        'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S',
        'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S',
        'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R',
        'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R',

        'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G',
        'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G',
        'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G',
        'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G',
        }

def rc(seq: str):
    """
    Convert the sequence to reverse complementary one.
    """
    return "".join([comp[s] for s in seq[::-1]])

def translate(seq: str, autostop=True):
    """
    Translate DNA sequence to protein sequence.
    """
    p = ""
    seq = seq.upper()
    for i in range(0, int(len(seq) / 3)):
        c = seq[3 * i : 3 * i + 3]
        if (c in codon):
            p += codon[c]
            if (codon[c] == "*" and autostop):
                break
        else:
            p += 'X'
    return p

def iter_findall(template: str, seq: str):
    """
    From template, iteratively find seq.
    """
    pos = 0
    while True:
        pos_ = template[pos:].find(seq)
        if (pos_ < 0):
            break
        pos += pos_
        yield pos
        pos += 1