from .dna import DNA


def extract(seq):
    """
    extract sequence information, such as amino acid sequence, from a string.
    """
    
    import re
    return "".join(re.findall("[a-zA-Z]+", seq))

def calc_bind(info, Kd):
    """
    Simple kinetic model: Kd = [A][B]/[AB].
    info = [conc1, conc2] or {name1:conc1, name2:conc2}
    """
    
    import numpy as np

    if (type(info) not in [dict, list, tuple]):
        raise TypeError("wrong info type.")
    if (len(info) != 2):
        raise ValueError("wrong info length.")
    n1 = "mole-1"
    n2 = "mole-2"
    if (type(info) == dict):
        items = iter(info.items())
        n1, c1 = next(items)
        n2, c2 = next(items)
    elif (type(info) in [list, tuple]):
        c1, c2 = info
    else:
        raise TypeError

    cb = (c1 + c2 + Kd - np.sqrt((c1 + c2 + Kd)** 2 - 4 * c1 * c2)) / 2

    return {f"[{n1}]_total": c1,
            f"[{n2}]_total": c2,
            f"[{n1}:{n2}]": cb,
            f"[{n1}]_free": c1 - cb,
            f"[{n2}]_free": c2 - cb
            }

def effconc(n_aa, c=9):
    """
    Timpe et al., 1995.
    Infinite plane model is supposed.
    -----------------------------------------------------
    For charged polypeptide, c=9
        (experimentally determined using poly-K and poly-E)
    For poly-G, c=2
    For poly-P, c=116
    """
    return 20*(n_aa*c)**-1.5

def extract_alignment(str_, n_seq=2, match="last"):
    """
    Extract alignment result, such as clustalw.
    """
    if (n_seq < 2):
        raise ValueError("n_seq must be >=2.")
    if (match == "last"):
        pass
    elif (match == "center"):
        if (n_seq == 2):
            pass
        else:
            raise ValueError("option match='center' can only be selected when n_seq=2.")
    else:
        raise ValueError("match should be either 'last' or 'center'")
    
    seqdict={} # e.g. {"protein1": "GGST...", "protein2": "VLKG..."}
    strs = str_.split("\n")
    group = strs[0:n_seq]
    for i in range(n_seq):
        protein_name, _ = group[i].split()
        seqdict[protein_name] = ""

    pos = 0
    while len(strs) - pos > 0:
        # a group will look like:
        # protein1  --GGST...
        # protein2  VLKG--...
        #             .*
        group = strs[pos:pos + n_seq + 2]
        if (match == "last"):
            for i in range(n_seq):
                protein_name, subsequence = group[i].split()
                seqdict[protein_name] += subsequence
        else:
            protein_name, subsequence = group[0].split()
            seqdict[protein_name] += subsequence
            protein_name, subsequence = group[2].split()
            seqdict[protein_name] += subsequence

        pos += n_seq + 2
    
    return seqdict