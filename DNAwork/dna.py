import os
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from .func import rc, translate, iter_findall
from .match import DNAmatch
from .feature import DNAfeature


def asDNA(obj, name=""):
    """
    Convert str to DNA object, if needed.
    """
    if (isinstance(obj, str)):
        obj = DNA(obj, name)
    elif (not isinstance(obj, DNA)):
        raise TypeError("Input must be either `string` or `DNA`")
    return obj

def read_seq_from_file(path: str):
    """
    Open file and read sequence and/or other information.
    """
    with open(path, "r") as f:
        raw = f.read()
    
    if (os.path.splitext(path)[1] == ".ape"):
        # get name and shape
        raw_line0 = raw.split("\n")[0]
        name = re.findall("LOCUS.*[0-9]+ bp", raw_line0)[0].split()[1]
        shape = "circular" if "circular" in raw_line0.split("bp")[-1] else "linear"
        # split info
        features0, seq0 = raw.split("ORIGIN")
        # get features
        if("misc_feature" in features0):
            feature_info_list = features0.split("misc_feature")[1:]
            feature_list = []
            for info in feature_info_list:
                try:
                    feature = DNAfeature(info)
                    if (len(feature) > 15):
                        feature_list.append(feature)
                except:
                    print("Exception happened during reading ape file.")
        else:
            feature_list = []

        # get sequence
        seq = "".join(re.sub("[0-9 \n/]+", "", seq0[6:]))
    else:
        name = ""
        shape = "linear"
        seq = raw
        feature_list = []
    
    return seq, name, feature_list, shape

class DNA:
    count = 1
    primer_color = "blue"
    digest_color = "red"
    def __init__(self, seq, name=""):
        feature_list = []
        shape = "linear"

        # import from txt or ape file
        if (os.path.isfile(seq)):
            seq, name0, feature_list, shape = read_seq_from_file(seq)
            if (name == ""):
                name = name0
        
        # import from FASTA format
        elif (seq.startswith(">")):
            name0, seq, *_ = seq.split("\n")
            name0 = name0[1:]
            if (name == ""):
                name = name0

        if (re.fullmatch(r"[ATGCatgc]+", seq) is None):
            print(seq)
            print(name)
            raise ValueError("Input has character(s) other than ATGC.")

        self.f = seq.upper()
        self.r = rc(seq).upper()
        self._shape = shape
        self.feature_list = feature_list
        self._set_name(name)
    
    def __getitem__(self, key):
        return self.f[key]
    
    @property
    def shape(self):
        return self._shape
    
    @shape.setter
    def shape(self, value):
        if (value in ["linear", "circular"]):
            self._shape = value
        else:
            raise ValueError("shape must be set to either `linear` or `circular`.")

    def __len__(self):
        return len(self.f)
    
    def __str__(self):
        return self.f
    
    def __repr__(self):
        return f"{self.name} (length = {len(self)}, {self.shape})"

    def __contains__(self, seq):
        f, r = self._findall(seq)
        return len(f) + len(r) > 0
    
    def PCR(self, forward, reverse, min_match=15):
        """
        Conduct PCR using `self` as template.
        min_match: The minimum length that defines two sequences are complementary.
        """
        # type check
        forward = asDNA(forward, name="fw primer")
        reverse = asDNA(reverse, name="rv primer")

        # get list of DNAmatch object
        match = self.find_match(forward, min_match) + self.find_match(reverse, min_match)
        ans = None

        # plot
        if(match):
            self._plot(match)

        # print result
        if (len(match) == 0):
            print(" No PCR product obtained. No match found.")
        elif (len(match) == 1):
            print(" No PCR product obtained. Only one match found.")
        elif (len(match) > 2):
            print(f" Too many matches ({len(match)} matches found).")
        elif (match[0].dir == match[1].dir):
            print(" Each primer binds to the template in the same direction.")
        else:
            ans = self._do_PCR(match)
            ans.name = f"{self.name}-amp"
            plt.text(match[0].start * 0.4 + match[1].end * 0.6, 0.22, f"{len(ans)} bp",
                     ha="center", size = 12, color = "orangered")
            plt.show()

        return ans
    
    def __add__(self, other):
        return self.ligate(other)
    
    def translate(self, startswith=0, endswith="until stop codon"):
        seq = self.f
        autostop = True
        if (isinstance(startswith, int)):
            aa0 = startswith
        elif (isinstance(startswith, str)):
            aa0 = seq.find(startswith.upper())
        else:
            raise TypeError("`startswith` must be int or str.")
        
        if (aa0 < 0 or aa0 >= len(seq)):
            ValueError("`startswith` is not appropriate.")

        if (endswith == "until stop codon"):
            aa1 = -1
        elif (isinstance(endswith, int)):
            aa1 = endswith
        elif (isinstance(endswith, str)):
            aa1 = seq[aa0:].find(endswith.upper()) + len(endswith) + aa0
        else:
            raise TypeError("`startswith` must be int or str.")
        
        if (aa1 > len(seq)):
            raise ValueError("Keyword arguments are not appropriate.")
            
        seq = seq[aa0:aa1]

        return translate(seq, autostop=autostop)

    def ligate(self, other):
        other = asDNA(other)
        l = len(self)
        self.f = self.f + other.f
        self.r = other.r + self.r
        for f in other.feature_list:
            f = f.copy()
            f.start += l
            f.end += l
            self.feature_list.append(f)
        return self
    
    def iter_feature(self, feature_name):
        for f in self.feature_list:
            if (f.tag == feature_name):
                yield f

    def crop(self, start, end):
        """
        crop out
        """
        if (start >= end and self.shape == "linear"):
            raise ValueError("start < end for linear DNA.")
        elif (start < end and self.shape == "linear"):
            _, seq = self.split(start)
            seq, _ = seq.split(end - start)
        elif (start < end and self.shape == "circular"):
            opened = self.open(start)
            seq, _ = opened.split(end - start)
        elif (start >= end and self.shape == "circular"):
            opened = self.open(start)
            seq, _ = opened.split(len(opened) - start + end)
        else:
            raise ValueError("Unknown error...")
        
        seq.name = self.name + "-frag"
        return seq


    def split(self, pos):
        """
        Split linear DNA into two fragments.
        Do NOT update self.
        """
        if (self.shape == "circular"):
            raise TypeError("Cannot split circular DNA.")

        frag_l = self.__class__(self.f[:pos], name="left")
        frag_r = self.__class__(self.f[pos:], name="right")
        for f in self.feature_list:
            f = f.copy()
            if (f.end <= pos):
                frag_l.feature_list.append(f)
            elif (pos <= f.start):
                f.start -= pos
                f.end -= pos
                frag_r.feature_list.append(f)
            else:
                f2 = f.copy()
                f.end = pos
                f2.start = 0
                f2.end -= pos
                frag_l.feature_list.append(f)
                frag_r.feature_list.append(f2)

        return frag_l, frag_r

    def open(self, pos):
        """
        Open circular DNA.
        Do NOT update self.
        """
        if (self.shape == "linear"):
            raise TypeError("Cannot open linear DNA.")
        self.shape = "linear"
        frag_l, frag_r = self.split(pos)
        self.shape = "circular"
        return frag_r.ligate(frag_l)

    def InFusion(self, insert):
        insert = asDNA(insert, name="insert")
        if (self.shape == "circular" or insert.shape == "circular"):
            raise TypeError("Both vector and insert must be linear DNA.")
        if (len(self) < 30):
            raise ValueError(f"`{self.name}` is too short.")
        if (len(insert) < 30):
            raise ValueError("insert is too short.")
        
        pos0 = len(self) // 2
        frag_l, frag_r = self.split(pos0)

        if (frag_l[:15] != insert[-15:]):
            raise ValueError("Mismatch! Check carefully:\n"\
                            f"--{insert[-20:]}\n"\
                            f"       {frag_l[:20]}--")
        if (frag_r[-15:] != insert[:15]):
            raise ValueError("Mismatch! Check carefully:\n"\
                            f"       {insert[:20]}--\n"\
                            f"--{frag_r[-20:]}")
        
        frag_r, _ = frag_r.split(len(frag_r) - 15)
        _, frag_l = frag_l.split(15)
        out = frag_r + insert + frag_l
        out.name = f"{insert.name} in {self.name}"
        out.shape = "circular"

        return out

    def digest(self, seq):
        """
        seq = G^GATCC
        """
        try:
            l, r = seq.split("^")
        except:
            raise ValueError("Input must be in a `G^GATCC`-like form.")
        if (rc(l + r) != l + r):
            raise ValueError("Input must be palindrome.")
        
        c = self.__class__.digest_color
        if (len(l) < len(r)):
            p = [-len(self) / 200, len(self) / 200]
        elif (len(l) > len(r)):
            p = [len(self) / 200, -len(self) / 200]
        else:
            p = [0, 0]

        self._plot([])
        match = self.find_match(l + r)
        x_list = []
        for m in match:
            if (m.dir == "->"):
                x_list.append((m.start + m.end) / 2)
        
        lens = []
        poss = []
        if(len(x_list) > 0):
            x_list.sort()
            edge = 0
            for i, xc in enumerate(x_list):
                x = [xc + p[0], xc + p[0], xc + p[1], xc + p[1]]
                y = [-0.03, -0.08, -0.08, -0.13]
                plt.plot(x, y, color=c, lw=2, zorder=1001)
                lens.append(xc - edge + (1 + (i > 0)) * abs(len(l) - len(r)))
                poss.append((xc + edge) / 2)
                edge = xc
            
            lens.append(len(self) - x_list[-1] + abs(len(l) - len(r)))
            poss.append((len(self) + x_list[-1]) / 2)
            
            if (self.shape == "circular"):
                if (lens[0] < lens[-1]):
                    poss[0] = poss[-1]
                else:
                    pass
                poss.pop()
                lens[0] += lens.pop()
                
            for i in range(len(lens)):
                plt.text(poss[i], 0, f"{int(lens[i])} bp", ha="center",
                         size = 12, color = "orangered")

        plt.show()
        return lens

    def find_match(self, seq, min_match=15):
        if (min_match <= 0):
            raise ValueError("`min_match` must be positive value")

        seq = asDNA(seq, name="temp")

        match = []
        fw_pos, rv_pos = self._findall(seq[-min_match:])

        if(min_match > len(seq)):
            min_match = len(seq)
            
        # 1: forward check
        for pos in fw_pos:
            prpos = len(seq) - min_match # position on seq
            while(pos > 0 and prpos > 0 and self.f[pos-1] == seq[prpos-1]):
                pos -= 1
                prpos -= 1
            
            match.append(DNAmatch("->", pos, len(seq) - prpos + pos, seq))
        
        # 2: reverse check
        for pos in rv_pos:
            prpos = min_match  # position on seq
            pos3 = pos + min_match
            while (pos3 < len(self) and prpos < len(seq) and
                    pos3 > min_match - 1 and prpos > 0 and
                    self.f[pos3] == seq.r[prpos]):
                pos3 += 1
                prpos += 1
            
            match.append(DNAmatch("<-", pos, pos3, seq))

        return match
    
    def plot(self):
        self._plot([])
        plt.show()
        return None
    
    def plot_match(self, seqs, min_match=15):
        if (not isinstance(seqs, (list, tuple))):
            seqs = [seqs]
         
        match = []
        for seq in seqs:
            match += self.find_match(seq, min_match)
        
        self._plot(match)
        plt.show()
        return None

    def _do_PCR(self, match):
        # order: match[0] < match[1]
        if (match[0].start > match[1].start):
            match[0], match[1] = match[1], match[0]
        
        features = []

        if (match[0].dir == "->" and match[1].dir == "<-"):
            mf, mr = match
            product_seq = self[mf.start:mr.end]

        elif (self.shape == "linear"):
            print(" No PCR product obtained.")
            return None
        else:
            mr, mf = match
            product_seq = (self.f * 2)[mf.start:mr.end + len(self)]
        
        for f in self.feature_list:
            new_feature = f.short_copy(mf.start, mr.end, len(self))
            if(isinstance(new_feature, DNAfeature)):
                features.append(new_feature)
            elif (isinstance(new_feature, list)):
                features += new_feature

        # deal with flanking regions
        out = len(mf.ref) - len(mf)
        if (out > 0):
            product_seq = mf.ref[:out] + product_seq
        out = len(mr.ref) - len(mr)
        if (out > 0):
            product_seq = product_seq + mr.ref.r[-out:]

        ans = DNA(product_seq, name="temp")
        ans.feature_list = features

        return ans

    def _plot(self, match):
        arrow_color = self.__class__.primer_color
        plt.figure(figsize=(10, max(0.3*len(match) + 1.2, 1.6)))
        ax = plt.subplot(111)

        for i, matchobj in enumerate(match):
            dr, pos0, pos1, seq = matchobj.items()  # match at pos0:pos1
            # x0 -> x1 or x1 <- x0
            if(dr == "->"):
                x0 = pos0
                x1 = max(pos1, min(pos0 + len(self) / 40, len(self)))
                mk = ">"
                x_ov = x0 - len(self) / 64
                ha_ = "right"
            else:
                x0 = pos1
                x1 = min(pos0, max(pos1 - len(self)/40, 0))
                mk = "<"
                x_ov = x0 + len(self) / 64
                ha_ = "left"

            # plot line
            plt.plot([x0, x1], [i * 0.1, i * 0.1], color=arrow_color, lw=3)
            # plot arrow head
            plt.scatter(x1, i * 0.1, c=arrow_color, marker=mk, s=30)
            # primer regions that do not bind to the template
            if(pos1 - pos0 < len(seq)):
                plt.plot([x0, x_ov], [i * 0.1, i * 0.1 + 0.05], color=arrow_color, lw=3)
                plt.text(x_ov, i * 0.1 + 0.01, seq.name, size=12, color="black", ha = ha_)    
            else:
                plt.text(pos0, i * 0.1 + 0.01, seq.name, size=12, color="black", ha = ha_)
            
        plt.plot([-0.05, len(self)], [-0.08, -0.08], "black", linewidth = 2)        
        for f in self.feature_list:
            try:
                r = patches.Rectangle(xy=(f.start, -0.11), width=len(f), height=0.06,
                                      fc=f.color, ec="black", label=f.tag, zorder=1000)
                ax.add_patch(r)
            except:
                print(f.tag, f.color, f.start, f.end)
        
        plt.text(len(self) * 0.5, -0.2, self.name, ha="center", size=13, color="black")
        plt.text(len(self), -0.2, self.shape, ha="right", size=11, color="brown")
        
        plt.grid(axis="x")
        plt.ylim(-0.24, len(match)*0.1 + 0.14)
        plt.legend(bbox_to_anchor = (1.05, 1), loc = "upper left", borderaxespad = 0)
        plt.tick_params(labelleft = False, bottom = False, left = False,
                        right = False, top = False)
        plt.tight_layout()
        
        return None

    def _set_name(self, name):
        if (not isinstance(name, str)):
            raise TypeError("type of `name` must be `string`")
        if (name == ""):
            self.name = f"seq-{self.__class__.count}"
        else:
            self.name = name
        self.__class__.count += 1
        return None    
    
    def _findall(self, seq):
        """
        Find all the match position, including reverse complement.
        e.g.
        Suppose self.f = TTTAGGAGGTTCCTT, seq = AGG, then
        fw_pos = [3, 6], rv_pos = [13]
        """
        # type check
        if (isinstance(seq, self.__class__)):
            seq = seq.f
        elif (not isinstance(seq, str)):
            raise TypeError("type of `seq` must be either `string` or `DNA`")
        seq = seq.upper()

        # all the matched position 
        fw_pos = []
        for pos in iter_findall(self.f, seq):
            fw_pos.append(pos)
        
        rv_pos = []
        for pos in iter_findall(self.f, rc(seq)):
            rv_pos.append(pos)

        return fw_pos, rv_pos

        