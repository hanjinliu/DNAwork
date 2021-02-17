import re

class DNAfeature:
    def __init__(self, info=""):
        self.start = 0
        self.end = 0
        self.tag = "tag"
        self.color = None
        reg = re.compile(r"(\s|\"|\n)*?")
        lines = info.split("/")
        for line in lines:
            if(".." in line):
                start, end = line.split("..")
                self.start = int(reg.sub("", start))
                self.end = int(reg.sub("", end))
            elif(line.startswith("label=")):
                self.tag = reg.sub("", line.lstrip("label="))
            elif(line.startswith("ApEinfo_fwdcolor=")):
                self.color = reg.sub("", line.lstrip("ApEinfo_fwdcolor="))
        
    
    def __len__(self):
        return self.end - self.start

    def copy(self):
        new = self.__class__()
        new.tag = self.tag
        new.color = self.color
        new.start = self.start
        new.end = self.end
        return new


    def short_copy(self, start, end, total_len):
        if (start < end):
            new = self._cover_inside(start, end)
        else:
            new = self._cover_outside(start, end, total_len)
        
        return new

    def _cover_inside(self, start, end):
        new = self.copy()

        if (self.start < start and start <= self.end and self.end < end):
            #  ====        <- feature
            # --[------]--
            #
            new.start = 0
            new.end = self.end - start
        elif (start <= self.start and self.start < end and end < self.end):
            #        ====  <- feature
            # --[------]--
            new.start = self.start - start
            new.end = self.end - start
        elif (start <= self.start and self.end <= end):
            #     ====     <- feature
            # --[------]--
            new.start = self.start - start
            new.end = self.end - start
        elif (self.start < start and end < self.end):
            #  ==========  <- feature
            # --[------]--
            new.start = 0
            new.end = self.end - self.start
        else:
            new = None
        return new
    
    def _cover_outside(self, start, end, total_len):
        new = self.copy()
        if (self.start < end and end <= self.end and self.end < start):
            #   ====       <- feature
            # ----]--[----
            new.start = self.start + total_len - start
            new.end = end + total_len - start
        elif (end <= self.start and self.start < start and start < self.end):
            #       ====   <- feature
            # ----]--[----
            new.start = 0
            new.end = self.end - start
        elif (self.start < start and end <= self.end):
            #    ======    <- feature
            # ----]--[----
            new1 = self.__class__()
            new1.tag = new.tag
            new1.color = new.color
            new2 = self.__class__()
            new2.tag = new.tag
            new2.color = new.color
            
            new1.start = 0
            new1.end = self.end - start
            new2.start = self.start + total_len - start
            new2.end = self.end + total_len - start
            new = [new1, new2]
        elif (self.end < end):
            #  ==          <- feature
            # ----]--[----
            new.start = total_len - start + self.start
            new.end = total_len - start + self.end
        elif (start < self.start):
            #          ==  <- feature
            # ----]--[----
            new.start = self.start - start
            new.end = self.end - start
        else:
            new = None
        return new