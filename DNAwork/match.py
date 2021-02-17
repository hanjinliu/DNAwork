class DNAmatch:
    def __init__(self, dir_, start, end, ref):
        self.dir = dir_     # direction. `->` or `<-`
        self.start = start
        self.end = end
        self.ref = ref      # reference sequence.
    
    def __len__(self):
        return self.end - self.start
    
    def items(self):
        return self.dir, self.start, self.end, self.ref