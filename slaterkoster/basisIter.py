import numpy as np

class BasisIterator:
    def __init__(self, R, functions):
        self.contents = [ (R,f)  for R in R  for f in functions ]
        self.pos = 0
        self.len = len(self.contents)

    def __next__(self):
        if self.pos >= self.len :
            raise StopIteration
        elt = self.contents[self.pos]
        self.pos += 1
        return elt


class Basis:

    def __init__(self, crystal, functions ):
        self.crystal = crystal
        self.functions = functions
        self.R = []
        for elts in crystal.nuclei.values():
            for R in elts.values():
                self.R.append(R)

    def __iter__(self):
        return BasisIterator(self.R, self.functions)

    def __len__(self):
        return len(self.crystal)*len(self.functions)