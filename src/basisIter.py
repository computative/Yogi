class BasisIterator:
    def __init__(self, basis):
        self.cont = [(R,f)  for R in basis.crystal.nuclei \
                            for f in basis.functions ]
        self.pos = 0
        self.len = len(self.cont)

    def __next__(self):
        if pos >= self.len :
            raise StopIteration
        elt = self.cont[pos]
        pos += 1
        return elt


class Basis:

    def __init__(self, crystal, functions ):
        self.crystal = crystal
        self.functions = functions

    def __iter__(self):
        return BasisIterator(self)

    def __len__(self):
        return len(self.crystal)*len(functions)