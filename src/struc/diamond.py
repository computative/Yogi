from struc.structure import Structure


class Diamond(Structure):

    def __init__(self,spp):
        sp1,sp2 = spp
        a = self.a[sp1][sp2]
        self.atoms = {sp1:[],sp2:[]}
        self.atoms[sp1].extend([
                [0.00,0.00,0.00],
                [0.00,a/2,a/2],
                [a/2,0.00,a/2],
                [a/2,a/2,0.00]
        ])
        self.atoms[sp2].extend([
                [a/4,a/4,a/4],
                [a/4,3*a/4,3*a/4],
                [3*a/4,a/4,3*a/4],
                [3*a/4,3*a/4,a/4]
        ])
        self.rep = [
            [a,0,0],
            [0,a,0],
            [0,0,a]
        ]
