class electron:
    def __init__(self):
        # hvilken atom tilhorer den
        # hvilke bond har den?
        return True
"""
class nucleus:
    def __init__(self):
        # hvilken atom tilhorer den
        # hvilken ladning har den
        return True
"""
class atom:
    def __init__(self,Q,):
        # hvilken atom tilhorer den
        # hvilken ladning har den
        return True

"""
class lattice:
    def __init__(self):
        # coords
        
        # hvilke elektroner har den?

        # hvilken id har den?
"""

class Solver:

    def __init__(self, N):
        # sette inn nuclei og elektroner
        self.data = [1,3,2,4]
        
    def __iter__(self):
        return iter(self.data)
    
    def __next__(self,itr):
        return next(itr)
    
    def perturb(self,Id,shift):
        # endre posisjon til elektroner og nucleus med Id.
        # endre nermeste naboer til den bevegede partikkelen
        #     og naboene selv
        return True

    def get_neighbors(self):
        return True
    
    def display(self,Ids = False):
        # vis frem krystall

        # vis Id-er dersom Ids = True
        # sykle gjennom forskjellige lag
        return True

    def hamiltonian(self):
        # lag hamiltonian
        # fa en liste over alle elektroner
        
        return True
        
obj = name()

for num in obj:
    print(num)

obj.data.append(7)

for num in obj:
    print(num)




N = 2
# make a 1x1-crystal object
crystal = Crystal(N)

# inside crystal
crystal.add_nuclei()
crystal.add_orbitals()


crystal.perturb([2,2,2],True,[0.15,0.15,0.15])

ø# output hamiltonian

k = array([0.,0.,0.])
H = crystal.hamiltonian(k)
