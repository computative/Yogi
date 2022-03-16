from numpy import *
import matplotlib.pyplot as plt
from param import *
from crystal import *

N = 2
# make a 1x1-crystal object
crystal = Crystal(N)
# distributes nuclei
#crystal.add_nuclei()
# shift an atom
crystal.perturb([2,2,2],True,[0.15,0.15,0.15])
# add orbitals
#crystal.add_orbitals()
# output hamiltonian

k = array([0.,0.,0.])
H = crystal.hamiltonian(k)

from matplotlib.pyplot import *



#imshow(absolute(H))
#show()
#exit()


coords = crystal.coords
states = crystal.states

####################### endres
_coords = coords[::5]

E,C = linalg.eigh(H)


H = H.reshape(N**3*2,5,N**3*2,5).sum(axis=(1,3))

for omega in [E[0]*0.9,0]:
    # greens function is to be written in a 1s 1p^3 2s basis
    R = zeros((N**3*2*5,N**3*2*5),dtype=complex128)

    # set omega to be close to energy eigenvalue
    
    for i in range(N**3*2*5):
        for j in range(N**3*2*5):
            _entry = [ C[i,k]*conj(C[j,k])/(omega - E[k]) for k in range(N**3*2*5)]
            R[i,j] = sum(_entry)

    R = R.reshape(N**3*2,5,N**3*2,5).sum(axis=(1,3)).real

    print("omega =",omega)
    for obj,label in zip([R.real, H.real],["Green","Hamiltonian"]):
        #print(label)
        #print(obj.diagonal())
        obj = ( obj - obj.min() )/(obj.max()-obj.min())/N**2
        # endres
        obj = obj/2
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # all the scatter-plots of individual dots
        s = [None for r in obj.diagonal()]
        alphas = array(obj.diagonal())

        for i, alpha in enumerate(alphas):
            s[i] = ax.scatter(_coords[i,0],_coords[i,1],_coords[i,2], alpha = alpha, linewidth=0,  c="k")
            s[i].set_edgecolors = s[i].set_facecolors = lambda *args:None

        # plotting of the lines
        alphas = obj
        #print(alphas)
        #exit()
        for i in range(len(alphas)):
            for j in range(i,len(alphas)):
                # endres
                if i == j:
                    continue                
                #if linalg.norm(_coords[i]-_coords[j]) == 0:
                #    continue
                xs = [_coords[i,0],_coords[j,0]]
                ys = [_coords[i,1],_coords[j,1]]
                zs = [_coords[i,2],_coords[j,2]]
                ax.plot(xs,ys,zs, alpha = alphas[i,j], linewidth=linalg.norm(_coords[i]-_coords[j])**(-1),  c="k")
        plt.show()
