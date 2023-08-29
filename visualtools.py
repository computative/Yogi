import numpy as np
from vector import c

class VisualTools:

    def line(self, idx):
        plt.plot([idx,idx],[-13.5,12.15], "k",linewidth=0.6)

    @staticmethod
    def reciprocal(a1,a2,a3):
        A = np.array([a1,a2,a3]) # this is the transpose of [a1 a2 a3]
        B = 2*np.pi*np.linalg.inv(A) # the columns are reciprocal vec.
        return B[:,0], B[:,1], B[:,2]   


    @staticmethod
    def bands(bands):
        import matplotlib.pyplot as plt 

        for band in bands:
            plt.plot(band, "k", linewidth=1)

        plt.title("Silicon crystal momentum vs energy", fontsize=11)
        plt.xlabel("Crystal momentum")
        plt.ylabel("Energy (eV)")
        plt.show()

    @staticmethod
    def kpath(a,b,_n):
        r = (b-a)#/np.linalg.norm(b-a)
        path = []
        for x in np.linspace(0, 1, int(_n)):
            path.append(a + x*r)
        return path

    @staticmethod
    def kpts(path, n):
        k = []
        if type(n) != list :
            n = np.ones(len(path)-1)*n

        for i in range(len(path)-1):
            a = np.array(path[i])
            b = np.array(path[i+1])
            _n = n[i]
            k.extend( VisualTools.kpath(a, b, _n) )

        past = k[0]
        k_norepeated = [k[0]]
        for item in k[1:]:
            if (item != past).any():
                k_norepeated.append(item)
            past = item
        return k_norepeated

    @staticmethod
    def fcc_sympts(b1,b2,b3):
        return {
            "Gamma": c(0,0,0),
            "X": b2,
            "L": 0.5*(b1 + b2 + b3),
            "W": 0.5*b1 + b2,
            "U": 0.25*( b1 + 4*b2 + b3 ),
            "K": 0.75*(b1 + b2)
        }

    @staticmethod
    def fcc_sympts2(a):
        return {
            "Gamma": c(0,0,0), 
            "X": c(0, 2*np.pi/a,0) ,
            "L": c(np.pi/a,np.pi/a,np.pi/a),
            "W": c(np.pi/a,2*np.pi/a,0),
            "U": c(np.pi/(2*a),2*np.pi/a,np.pi/(2*a)),
            "K": c(3*np.pi/(2*a),3*np.pi/(2*a),0)
            }

    @staticmethod
    def diamond_sympts2(a):
        return {
            "Gamma": c(0,0,0), 
            "X": c(2*np.pi/a, 0,0) ,
            "L": c(np.pi/a,np.pi/a,np.pi/a),
            "W": c(2*np.pi/a,np.pi/a,0),
            "U": c(2*np.pi/a,np.pi/(2*a),np.pi/(2*a)),
            "K": c(3*np.pi/(2*a),3*np.pi/(2*a),0)
            }
    @staticmethod
    def matrix_formatting():
        np.set_printoptions(threshold=np.inf)
        np.set_printoptions(linewidth=np.inf)
        np.set_printoptions(suppress=True)
