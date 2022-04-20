import numpy as np

class PlotTools:

    def line(self, idx):
            plt.plot([idx,idx],[-13.5,12.15], "k",linewidth=0.6)


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
            k.extend( PlotTools.kpath(a, b, _n)[:-2] )
        return np.array(k)

    
