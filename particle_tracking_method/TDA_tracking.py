

R"""
Introduction:
    in order to design an algorithm to recognize particle without manually set parameters, 
    which can hence output particle density and distance immediately after taking photos.
    the program should be set in Windows system.
Process:
    brightness Ranking LINK particle groups

    let the num of groups increase till the 1st maximum peak, see these groups as particle.

    calculate their center as particle center
MATLAB resources:
    https://www.mathworks.com/matlabcentral/fileexchange/47009-topological-data-analysis-learning-code
Python resources:
 



"""
import matplotlib as mpl, matplotlib.pyplot as plt, matplotlib.image as imread
import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import pims
import trackpy as tp

class particle_track_TDA:
    def __init__(self):
        pass

    def single_frame(self,filename,i=1):
        R"""
        from github.soft-matter.trackpy
        http://soft-matter.github.io/trackpy/v0.5.0/tutorial/walkthrough.html

        example:
            import TDA_tracking as pt
            filename='/home/tplab/xiaotian_file/data/20220924/DefaultVideo_5.tif'#'/home/remote/Downloads/3-times-size.png'# 
            track = pt.particle_track_TDA()
            track.single_frame(filename,i=12)
        """
        
        
        # change the following to %matplotlib notebook for interactive plotting
        #%matplotlib inline

        # Optionally, tweak styles.
        #mpl.rc('figure',  figsize=(10, 5))
        #mpl.rc('image', cmap='gray')

        #read a tiff file
        
        frames = pims.open(filename)
        f0= frames[i]
        f0[:]= 255-f0[:] # reverse brightness
        sz=np.shape(f0)
        print(sz)
        plt.imshow(f0,'Greys')
        plt.show()
        #diameter of particles should includes dark edge! 35
        
        """
        tp.annotate(f, frames[i])
        print(f.head())
        #plt.imshow(frames[0])

        plt.figure()
        plt.hist(f['mass'], bins=20)#here mass is the brightness of a particle
        plt.show()
        # Optionally, label the axes.
        #ax.set(xlabel='mass', ylabel='count');
        """
    def tda_process(self):
        R"""
        Introduction:
        Libraries:
            scikit-tda: https://scikit-tda.org/libraries.html
            giotto-learn
        Knowledge:
            https://zhuanlan.zhihu.com/p/41292296

        """
        #import scikit
        import ripser
        import numpy as np
        from ripser import ripser
        from persim import plot_diagrams
        import tadasets
        #tutorial
        #https://ripser.scikit-tda.org/en/latest/notebooks/Representative%20Cocycles.html
        np.random.seed(9)
        x = tadasets.dsphere(n=12, d=1, noise=0.1)
        #data = np.random.random((100,2))

        plt.figure(1)
        plt.scatter(x[:, 0], x[:, 1])
        plt.axis('equal')
        #plt.show()

        plt.figure(2)
        result = ripser(x, coeff=17, do_cocycles=True)
        diagrams = result['dgms']
        cocycles = result['cocycles']
        D = result['dperm2all']

        dgm1 = diagrams[1]
        idx = np.argmax(dgm1[:, 1] - dgm1[:, 0])
        plot_diagrams(diagrams, show = False)
        plt.scatter(dgm1[idx, 0], dgm1[idx, 1], 20, 'k', 'x')
        plt.title("Max 1D birth = %.3g, death = %.3g"%(dgm1[idx, 0], dgm1[idx, 1]))
        
        plt.figure(3)
        cocycle = cocycles[1][idx]
        thresh = dgm1[idx, 1] #Project cocycle onto edges less than or equal to death time
        self.plotCocycle2D(D, x, cocycle, thresh)
        plt.title("1-Form Thresh=%g"%thresh)
        plt.show()

    def drawLineColored(self, X, C):
        for i in range(X.shape[0]-1):
            plt.plot(X[i:i+2, 0], X[i:i+2, 1], c=C[i, :])#, lineWidth = 3

    def plotCocycle2D(self, D, X, cocycle, thresh):
        """
        Given a 2D point cloud X, display a cocycle projected
        onto edges under a given threshold "thresh"
        """
        #Plot all edges under the threshold
        N = X.shape[0]
        t = np.linspace(0, 1, 10)
        c = plt.get_cmap('Greys')
        C = c(np.array(np.round(np.linspace(0, 255, len(t))), dtype=np.int32))
        C = C[:, 0:3]

        for i in range(N):
            for j in range(N):
                if D[i, j] <= thresh:
                    Y = np.zeros((len(t), 2))
                    Y[:, 0] = X[i, 0] + t*(X[j, 0] - X[i, 0])
                    Y[:, 1] = X[i, 1] + t*(X[j, 1] - X[i, 1])
                    self.drawLineColored(Y, C)
        #Plot cocycle projected to edges under the chosen threshold
        for k in range(cocycle.shape[0]):
            [i, j, val] = cocycle[k, :]
            if D[i, j] <= thresh:
                [i, j] = [min(i, j), max(i, j)]
                a = 0.5*(X[i, :] + X[j, :])
                plt.text(a[0], a[1], '%g'%val, color='b')
        #Plot vertex labels
        for i in range(N):
            plt.text(X[i, 0], X[i, 1], '%i'%i, color='r')
        plt.axis('equal')
      