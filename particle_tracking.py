#particle tracking

import matplotlib as mpl, matplotlib.pyplot as plt#, matplotlib.image as imread
import numpy as np
#import pandas as pd
#from pandas import DataFrame, Series  # for convenience
import pims
import trackpy as tp

class particle_track:
    R"""
        algorithm from Eric Weeks
            https://physics.emory.edu/faculty/weeks/idl/
        library from github.soft-matter.trackpy
            http://soft-matter.github.io/trackpy/v0.5.0/tutorial/walkthrough.html

        f collumn name: y 	x 	mass 	size 	ecc 	signal 	raw_mass 	ep 	frame
        """
    def __init__(self):
        pass
    
    def single_frame_particle_tracking(self,filename,D=19,minmass=1000,calibration=False):
        R"""
        introduction:
            algorithm from Eric Weeks
                https://physics.emory.edu/faculty/weeks/idl/
            library from github.soft-matter.trackpy
                http://soft-matter.github.io/trackpy/v0.5.0/tutorial/walkthrough.html
        
        instruction:
            D should be the the diameter of dark ring of the smallest particle in image.
            minmass should be the maximum mass of the darker cluster.
            separation should be 0.9*D if dimer particles occur.

            f collumn name: y 	x 	mass 	size 	ecc 	signal 	raw_mass 	ep 	frame
                DataFrame([x, y, mass, size, ecc, signal, raw_mass])
                where "x, y" are appropriate to the dimensionality of the image,
                mass means total integrated brightness of the blob,
                size means the radius of gyration of its Gaussian-like profile,
                ecc is its eccentricity (0 is circular),
                and raw_mass is the total integrated brightness in raw_image.
                
        example:
            import particle_tracking as pt
            filename= '/home/tplab/xiaotian_file/data/20220924/DefaultVideo_5.tif'
            track = pt.particle_track()
            track.single_particle_tracking(filename,i=12,D=31,minmass=5000)
            track.multiple_particle_tracking(filename) 
        """
        f0 = plt.imread(filename)
        f0 = f0[:,:,0]#3 channel to 1 channel
        sz=np.shape(f0)
        #print(sz)
        #plt.imshow(frames[0])
        #diameter of particles should includes dark ring!
        self.f = tp.locate(f0, D, minmass= minmass,invert=False,separation=0.9*D)#diameter must be odd in pixels
        #f is feature.
        #collumn name: y 	x 	mass 	size 	ecc 	signal 	raw_mass 	ep 	frame
        if calibration:
            #print(self.f.head())#print(self.f['x'])
            plt.figure(1)
            plt.scatter(self.f['mass'],self.f['size'])
            plt.xlabel('mass')
            plt.ylabel('size')

            plt.figure(2)
            plt.scatter(self.f['mass'],self.f['ecc'])
            plt.xlabel('mass')
            plt.ylabel('ecc')
            #plt.show()
            
            #selected_feature
            #sf=self.f[self.f['mass'] > 6000 ]
            #sf2=sf[sf['mass'] < 3000]

            plt.figure(3)
            tp.annotate(self.f, f0)#,imshow_style=
            #plt.imsave(png_filename)#savefig(png_filename)   miss array
            #plt.close()
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
        self.xy = np.array(self.f[['x','y']])# unit: pixel
        self.Diameter = D
        self.minmass = minmass

    def single_stack_particle_tracking(self,filename,i=7,D=35,minmass=7000):
        R"""
        from github.soft-matter.trackpy
        http://soft-matter.github.io/trackpy/v0.5.0/tutorial/walkthrough.html

        example:
            import particle_tracking as pt
            filename= '/home/tplab/xiaotian_file/data/20220924/DefaultVideo_5.tif'
            track = pt.particle_track()
            track.single_particle_tracking(filename,i=12,D=31,minmass=5000)
            track.multiple_particle_tracking(filename) 
        """
        
        
        # change the following to %matplotlib notebook for interactive plotting
        #%matplotlib inline

        # Optionally, tweak styles.
        #mpl.rc('figure',  figsize=(10, 5))
        #mpl.rc('image', cmap='gray')

        #read a tiff file
        
        frames = pims.open(filename)
        f0= frames[i]
        sz=np.shape(f0)
        #print(sz)
        #plt.imshow(frames[0])
        #diameter of particles should includes dark edge! 35
        f = tp.locate(frames[i], D, minmass= minmass,invert=False)#diameter must be odd in pixels
        
        #f means feature.
        #collumn name: y 	x 	mass 	size 	ecc 	signal 	raw_mass 	ep 	frame
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
        
        self.Diameter = D
        self.minmass = minmass

    def multiple_particle_tracking(self,filename):
        frames = pims.open(filename)
        """
        sz=np.shape(frames)
        f=[]
        #for i in range(sz[0]):
        for i in range(3):
            fi = tp.locate(frames[i], 35, minmass= 2000,invert=False)
            f.append(fi)
        """
        f = tp.batch(frames[:100],35,minmass=5000)# self.Diameter self.minmass what the fuck corrupted tag list
        t = tp.link(f, 5, memory=3)
        t1 = tp.filter_stubs(t)
        # Compare the number of particles in the unfiltered and filtered data.
        print('Before:', t['particle'].nunique())
        print('After:', t1['particle'].nunique())
        plt.figure()
        tp.plot_traj(t);

    def track_bright_field(self,filename):
        R"""
        import particle_tracking as pt
        filename= '/home/tplab/xiaotian_file/data/20220924/DefaultVideo_5.tif'
        track = pt.particle_track()
        track.track_bright_field(filename)
        """
        frames = pims.open(filename)
        micron_per_pixel = 3.0/32.0
        feature_diameter = 2.1 # um
        radius = int(np.round(feature_diameter/2.0/micron_per_pixel))
        if radius % 2 == 0:
            radius += 1
        print('Using a radius of {:d} px'.format(radius))

        # we use a slightly larger radius
        f_locate = tp.locate(frames[0], radius+2, minmass=2000)
        tp.annotate(f_locate, frames[0], plot_style={'markersize': radius});

        plt.figure()
        plt.hist(f_locate['mass'], bins=20)#here mass is the brightness of a particle
        plt.show()


    

