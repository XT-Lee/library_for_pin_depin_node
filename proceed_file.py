import numpy as np
import os
#import points_analysis_2D as pa
R"""
CLASS list:
    proceed_gsd_file:
    proceed_exp_file:
"""
class proceed_file_shell:
    R"""
    """
    def __init__(self):
        pass
    
    def create_folder(self,prefix='/home/remote/Downloads/',folder='folder_new'):
        R"""
        introduction:
            check and create a new folder.
        input:
            prefix: must end with '/'
            folder: must not end with '/'
        output:
            prefix_new: end with '/'
        """
        folder_name=prefix+folder
        #check if the folder exists
        isExists=os.path.exists(folder_name)
        if isExists:
            pass
        else:
            os.makedirs(folder_name)
        prefix_new=folder_name+'/'
        return prefix_new

    def create_prefix_results(self,account='remote',simu_index=0,seed=9):
        R"""
        input:
            account: (string);
            simu_index: (int)index;
            seed: (int)0-9
            io_only: just return results, not proceeding data.
        return:
            prefix_simu: (str)directory end with '/';
            str_simu_index: (str)'index_seed' for example.
        """
        str_simu_index = str(int(simu_index))+'_'+str(seed)
        prefix0 = '/home/'+account+'/Downloads/'
        prefix_simu = self.create_folder(prefix0,str_simu_index)

        return prefix_simu,str_simu_index

class proceed_gsd_file:
    R"""
    Introduction:
        the class is designed to preproceed the motion of 2D points(.gsd file), 
        to analyze dynamic properties.

    Parameters:
        filename_gsd: a path to gsd file;
        trajectory: trajectory data read from a hoomd gsd file;
        num_of_frames: the number of frames in the trajectory;

    Methods:
        open_gsd:
        read_a_frame:
        get_displacement_field: draw the displacements of final state from initial state.

    
    Examples:
    """
    
    def __init__(self,filename_gsd_seed=None, account="tplab",simu_index=None,seed=None):
        #load positions of particles
        if simu_index is None:
            if filename_gsd_seed is None:
                    print("Error: input a correct path to gsd file please!\n")
            else:
                self.filename_gsd = filename_gsd_seed
                #filename_gsd="/home/tplab/hoomd-examples_0/trajectory_auto619.gsd"
                #self.data_gsd = np.loadtxt(self.filename_gsd)
                self.__open_gsd()

                prefix_gsd = '/home/'+account+'/hoomd-examples_0/trajectory_auto'
                simu_index = filename_gsd_seed.strip(prefix_gsd)
                id=simu_index.index('_')
                self.simu_index = simu_index[0:id]
        else :
            self.simu_index = simu_index
            if not seed is None:
                simu_index = str(int(simu_index))+'_'+str(int(seed))
            prefix_gsd = '/home/'+account+'/hoomd-examples_0/trajectory_auto'
            postfix_gsd = '.gsd'
            self.filename_gsd = prefix_gsd+str(simu_index)+postfix_gsd
            self.__open_gsd()
        
        self.box = self.trajectory.read_frame(-1).configuration.box

    def __open_gsd(self):
        import gsd.hoomd
        self.trajectory=gsd.hoomd.open(self.filename_gsd)#open a gsd file
        self.num_of_frames=len(self.trajectory)
        
    def read_a_frame(self,frame_num):
        snap=self.trajectory.read_frame(frame_num)#take a snapshot of the N-th frame
        positions=snap.particles.position[:,0:2]#just record [x,y] ignoring z
        #self.N = self.snap.particles.N
        return positions
        
    def get_trajectory_data(self,save_prefix = None):
        R"""
        introduction:
            transform gsd file into an array [Nframes,Nparticles,3],
            recording the trajectory of particles.
        input:
            gsd_file
        return:
            txyz [Nframes,Nparticles,3] or
            (npy file)[Nframes,Nparticles,3]
        """
        frame = 0
        snapi = self.trajectory.read_frame(frame)
        pos_list = np.zeros([self.num_of_frames,snapi.particles.N,3])#gsd_data.trajectory[0].particles.N,
        while frame < self.num_of_frames:
            pos_list[frame] = self.trajectory.read_frame(frame).particles.position
            #print(self.trajectory.read_frame(iframe).configuration.box)
            frame = frame + 1
        
        self.txyz = pos_list

        if not save_prefix is None:
            file_txyz_npy = save_prefix+'txyz'
            np.save(file = file_txyz_npy,arr = self.txyz)
        
    def get_trajectory_stable_data(self,save_prefix = None):
        R"""
        introduction:
            transform trajectory data from simulation with periodic boundary condition 
            into trajectories of which particles never move across boundary(box).
        return:
            txyz_stable: (N_frames,N_particles,3)

        """
        #dedrift?
        frames,particles,dimensions=self.txyz.shape
        if hasattr(self,'box'):
            #print(locals())#local variable not of class
            self.dtxyz = self.txyz[1:,:,:] - self.txyz[:frames-1,:,:]
            #cross is true
            list_crossleft = self.dtxyz[:,:,0] > 0.9*self.box[0]
            list_crossbottom = self.dtxyz[:,:,1] > 0.9*self.box[1]
            list_crossright = self.dtxyz[:,:,0] < -0.9*self.box[0]
            list_crosstop = self.dtxyz[:,:,1] < -0.9*self.box[1]
            #mark all the frames where cross event occur as True
            list_crossx = np.logical_or(list_crossleft,list_crossright)
            list_crossy = np.logical_or(list_crossbottom,list_crosstop)
            list_cross = np.logical_or(list_crossx,list_crossy)
            #mark all the particles who have experienced cross event as True
            list_cross_true = np.array(list_cross[0,:]) 
            #list_cross_true = list_cross_true[0]#remove empty extra dimension
            #print(list_cross_true.shape)
            i=0
            while i<particles:
                list_cross_true[i] = np.max(list_cross[:,i])
                i = i + 1
            list_stable_id = np.where(list_cross_true[:]==False)
            list_stable_id = list_stable_id[0]#remove empty extra dimension
            #print(list_stable_id.shape)
            
            self.txyz_stable = self.txyz[:,list_stable_id,:]
            self.particles = list_stable_id.shape[0]

            if not save_prefix is None:
                file_txyz_npy = save_prefix+'txyz_stable'
                np.save(file = file_txyz_npy,arr = self.txyz_stable)
    
    def get_trajectory_data_large(self,save_prefix = None):
        R"""
        introduction:
            transform gsd file into an array [Nframes,Nparticles,3],
            recording the trajectory of particles.
        input:
            gsd_file
        return:
            txyz [Nframes,Nparticles,3] or
            (npy file)[Nframes,Nparticles,3]
        """
        frame = 0
        snapi = self.trajectory.read_frame(frame)
        pos_list = np.zeros([int(self.num_of_frames/10+1),snapi.particles.N,3])#gsd_data.trajectory[0].particles.N,
        i=0
        while frame < self.num_of_frames:
            pos_list[i] = self.trajectory.read_frame(frame).particles.position
            #print(self.trajectory.read_frame(iframe).configuration.box)
            frame = frame + 10
            i = i+1
        
        self.txyz = pos_list

        if not save_prefix is None:
            file_txyz_npy = save_prefix+'txyz'
            np.save(file = file_txyz_npy,arr = self.txyz)

class proceed_exp_file:
    R"""
    see particle_tracking.py to get trajectories of particles from a video.
    """
    pass
