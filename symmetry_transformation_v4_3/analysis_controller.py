import points_analysis_2D as pa
import gsd.hoomd
import numpy as np
class from_setup_to_a_gsd:
    def __init__(self):
        pass

    def get_cn3_from_a_gsd(self):
        self.gsd_data = gsd.hoomd.open(self.input_file_gsd)
        self.gsd_data[-1]

    def set_file_parameters(self,simu_index,seed):
        self.prefix_read =  "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
        self.prefix_write = "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
        self.seed = seed
        self.simu_index = simu_index
        self.input_file_gsd = self.prefix_write+'trajectory_auto'+str(int(self.simu_index))+'_'+str(int(self.seed))+'.gsd'
        #self.input_file_gsd = self.prefix_read+'particle_and_trap.gsd'
        #self.output_file_gsd = self.prefix_write+'trajectory_auto'+str(int(self.simu_index))+'_'+str(int(self.seed))+'.gsd'
        self.snap_period = 1000
        self.gsd_data = gsd.hoomd.open(self.input_file_gsd)

    def file_to_data(self):
        self.gsd_data = gsd.hoomd.open(self.input_file_gsd)
        self.gsd_data[0]


class from_a_gsd_to_data:
    def __init__(self):
        pass

    def get_cn3_from_a_gsd_frame(self,last_frame):
        #get the data
        #last_frame = gsd_data[-1]
        xy_particles_traps = last_frame.particles.position[:,:2]

        ids = np.array(last_frame.particles.typeid)
        list_p = ids == 0
        list_t = ids == 1
        points = xy_particles_traps[list_p]
        traps = xy_particles_traps[list_t]
        obj_of_simu_index = pa.static_points_analysis_2d(points,hide_figure=False)
        tune_dis = 2.4#lattice_a*lcr?
        xmax = max(traps[:,0]) - tune_dis
        ymax = max(traps[:,1]) - tune_dis
        xmin = min(traps[:,0]) + tune_dis
        ymin = min(traps[:,1]) + tune_dis
        obj_of_simu_index.cut_edge_of_positions_by_xylimit(xmin,xmax,ymin,ymax)
        obj_of_simu_index.get_coordination_number_conditional()
        ccn=obj_of_simu_index.count_coordination_ratio
        #print(ccn[3])

        """import matplotlib.pyplot as plt
        fig,ax = plt.subplots()
        ax.scatter(points[:,0],points[:,1],color='k')#
        ax.scatter(traps[:,0],traps[:,1],color='r',marker = 'x')#
        fence = np.array([[xmin,ymin],[xmin,ymax],[xmax,ymax],[xmax,ymin],[xmin,ymin]])
        ax.plot(fence[:,0],fence[:,1])
        #ax.scatter(dula[:,0],dula[:,1],facecolors='none',edgecolors='k')#,marker = 'x'
        ax.set_xlabel('x label')  # Add an x-label to the axes.
        ax.set_ylabel('y label')  # Add a y-label to the axes.
        ax.set_title("Simple Plot")  # Add a title to the axes
        ax.set_aspect('equal','box')
        plt.show()
        plt.close('all')"""

        return ccn[3]

