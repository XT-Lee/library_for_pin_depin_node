import points_analysis_2D as pa
import gsd.hoomd
import numpy as np
class get_a_gsd_from_setup:
    def __init__(self):
        pass

    def set_file_parameters(self,simu_index,seed):
        self.prefix_read =  "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_1/"
        self.prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_1/"
        self.seed = seed
        self.simu_index = simu_index
        self.input_file_gsd = self.prefix_write+'trajectory_auto'+str(int(self.simu_index))+'_'+str(int(self.seed))+'.gsd'
        #self.input_file_gsd = self.prefix_read+'particle_and_trap.gsd'
        #self.output_file_gsd = self.prefix_write+'trajectory_auto'+str(int(self.simu_index))+'_'+str(int(self.seed))+'.gsd'
        self.snap_period = 1000
        self.gsd_data = gsd.hoomd.open(self.input_file_gsd)

    def get_gsd_data_from_file(self):
        self.gsd_data = gsd.hoomd.open(self.input_file_gsd)
    
    def get_gsd_data_from_filename(self,input_file_gsd):
        self.gsd_data = gsd.hoomd.open(input_file_gsd)

class get_gsds_from_csv:
    def __init__(self):
        pass

    def get_record_from_csv(self,csv_filename):
        R"""
        table_name = 'pin_hex_to_honeycomb_klt_2m_gauss'
        | simu_index | seed | lcr  | trap_gauss_epsilon | temperature |
        list_simu = [colums of table]
        """
        import symmetry_transformation_v4_3.simulation_core as sc
        import pandas as pd
        record = pd.read_csv(csv_filename)
        simu_index = record['simu_index'].values
        seed = record['seed'].values
        n_simu = len(seed)
        gsds = []
        for i in range(n_simu):
            sct = sc.simulation_core_traps(simu_index[i],seed[i])
            gsds.append(sct.output_file_gsd)#record[i].append(sct.output_file_gsd) 
            #print(record[i])
            del sct
        self.record = record
        return gsds

    def get_last_frame_from_gsd(self,gsd_filename):
        gsd_data = gsd.hoomd.open(gsd_filename)
        return gsd_data[-1]

class get_data_from_a_gsd_frame_with_traps:
    def __init__(self,last_frame):
        self.last_frame = last_frame
        xy_particles_traps = last_frame.particles.position[:,:2]
        ids = np.array(last_frame.particles.typeid)
        list_p = ids == 0
        list_t = ids == 1
        self.points = xy_particles_traps[list_p]
        self.traps = xy_particles_traps[list_t]

    def get_cn_k_from_a_gsd_frame(self,tune_dis = 2.4,k=None):
        points = self.points
        traps = self.traps
        obj_of_simu_index = pa.static_points_analysis_2d(points)#,hide_figure=False
        #tune_dis = 2.4#lattice_a*lcr?
        xmax = max(traps[:,0]) - tune_dis
        ymax = max(traps[:,1]) - tune_dis
        xmin = min(traps[:,0]) + tune_dis
        ymin = min(traps[:,1]) + tune_dis
        obj_of_simu_index.cut_edge_of_positions_by_xylimit(xmin,xmax,ymin,ymax)
        obj_of_simu_index.get_coordination_number_conditional(tune_dis)
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
        if k is None:
            return ccn
        else:
            return ccn[k]
    
    def get_bonds_png_from_a_gsd_frame(self,png_filename, tune_dis=2.4):
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots()
        p2d = pa.static_points_analysis_2d(self.points,hide_figure=False)#dis_edge_cut=
        
        p2d.get_first_minima_bond_length_distribution(
            lattice_constant=tune_dis, png_filename='bond_hist.png')
        #draw bonds selected
        bpm = pa.bond_plot_module(fig,ax)#
        bpm.restrict_axis_property_relative('(sigma)')#'(sigma)',hide_axis=True
        list_bond_index = bpm.get_bonds_with_conditional_bond_length(
            p2d.bond_length,[1.2,p2d.bond_first_minima_left])
        
        bpm.draw_points_with_given_bonds(self.points,list_bond_index,bond_color='k',particle_color='k')#p2d.bond_length[:,:2].astype(int)
        bpm.plot_traps(self.traps)
        #bpm.restrict_axis_limitation([-10,10],[-10,10])#[-20,0],[-5,15]
        bpm.save_figure(png_filename)
        del bpm
    
    def get_given_bonds_png_from_a_gsd_frame(self,png_filename,bond_first_minima_left):
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots()
        p2d = pa.static_points_analysis_2d(self.points,hide_figure=False)#dis_edge_cut=
        
        p2d.get_first_minima_bond_length_distribution(png_filename='bond_hist.png')
        #draw bonds selected
        bpm = pa.bond_plot_module(fig,ax)#
        bpm.restrict_axis_property_relative(hide_axis=True)#'(sigma)'
        bpm.restrict_axis_limitation([-10,10],[-10,10])
        list_bond_index = bpm.get_bonds_with_conditional_bond_length(p2d.bond_length,[1.2,bond_first_minima_left])
        
        bpm.draw_points_with_given_bonds(self.points,list_bond_index,bond_color='k',particle_color='k')#p2d.bond_length[:,:2].astype(int),,particle_size=20
        bpm.plot_traps(self.traps)
        #bpm.restrict_axis_limitation([-10,10],[-10,10])
        bpm.save_figure(png_filename)
        del bpm

class get_data_from_a_gsd_frame:
    def __init__(self,last_frame):
        self.last_frame = last_frame
        xy_particles = last_frame.particles.position[:,:2]
        self.points = xy_particles

    def get_cn_k_from_a_gsd_frame(self,tune_dis = 2.4,k=3):
        points = self.points
        obj_of_simu_index = pa.static_points_analysis_2d(points)#,hide_figure=False
        #tune_dis = 2.4#lattice_a*lcr?
        xmax = max(points[:,0]) - tune_dis
        ymax = max(points[:,1]) - tune_dis
        xmin = min(points[:,0]) + tune_dis
        ymin = min(points[:,1]) + tune_dis
        obj_of_simu_index.cut_edge_of_positions_by_xylimit(xmin,xmax,ymin,ymax)
        obj_of_simu_index.get_coordination_number_conditional(lattice_constant=tune_dis)
        ccn=obj_of_simu_index.count_coordination_ratio
        #print(ccn[3])
        return ccn[k]
    
    def get_grs_from_a_gsd_frame(self,tune_dis = 2.4,png_filename=None,normalize_r0=None):
        points = self.points
        obj_of_simu_index = pa.static_points_analysis_2d(points,dis_edge_cut=tune_dis)#,hide_figure=False
        box = self.last_frame.configuration.box[:2]
        obj_of_simu_index.get_radial_distribution_function(box,r_max=tune_dis*4,dr=tune_dis/11,png_filename=png_filename,normalize_r0=normalize_r0)#