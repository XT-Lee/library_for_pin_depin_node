import points_analysis_2D as pa
import gsd.hoomd
import numpy as np
class get_a_gsd_from_setup:
    def __init__(self):
        pass

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

    def get_gsd_data_from_filename(self):
        self.gsd_data = gsd.hoomd.open(self.input_file_gsd)

class get_gsds_from_mysql:
    def __init__(self):
        pass

    def get_record_from_sql_by_lcr(self,lcr1=0.81,table_name = 'pin_hex_to_honeycomb_klt_2m_gauss'):
        R"""
        table_name = 'pin_hex_to_honeycomb_klt_2m_gauss'
        | simu_index | seed | lcr  | trap_gauss_epsilon | temperature |
        list_simu = [colums of table]
        """
        import opertateOnMysql as osql
        #import data_retriever as dr
        #se = dr.search_engine_for_simulation_database()
        
        #se.search_single_simu_by_lcr_k(table_name,0.81)
        lcr_step = 0.0001
        lcr_min=lcr1 - 0.5*lcr_step
        lcr_max=lcr1 + 0.5*lcr_step
        con = ' where lcr >'+str(lcr_min)+' and lcr <'+str(lcr_max)
        list_simu = osql.getDataFromMysql(table_name=table_name,search_condition=con)
        #simu_index,seed,lcr,k,kT
        return list_simu
    
    def get_gsds_from_mysql_record(self,record):
        R"""
        table_name = 'pin_hex_to_honeycomb_klt_2m_gauss'
        | simu_index | seed | lcr  | trap_gauss_epsilon | temperature |
        list_simu = [colums of table]
        """
        import symmetry_transformation_v4_3.simulation_core as sc
        n_simu = len(record)
        gsds = []
        for i in range(n_simu):
            sct = sc.simulation_core_traps(record[i][0],record[i][1])
            gsds.append(sct.output_file_gsd)#record[i].append(sct.output_file_gsd) 
            #print(record[i])
        return gsds

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
        self.record = record
        return gsds

    def get_last_frame_from_gsd(self,gsd_filename):
        gsd_data = gsd.hoomd.open(gsd_filename)
        return gsd_data[-1]

class get_data_from_a_gsd_frame:
    def __init__(self,last_frame):
        self.last_frame = last_frame
        xy_particles_traps = last_frame.particles.position[:,:2]
        ids = np.array(last_frame.particles.typeid)
        list_p = ids == 0
        list_t = ids == 1
        self.points = xy_particles_traps[list_p]
        self.traps = xy_particles_traps[list_t]

    def get_cn3_from_a_gsd_frame(self,last_frame,tune_dis = 2.4):
        #get the data
        #last_frame = gsd_data[-1]
        xy_particles_traps = last_frame.particles.position[:,:2]
        ids = np.array(last_frame.particles.typeid)
        list_p = ids == 0
        list_t = ids == 1
        points = xy_particles_traps[list_p]
        traps = xy_particles_traps[list_t]
        obj_of_simu_index = pa.static_points_analysis_2d(points)#,hide_figure=False
        #tune_dis = 2.4#lattice_a*lcr?
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
    
    def get_bonds_from_a_gsd_frame(self,tune_dis = 2.4):
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots()
        p2d = pa.static_points_analysis_2d(self.points,hide_figure=False)#
        
        p2d.get_first_minima_bond_length_distribution(png_filename='bond_hist.png')
        #draw bonds selected
        bpm = pa.bond_plot_module(fig,ax)#
        bpm.restrict_axis_property_relative('(sigma)')
        list_bond_index = bpm.get_bonds_with_conditional_bond_length(p2d.bond_length,[2,p2d.bond_first_minima_left])
        
        bpm.draw_points_with_given_bonds(self.points,list_bond_index,bond_color='k',particle_color='k')#p2d.bond_length[:,:2].astype(int)
        bpm.plot_traps(self.traps)
        plt.show()
        """traps=self.traps
        #tune_dis = 2.4#lattice_a*lcr?
        xmax = max(traps[:,0]) - tune_dis
        ymax = max(traps[:,1]) - tune_dis
        xmin = min(traps[:,0]) + tune_dis
        ymin = min(traps[:,1]) + tune_dis
        p2d.cut_edge_of_positions_by_xylimit(xmin,xmax,ymin,ymax)"""
    
    def get_bonds_png_from_a_gsd_frame(self,png_filename):
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots()
        p2d = pa.static_points_analysis_2d(self.points,hide_figure=False)#
        
        p2d.get_first_minima_bond_length_distribution(png_filename='bond_hist.png')
        #draw bonds selected
        bpm = pa.bond_plot_module(fig,ax)#
        bpm.restrict_axis_property_relative('(sigma)')
        list_bond_index = bpm.get_bonds_with_conditional_bond_length(p2d.bond_length,[2,p2d.bond_first_minima_left])
        
        bpm.draw_points_with_given_bonds(self.points,list_bond_index,bond_color='k',particle_color='k')#p2d.bond_length[:,:2].astype(int)
        bpm.plot_traps(self.traps)
        bpm.save_figure(png_filename)
        del bpm

