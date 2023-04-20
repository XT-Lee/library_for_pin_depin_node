import numpy as np
import matplotlib.pyplot as plt
from prometheus_client import Summary
from scipy.spatial import Voronoi
from scipy.spatial import  Delaunay
from scipy.spatial import  distance
import pandas as pd
import matplotlib 
import os
import freud

R"""
CLASS list:
    static_points_analysis_2d: old_class_name: PointsAnalysis2D
    proceed_gsd_file:
    proceed_exp_file:
    dynamic_points_analysis_2d: old_class_name: msd
"""
class static_points_analysis_2d:#old_class_name: PointsAnalysis2D
    R"""
    Introduction:
        the class is designed to analyze the 
        geometric properties(listed as follows) of 2D points.
    
    Parameters:
        points: n rows of (x,y), coordinations of given particles.
        ---
        voronoi: scipy.spatial.Voronoi
        voronoi.ridge_length: num-of-ridges row of (ridge_length)
        voronoi.cut_ridge_length: a length that 
            if a ridge whose length is smaller than which, the ridge will be cut.
        voronoi.cluster_vertices: list the row indices of vertices that may form clusters.
        ---
        delaunay: scipy.spatial.Delaunay
        bond_length: n rows of (start_point_index, end_point_index, bond_length)
        cutted_bonds: n rows of (start_point_index, end_point_index)
        vertices_triangle: input row indices of vertices,
            and get row indices of triangles which cover these vertices.[x]
        count_coordination: record the count of the coordination number who has i neighbours [0n,1n,2n,...,9n]
            that is, 0n is the count of coordination number 0, 1n is the count of coordination number 1,
            5n is the count of coordination number 5,6n is the count of coordination number 6
        count_coordination_ratio:[x]count_coordination whose sum is normalized to 1.
        ---
        Psi_k: n rows of complex number, Bond Orientaional Order(k fold disclination) 
            for each particle.
        Psi_k_global: arithmetic mean of absolute value of Psi_k s.
        Psi_k_rate: A is the num of particles whose psi_k is larger than 0.9, while
            B is the num of all the particles. Psi_k_rate = A/B.

    Methods:
        points_analysis: get the parameters listed before
        get_ridge_length: 
        get_bond_length:
        get_first_minima_bond_length_distribution: 
        get_coordination_number_conditional:
        get_cluster_vertices:
        get_vertices_triangle:[x]
        get_bond_orientational_order: get bond orientational order, Psi_n, for each particle.
        get_bond_orientational_order_selected:
        get_neighbor_cloud: place all the neighbor particles into a plot 
                            with the same position of central particles
        ----
        draw_bonds_conditional_bond: given bond_length condition, draw bonds
        draw_bonds_conditional_ridge: given ridge_length condition, draw bonds
        draw_bond_length_distribution_and_first_minima
        ----
        print_benchmark: print (min,median,max) data for bond or ridge
        ----
        count_n_polygon_fast: vertices cluster method

    Examples:
        import points_analysis_2D
        obj_of_simu_index = points_analysis_2D.PointsAnalysis2D(filename=data_filename)
    """
    def __init__(self,points=None,filename=None,hide_figure = True):
        #load positions of particles
        if points is None:
            if filename is None:
                print("Error: input points or file please!\n")
            else:
                self.filename = filename#[x]
                #filename='/home/tplab/Downloads/index93.0'
                self.points = np.loadtxt(filename)
                self.points = self.points[:,0:2]
                self.basic_points_analysis()
        else :
            self.points = points[:,:2]
            self.basic_points_analysis()

        #not show figures
        if hide_figure:
            matplotlib.use(backend="agg")#Backend agg is non-interactive backend. Turning interactive mode off. 'QtAgg' is interactive mode
        #set basic parameters
        #self.prefix='/home/tplab/Downloads/'
            
    def basic_points_analysis(self):
        self.voronoi = Voronoi(self.points)
        self.delaunay = Delaunay(self.points)
        self.get_ridge_length()
        self.get_bond_length()
        #cut the edge of positions
        self.__cut_edge_of_positions()#effective lattice constant is given 3 defaultly

    def get_ridge_length(self):
        #where the ridge length of vertex-i & vertex-j is ridge_length[i,j]
        self.voronoi.ridge_vertices=np.array(self.voronoi.ridge_vertices)#选不了ridge_vertices,不是数组格式
        self.voronoi.ridge_length=distance.cdist(self.voronoi.vertices[self.voronoi.ridge_vertices[:,0],:],self.voronoi.vertices[self.voronoi.ridge_vertices[:,1],:],'euclidean')
        self.voronoi.ridge_length=self.voronoi.ridge_length.diagonal()

    def get_bond_length(self):
        #self.bond_length: (start_point_index,end_point_index,bond_length)
        shp=np.shape(self.voronoi.ridge_points)
        self.bond_length = np.zeros((shp[0],shp[1]+1))#(start_point,end_point,length)
        self.bond_length_temp=distance.cdist(self.voronoi.points[self.voronoi.ridge_points[:,0],:], self.voronoi.points[self.voronoi.ridge_points[:,1],:], 'euclidean')
        self.bond_length[:,0:2]=self.voronoi.ridge_points
        self.bond_length[:,2]=self.bond_length_temp.diagonal()
        #should I list point pair like small-to-large？
        del self.bond_length_temp #delete temp data

    def get_first_minima_bond_length_distribution(self,lattice_constant=3.0,hist_cutoff=2,
                                    method="global_compare_method", png_filename=None,x_unit='(1)'):
        R"""  
        Introduction:
            It's likely that bonds whose length are larger than 2A are not bonded. 
            Typical bond-lengths of honeycomb are A and 1.73A (bond_first_neighbour & bond_second_neighbour namely),
            hence the range=[0,2] times of lattice constant should be set.
        Parameters:
            lattice_constant: particle diameter or lattice constant is recommended(also the unit A in x-axis in hist)
            
            method:"local_compare_nethod" to get the 1st minima through comparing local minima,
                suitable for continuous peak-valley hist;
                "global_compare_method" to get the 1st minima through comparing all the peaks, 
                selecting the 1st main peak(ignoring tiny peak before which), 
                finding the 1st main minima just after the 1st main peak, 
                suitable systems with powerful perturbation.
            
            hist_cutoff: plot hist of bond_length till n times lattice_constant where n is hist_cutoff.
            
            png_filename="prefix/bond_hist_index1512"

            bond_first_minima_left:the upper limit of 1st neighbour bonds comming from the first minima 
                                    of bond length distribution, with sigma as lower limit.
            
        Warning:
            [x]What if there is no minima found? 
            count_1st_max / count_1st_minima > 10 used to check the effective of minima? 

        Examples:
            s = "/home/tplab/Downloads/"
            index_num = 1387
            index_name = "index"+str(int(index_num))
            fname = s+index_name
            bb = pa.PointsAnalysis2D(filename=fname)
            oname1 = s+"bond_hist_"+index_name
            bb.draw_bond_length_distribution_and_first_minima(png_filename=oname1)
            oname2 = s+"bond_plot_"+index_name
            bb.draw_bonds_conditional_bond(check=[0.9, bb.bond_first_minima_left], png_filename=oname2)
        """
        self.lattice_constant=lattice_constant
        
        #locate the 1st minima of bond-length distribution
        plt.figure()
        count_bins=plt.hist(self.bond_length[:,2]/self.lattice_constant,bins=20,range=[0,hist_cutoff])
        
        self._count=count_bins[0]
        self._bins=count_bins[1]
        self.bond_sorted=np.sort(self.bond_length[:,2]/self.lattice_constant)
        #find the minimum bin, then set the left side as bond_first_neighbour
        i_max=self._bins.size
        i=0
        if method=="local_compare_method":
            while i < i_max-3:#i start from 0, which have to -1;compares i,i+1 and i+2, which have to -2, hence -3
                #np.where( self.bond_sorted[:]>res.bins[i] & self.bond_sorted[:]< res.bins[i+1] ) 
                #print(self._count[i])
                if self._count[i] > self._count[i+1]:
                    if self._count[i+1] <= self._count[i+2]:
                        if self._bins[i+1] > 1/self.lattice_constant:
                            # bond_length should be larger than sigma of particle
                            i+=1
                            break
                i+=1

        elif method=="global_compare_method":
            self._count = np.array(self._count)
            self.count_sorted=np.sort(self._count)
            """
            print(self.count_sorted)
            print(self.count_sorted[-1])#max1
            print(self.count_sorted[-2])#max2
            """
            i_bins_for_count_peak_1 = np.where(self._count[:]==self.count_sorted[-1])
            i_bins_for_count_peak_2 = np.where(self._count[:]==self.count_sorted[-2])
            i_bins_for_count_peak_1 = i_bins_for_count_peak_1[0]
            i_bins_for_count_peak_2 = i_bins_for_count_peak_2[0]
            #if there are more than one bins share the same count, select the smallest bin number.
            if np.shape(i_bins_for_count_peak_1)[0]>1:
                i_bins_for_count_peak_1 = min(i_bins_for_count_peak_1)
            if np.shape(i_bins_for_count_peak_2)[0]>1:
                i_bins_for_count_peak_2 = min(i_bins_for_count_peak_2)
            i_bins_for_count_1st_peak = min(i_bins_for_count_peak_1,i_bins_for_count_peak_2)
            i_bins_for_count_1st_peak = int(i_bins_for_count_1st_peak)#force i be a simple data_type(int)
            """
            print(i_bins_for_count_peak_1) 
            print(i_bins_for_count_peak_2)
            print(i_bins_for_count_1st_peak) 
            """
            #list_xy = np.logical_and(list_x,list_y)
            #self.edge_cut_positions_list = np.where(==)
            i = i_bins_for_count_1st_peak
            while i < i_max-3:#i start from 0, which have to -1;compares i,i+1 and i+2, which have to -2, hence -3
                if self._count[i] > self._count[i+1]:
                    if self._count[i+1] <= self._count[i+2]:
                        if self._bins[i+1] > 1/self.lattice_constant:
                            # bond_length should be larger than sigma of particle
                            i+=1
                            break
                i+=1
        else:
            print("Err: input 'local_compare_method' or 'global_compare_method', please!")
            print("x")
            
        self.bond_first_minima_left=self._bins[i]
        self.bond_first_neighbour=self.bond_sorted[np.where(self.bond_sorted[:]<self.bond_first_minima_left)]
        
        self.bond_length_median = np.around( np.median(self.bond_first_neighbour),2)
        self.bond_length_std = np.around(np.std(self.bond_first_neighbour),2)

        if not png_filename  is None:
            #plot divide line
            plt.plot([self.bond_first_minima_left,self.bond_first_minima_left],[0,self._count[i]],linestyle='--')
            plt.title("1st peak:"+str( self.bond_length_median )+"+-"+str( self.bond_length_std )+x_unit)
            plt.xlabel("bond length "+x_unit)
            plt.ylabel("count (1)")
            #plt.show()
            #png_filename=prefix+"bond_hist_index"+str(index_num)
            plt.savefig(png_filename)
        plt.close()
        #let bond_first_minima_left be a absolute length, not normalized one
        self.bond_first_minima_left=self.bond_first_minima_left*self.lattice_constant
    
    def get_first_minima_ridge_length_distribution(self,hist_cutoff=5,
                                    method="median_minima_method", png_filename=None,x_unit='(1)'):
        R"""  
        Introduction:
            It's likely that bonds whose voronoi ridges are extremely small are not bonded. 
        Parameters:
            lattice_constant: particle diameter or lattice constant is recommended(also the unit A in x-axis in hist)
            
            method:"local_compare_method" to get the 1st minima through comparing local minima,
                suitable for continuous peak-valley hist; given here exit several local minima, the method has been given up.
                "global_compare_method" to get the 1st minima through comparing all the peaks, 
                selecting the 1st main peak(ignoring tiny peak before which), 
                finding the 1st main minima just after the 1st main peak, 
                suitable systems with powerful perturbation.
            
            hist_cutoff: plot hist of ridge_length till n times lattice_constant where n is hist_cutoff.
            
            png_filename="prefix/bond_hist_index1512"

            bond_first_minima_left:the upper limit of 1st neighbour bonds comming from the first minima 
                                    of bond length distribution, with sigma as lower limit.
            
        Warning:
            [x]What if there is no minima found? 
            count_1st_max / count_1st_minima > 10 used to check the effective of minima? 

        Examples:
            
        """
        #locate the 1st minima of bond-length distribution
        plt.figure()
        count_bins=plt.hist(self.voronoi.ridge_length[:],bins=20,range=[0,hist_cutoff])#
        _count=count_bins[0]
        _bins=count_bins[1]
        ridge_sorted=np.sort(self.voronoi.ridge_length[:])
        #find the minimum bin, then set the left side as bond_first_neighbour
        i_max=_bins.size
        i=0
        """1st half minimum count"""
        if method=="local_compare_method":
            while i < i_max-3:#i start from 0, which have to -1;compares i,i+1 and i+2, which have to -2, hence -3
                #np.where( self.bond_sorted[:]>res.bins[i] & self.bond_sorted[:]< res.bins[i+1] ) 
                #print(self._count[i])
                if _count[i] > _count[i+1]:
                    if _count[i+1] <= _count[i+2]:
                        i+=1
                        break
                i+=1 
        elif method=="median_minima_method":
            R"""
            n_ridges_half: it is believed that extremely short ridges are rare, 
            so the total number of short ridges should be less than half the num of ridges.
            hence n_ridges_half is set as a limit.
            """
            n_ridges_half = int(sum(_count)/2) 
            #i=1
            count_ridges_half = _count[0]
            count_min_i = 0
            while i < i_max-3:
                if count_ridges_half<n_ridges_half:
                    count_ridges_half = count_ridges_half + _count[i+1]
                else:
                    i=count_min_i#+1
                    break
                if _count[count_min_i] > _count[i+1]:
                    count_min_i = i+1
                i+=1

        elif method=="global_compare_method":
            _count = np.array(_count)
            count_sorted=np.sort(_count)
            """
            print(count_sorted)
            print(count_sorted[-1])#max1
            print(count_sorted[-2])#max2
            """
            i_bins_for_count_peak_1 = np.where(_count[:]==count_sorted[-1])
            i_bins_for_count_peak_2 = np.where(_count[:]==count_sorted[-2])
            i_bins_for_count_peak_1 = i_bins_for_count_peak_1[0]
            i_bins_for_count_peak_2 = i_bins_for_count_peak_2[0]
            #if there are more than one bins share the same count, select the smallest bin number.
            if np.shape(i_bins_for_count_peak_1)[0]>1:
                i_bins_for_count_peak_1 = min(i_bins_for_count_peak_1)
            if np.shape(i_bins_for_count_peak_2)[0]>1:
                i_bins_for_count_peak_2 = min(i_bins_for_count_peak_2)
            i_bins_for_count_1st_peak = min(i_bins_for_count_peak_1,i_bins_for_count_peak_2)
            i_bins_for_count_1st_peak = int(i_bins_for_count_1st_peak)#force i be a simple data_type(int)
            """
            print(i_bins_for_count_peak_1) 
            print(i_bins_for_count_peak_2)
            print(i_bins_for_count_1st_peak) 
            """
            #list_xy = np.logical_and(list_x,list_y)
            #edge_cut_positions_list = np.where(==)
            i = i_bins_for_count_1st_peak
            while i < i_max-3:#i start from 0, which have to -1;compares i,i+1 and i+2, which have to -2, hence -3
                if _count[i] > _count[i+1]:
                    if _count[i+1] <= _count[i+2]:
                        i+=1
                        break
                i+=1
        else:
            print("Err: input 'local_compare_method' or 'global_compare_method', please!")
            print("x")
        
        #check the length of ridge
        R"""
        short_ridge_rate: the variable is set to ensure that the ridge_length is not popular. 
        """
        n_ridges = sum(_count)
        short_ridge_rate = _count[i-1]/n_ridges
        if short_ridge_rate<0.1:#short ridge must be rare
            self.ridge_first_minima_left=_bins[i]
        else:
            self.ridge_first_minima_left=_bins[1]

        if not png_filename  is None:
            #plot divide line
            plt.plot([self.ridge_first_minima_left,self.ridge_first_minima_left],[0,_count[i]],linestyle='--')
            plt.title("1st minima:"+str( self.ridge_first_minima_left )+x_unit)
            plt.xlabel("ridge length "+x_unit)
            plt.ylabel("count (1)")
            plt.savefig(png_filename)
        plt.close()

    def get_conditional_bonds_and_simplices(self,long_bond_cutoff=6):
        R"""
        return:
            vertex_bonds_index: n rows of [start_point_index, end_point_index]
            list_simplex_cluster: n_vertices rows of [simplex_id, cluster_id],
                where delaunay.simplices[simplex_id], 
                and cluster_id is the cluster the simplex belonging to
        method:
            vertices cluster method: find the shortest ridges and link the related vertices.
        """

        list_short_ridge_bool = self.voronoi.ridge_length[:] <= self.ridge_first_minima_left
        list_short_bonds = self.voronoi.ridge_points[np.logical_not(list_short_ridge_bool)]#[start_point_index, end_point_index]
        """
        select_short_bonds_bool = self.bond_length[list_short_bonds[:,0],2] < long_bond_cutoff
        #select_short_bonds_bool = self.bond_length[list_short_bonds_index,2] < long_bond_cutoff
        list_short_bonds = list_short_bonds[select_short_bonds_bool]
        """
        list_ridge_vertices_index = self.voronoi.ridge_vertices[list_short_ridge_bool]
        list_vertices_index = np.unique(list_ridge_vertices_index)
        n_vertices = np.shape(list_vertices_index)[0]
        '''
        method2:
        record_vertices_cluster = np.ones((n_vertices,10))
        record_vertices_cluster[:]=-record_vertices_cluster[:]
        '''
        #record_vertices_cluster: n_vertices rows of [vertex_id, cluster_id],
        #where vertex_id belongs to list_vertices
        #when cluster_id=-1, the vertex has not been included by any cluster.
        record_vertices_cluster = np.ones((n_vertices,2))
        record_vertices_cluster[:]=-record_vertices_cluster[:]
        record_vertices_cluster[:,0]=list_vertices_index
        cluster_id = 0
        for i in range(n_vertices):
            vertex_id = list_vertices_index[i]
            #print("cluster_id\n",cluster_id)
            list_linked = np.where(list_ridge_vertices_index[:,:]==vertex_id)
            list_linked_vertices_pair = list_ridge_vertices_index[list_linked[0]]
            #print("vertices_pair\n",list_linked_vertices_pair)
            list_vertex_id_of_cluster1 = np.unique(list_linked_vertices_pair)
            for j in list_vertex_id_of_cluster1:
                record_id = record_vertices_cluster[:,0]==j#find the row_index for vertex_id
                if record_vertices_cluster[record_id,1]==-1:
                    record_vertices_cluster[record_id,1]=cluster_id
                elif record_vertices_cluster[record_id,1]==cluster_id:
                    pass
                else:
                    #new cluster merge old one by refreshing a list of vertices in old cluster
                    contradictory_cluster_id = record_vertices_cluster[record_id,1]
                    list_bool_of_cluster_to_merge= (record_vertices_cluster[:,1]==contradictory_cluster_id)
                    record_vertices_cluster[list_bool_of_cluster_to_merge,1]=cluster_id
            cluster_id +=1
        #statistics for clusters
        list_cluster_id = np.unique(record_vertices_cluster[:,1])
        count_polygon = np.zeros((20,2))#(10,2)
        for i in list_cluster_id:
            list_cluster_i = record_vertices_cluster[record_vertices_cluster[:,1]==i,0]
            n_complex_i = np.shape(list_cluster_i)[0]
            count_polygon[n_complex_i,0]=n_complex_i+2#[x]frame=1421,polygon=10;1609,12
            count_polygon[n_complex_i,1]+=1
        count_polygon[1,0]=3
        count_polygon[1,1]=np.shape(list_short_bonds)[0]
        #print("polygon_n, count\n",count_polygon)
        #
        count_polygon_relative = np.array(count_polygon)
        count_polygon_relative[:,1] = count_polygon_relative[:,1]/sum(count_polygon_relative[:,1]) \
                                        *(count_polygon_relative[:,0]-2)*100#see one simplex as one weight
        #print("polygon_n, count%\n",count_polygon_relative)
        """
        nm = np.maximum(count_polygon[:,1])
        np.shape(list_cluster_id)
        list_polygon_simplex = -np.ones((10,nm))
        for i in list_cluster_id:
            list_cluster_i = record_vertices_cluster[record_vertices_cluster[:,1]==i,0]
            list_polygon_simplex[]
        """
       
        self.vertex_bonds_index = list_short_bonds#short bonds index offered by ridge comapre method 
        self.list_simplex_cluster = record_vertices_cluster
        return count_polygon_relative

    def draw_bonds_simplex_conditional_oop(self,png_filename=None,xy_stable=None,nb_change=None,x_unit='(um)',
                                    LinearCompressionRatio=None,trap_filename=None,
                                    axis_limit=None,fig=None,ax=None):
        R"""
        xy(self.points) must be a frame of txyz( not txyz_stable)!
        xy_stable: a frame of txyz_stable
        LinearCompressionRatio=0.79
        trap_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12-part1"
        """
        if ax is None:
            self.draw_bonds = bond_plot_module()
        else:
            self.draw_bonds = bond_plot_module(fig,ax)    
        self.draw_bonds.restrict_axis_property_relative(self.points,x_unit=x_unit)
        if not axis_limit is None:
            xlim = [0,axis_limit[0]]
            ylim = [0,axis_limit[1]]
            self.draw_bonds.restrict_axis_limitation(xlim,ylim)
        self.draw_bonds.draw_points_with_conditional_vertices(self.points,self.bond_length,self.vertex_bonds_index)
        #self.draw_bonds.draw_bonds_conditional_bond(self.points,self.bond_length,bond_length_limmit=check,x_unit=x_unit)
        if not trap_filename is None:
            self.draw_bonds.plot_traps(trap_filename=trap_filename,LinearCompressionRatio=LinearCompressionRatio)
        if not nb_change is None:
            self.draw_bonds.plot_neighbor_change(xy_stable,nb_change)
        if not png_filename is None:
            self.draw_bonds.save_figure(png_filename)
        #else:
        #    return ax

    def draw_polygon_patch_oop(self,fig=None,ax=None):
        R"""
        parameters:
            vertex_bonds_index: n rows of [start_point_index, end_point_index]
            list_simplex_cluster: n_vertices rows of [simplex_id, cluster_id],
                where delaunay.simplices[simplex_id], 
                and cluster_id is the cluster the simplex belonging to
                (from get_conditional_bonds_and_simplices())
            is_polygon_n: mark the num n of edges of the polygon.
        
        """
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        if ax is None:
            fig,ax = plt.subplots()

        vertices_cluster = self.list_simplex_cluster
        points = self.points

        list_cluster_id = np.unique(vertices_cluster[:,1])
        #patches=[]
        for cluster_id in list_cluster_id:
            list_simplex_ids_in_cluster_i = vertices_cluster[vertices_cluster[:,1]==cluster_id,0]
            n_simplex_in_cluster_i = np.shape(list_simplex_ids_in_cluster_i)[0]
            cluster_i_is_polygon_n=n_simplex_in_cluster_i+2
            if cluster_i_is_polygon_n==6:
                #patches=[]
                for simplex_id in list_simplex_ids_in_cluster_i.astype(int):
                    list_points_ids_in_simplex = self.delaunay.simplices[simplex_id]
                    list_points_xy = points[list_points_ids_in_simplex]
                    """
                    polygon = Polygon(list_points_xy, closed=True)
                    polygon.set(color='r')
                    patches.append(polygon)
                    p = PatchCollection(patches)#, alpha=0.4
                    ax.add_collection(p)
                    """
                    set_color='r'
                    ax.fill(list_points_xy[:,0],list_points_xy[:,1],facecolor=set_color,edgecolor=set_color)
                    
        """
        #rearrange points anti-clockwise
        """
        
        return fig,ax

    def get_first_minima_radial_distribution_function(self,rdf,lattice_constant=3.0,png_filename=None):#[x]
        #print(rdf.bin_centers) print(rdf.bin_counts)
        self.lattice_constant=lattice_constant
        #locate the 1st minima of bond-length distribution
        self._count=rdf.bin_counts
        self._bins=rdf.bin_centers
        self.bond_sorted=np.sort(self.bond_length[:,2]/self.lattice_constant)
        #find the minimum bin, then set the left side as bond_first_neighbour
        
        "global_compare_method"
        self._count = np.array(self._count)
        self._bins = np.array(self._bins)
        self.count_sorted=np.sort(self._count)

        i_max=self._bins.size
        i=0
        """
        print(self.count_sorted)
        print(self.count_sorted[-1])#max1
        print(self.count_sorted[-2])#max2
        """
        #locate the 1st minima of bond-length distribution
        plt.figure()
        rdf.plot()
        """
        print(self.count_sorted)
        print(self.count_sorted[-1])#max1
        print(self.count_sorted[-2])#max2
        """
        i_bins_for_count_peak_1 = np.where(self._count[:]==self.count_sorted[-1])
        i_bins_for_count_peak_2 = np.where(self._count[:]==self.count_sorted[-2])
        i_bins_for_count_peak_1 = i_bins_for_count_peak_1[0]
        i_bins_for_count_peak_2 = i_bins_for_count_peak_2[0]
        #if there are more than one bins share the same count, select the smallest bin number.
        if np.shape(i_bins_for_count_peak_1)[0]>1:
            i_bins_for_count_peak_1 = min(i_bins_for_count_peak_1)
        if np.shape(i_bins_for_count_peak_2)[0]>1:
            i_bins_for_count_peak_2 = min(i_bins_for_count_peak_2)
        i_bins_for_count_1st_peak = min(i_bins_for_count_peak_1,i_bins_for_count_peak_2)
        i_bins_for_count_1st_peak = int(i_bins_for_count_1st_peak)#force i be a simple data_type(int)
        """
        print(i_bins_for_count_peak_1) 
        print(i_bins_for_count_peak_2)
        print(i_bins_for_count_1st_peak) 
        """
        #list_xy = np.logical_and(list_x,list_y)
        #self.edge_cut_positions_list = np.where(==)
        i = i_bins_for_count_1st_peak
        while i < i_max-3:#i start from 0, which have to -1;compares i,i+1 and i+2, which have to -2, hence -3
            if self._count[i] > self._count[i+1]:
                if self._count[i+1] <= self._count[i+2]:
                    if self._bins[i+1] > 1/self.lattice_constant:
                        # bond_length should be larger than sigma of particle
                        i+=1
                        break
            i+=1
            
        self.bond_first_minima_left=self._bins[i]
        self.bond_first_neighbour=self.bond_sorted[np.where(self.bond_sorted[:]<self.bond_first_minima_left)]
        
        if not png_filename  is None:
            #plot divide line
            plt.plot([self.bond_first_minima_left,self.bond_first_minima_left],[0,0.5],linestyle='--')#self._count[i]
            plt.title("1st peak:"+str( np.around( np.median(self.bond_first_neighbour),2) )+"+-"+str( np.around(np.std(self.bond_first_neighbour),2) ))
            plt.xlabel("interparticle distance (unit=A)")
            plt.ylabel("count (1)")
            #plt.show()
            #png_filename=prefix+"bond_hist_index"+str(index_num)
            plt.savefig(png_filename)
        plt.close()
        #let bond_first_minima_left be a absolute length, not normalized one
        self.bond_first_minima_left=self.bond_first_minima_left*self.lattice_constant
        
    def get_coordination_number_conditional(self,lattice_constant=3):
        #cut edge to remove CN012
        R"""
        Introduction:
            get_coordination_number_with given bond_length_limit [min,max].
        parameters:

        Variables:
            __coordination_bond: n rows of (start_point_index, end_point_index, bond_length)
            count_coordination: 10 rows of [count]. The i-th [count_i] represents the count 
                of points whose coordination number is i, where i from 0 to 9. 
            count_coordination_ratio: count_coordination whose sum is normalized to 1.
        """
        #check whether bond_first_minima_left exist
        try:
            bond_length_limit=self.bond_first_minima_left
        except AttributeError:
            self.get_first_minima_bond_length_distribution(lattice_constant=lattice_constant)
            bond_length_limit=self.bond_first_minima_left

        #select the bonds within the limit(bond length of 1st neighbour)
        self._bond_min=0.0
        self._bond_max=bond_length_limit
        place=np.where((self.bond_length[:,2]>self._bond_min) & (self.bond_length[:,2]<self._bond_max))
        self.__coordination_bond=self.bond_length[place,0:2]
        self.__coordination_bond=self.__coordination_bond[0]#remove a [] in [[]] structure
        
        
        #self._particle_id_max=np.max(self.__coordination_bond[:])#get the largest id of particle to restrain while loop
        #print(self.__coordination_bond)

        cn_max = 12#for honeycomb, cn10 truly appears in Index4340_seed9_frame2!
        self.count_coordination=np.zeros([cn_max,1])
        
        for id in self.edge_cut_positions_list[0]:
            place=np.where(self.__coordination_bond[:]==id)
            """
            Exception has occurred: IndexError       (note: full exception trace is shown but execution is paused at: <module>)
index 10 is out of bounds for axis 0 with size 10
  File "/home/tplab/xiaotian_file/lxt_code_py/points_analysis_2D.py", line 340, in get_coordination_number_conditional
    self.count_coordination[i]+=1
  File "/home/tplab/xiaotian_file/lxt_code_py/data_analysis_cycle.py", line 1024, in save_from_gsd
    a_frame.get_coordination_number_conditional()#cut edge to remove CN012
  File "/home/tplab/xiaotian_file/lxt_code_py/workflow_part.py", line 1715, in workflow_simu_to_mysql_pin_hex_to_kagome_cycle_oop_kT
    account=account)
  File "/home/tplab/xiaotian_file/lxt_code_py/control_table_multi_thread.py", line 13, in <module> (Current frame)
    index_end=tt.workflow_simu_to_mysql_pin_hex_to_kagome_cycle_oop_kT(index1=index1,lcr=lcr1,kT=kT,seed=seed)
            """
            i=len(place[0])
            #print(place[0])
            self.count_coordination[i]+=1
        #print(self.count_coordination)
        self.count_coordination_ratio=self.count_coordination/sum(self.count_coordination)
        #print(self.count_coordination_ratio)
        #return self.__coordination_bond

    def search_neighbor_id(self,id):
        #cut edge to remove CN012
        R"""
        Introduction:
            input a particle id and get the ids of whose neighbors
        Variables:
            __coordination_bond: n rows of (start_point_index, end_point_index)
            place: (2,Nneighbors)[neighbor_id,column]
        """
        place=np.where(self.__coordination_bond[:]==id)
        n_nb=len(place[0])#number of neighbors
        place_nb=np.array(place,dtype=np.int32) 
        #print(place_nb)
        place_nb[1]=1-place_nb[1]
        #print(place_nb)
        id_nb=self.__coordination_bond[place_nb[0],place_nb[1]].astype(np.int32)
        return id_nb

    def get_nearest_neighbor_dispalcement(self):
        R"""
        Introduction:
            lindemann_dispalcement: r(i,t)-<r(j,t)>
        Return: 
            list_id_dxy_nb: 
                list[particle_id,num_neighbors, sum_id_neighbors,dx,dy], 
                id of center particles and their displacements relative to neighbors.
                if the particle has new neighbors, unless inter-frame time is too long to
                catch the event of breaking old bond and forming new bond, number of neighbors
                will change together with the event.  
            dict_c_nb: 
                dict{'id_c':id_nb}, center particles and their neighbors.
                To monitor the event of breaking old bond and forming new bond and hence
                get the life span of a caging state.
        """
        
        #check whether bond_first_minima_left exist
        try:
            bond_length_limit=self.bond_first_minima_left
        except AttributeError:
            self.get_first_minima_bond_length_distribution()
            bond_length_limit=self.bond_first_minima_left

        #select the bonds within the limit(bond length of 1st neighbour)
        self._bond_min=0.0
        self._bond_max=bond_length_limit
        place=np.where((self.bond_length[:,2]>self._bond_min) & (self.bond_length[:,2]<self._bond_max))
        self.__coordination_bond=self.bond_length[place,0:2]
        self.__coordination_bond=self.__coordination_bond[0]#remove a [] in [[]] structure
        
        
        for id in self.edge_cut_positions_list[0]:
            place=np.where(self.__coordination_bond[:]==id)
            n_nb=len(place[0])#number of neighbors
            place_nb=np.array(place,dtype=np.int32) 
            #print(place_nb)
            place_nb[1]=1-place_nb[1]
            #print(place_nb)
            id_nb=self.__coordination_bond[place_nb[0],place_nb[1]].astype(np.int32)
            sum_nb=np.sum(id_nb)
            #print(id_nb)#id s 
            
            #print(self.points[id_nb])
            center_nb = np.average(self.points[id_nb],axis = 0)
            dxy_id = self.points[id]-center_nb
            #print(center_nb)
            if 'list_id_dxy_nb' in locals():
                list_id_dxy_nb.append([id,n_nb,sum_nb,dxy_id[0],dxy_id[1]]) 
                dict_c_nb[str(id)] = id_nb
            else:
                list_id_dxy_nb = []
                list_id_dxy_nb.append([id,n_nb,sum_nb,dxy_id[0],dxy_id[1]])
                dict_c_nb = {str(id):id_nb}
        
        return list_id_dxy_nb,dict_c_nb

    def get_neighbor_cloud_method_voronoi(self,png_filename=None):#checked right
        R"""
        Introduction:
            place all the (voronoi) neighbor particles into a plot 
            with the same position of central particles.
        Example:
            import points_analysis_2D
            prefix="/home/tplab/hoomd-examples_0/"
            filename = prefix +'testhex3-16-8'#index2572~cairo
            obj_of_simu_index = points_analysis_2D.PointsAnalysis2D(filename=filename)
            obj_of_simu_index.get_neighbor_cloud_method_voronoi()
        """
        (indptr, indices) = self.delaunay.vertex_neighbor_vertices

        plt.figure()
        for k in self.edge_cut_positions_list[0]:#range(num):
            nbr = indices[indptr[k]:indptr[k+1]]
            plt.scatter(self.points[nbr,0]-self.points[k,0],self.points[nbr,1]-self.points[k,1],c='k',marker='.') 

        plt.axis('equal')
        plt.xlabel('x(sigma)')
        plt.ylabel('y(sigma)')
        """
        the displayed image size will be the smaller one 
        between axis limitation for xlim/ylim or data itself. 
        Hence the effective xlim/ylim should be set smaller than data.

        To ensure that xlim/ylim < data region, 
        we add physically meaningless points.
        """
        dis=4*1.5#the distance from origin to xmax or ymax
        #to ensure that xlim/ylim < data region
        plt.scatter([dis,-dis],[dis,-dis],c='k',marker='.')
        #restrict ticks to show
        #(let all the images share the same size)
        new_ticks = np.linspace(-dis,dis,int(2*dis+1))
        new_ticks = new_ticks.astype(int)
        #print(new_ticks.astype(str))
        #print((new_ticks))
        plt.xticks(new_ticks,new_ticks.astype(str))
        plt.yticks(new_ticks,new_ticks.astype(str))
        #restrict data region to calculate
        plt.xlim(-dis,dis)
        plt.ylim(-dis,dis)

        if png_filename is None:
            plt.show()
        else:
            plt.savefig(png_filename)
            plt.close()

    def get_neighbor_cloud_method_1st_minima_bond(self,png_filename=None):#checked right
        R"""
        Introduction:
            place all the (1st bond minima) neighbor particles into a plot 
            with the same position of central particles.
        Example:
            import points_analysis_2D
            prefix="/home/tplab/hoomd-examples_0/"
            filename = prefix +'testhex3-16-8'#index2572~cairo
            obj_of_simu_index = points_analysis_2D.PointsAnalysis2D(filename=filename)
            obj_of_simu_index.get_neighbor_cloud_bond_minima()
        """
        #check whether bond_first_minima_left exist
        try:
            bond_length_limit=self.bond_first_minima_left
        except AttributeError:
            self.get_first_minima_bond_length_distribution()
            bond_length_limit=self.bond_first_minima_left
        
        #find pairs of particles whose distance < bond_length_limit
        import scipy
        kdt = scipy.spatial.KDTree(self.points)
        pairs = kdt.query_pairs(bond_length_limit)
        
        plt.figure()
        #count=np.zeros((self.points.shape[0],1))
        for i,j in pairs:
            if self.edge_cut_positions_bool[i]:
                #np.np.concatenate([array1,array2]) 
                # would be useful to accelerate [x]
                #edge should be cut [x]
                plt.scatter(self.points[j,0]-self.points[i,0],self.points[j,1]-self.points[i,1],c='k',marker='.')
                #count[i,0]=count[i,0]+1
                plt.scatter(self.points[i,0]-self.points[j,0],self.points[i,1]-self.points[j,1],c='k',marker='.')
                #count[j,0]=count[j,0]+1
        
        plt.axis('equal')
        plt.xlabel('x(sigma)')
        plt.ylabel('y(sigma)')
        """
        the displayed image size will be the smaller one 
        between axis limitation for xlim/ylim or data itself. 
        Hence the effective xlim/ylim should be set smaller than data.

        To ensure that xlim/ylim < data region, 
        we add physically meaningless points.
        """
        dis=4*1.5#the distance from origin to xmax or ymax
        #to ensure that xlim/ylim < data region
        plt.scatter([dis,-dis],[dis,-dis],c='k',marker='.')
        #restrict ticks to show
        #(let all the images share the same size)
        new_ticks = np.linspace(-dis,dis,int(2*dis+1))
        new_ticks = new_ticks.astype(int)
        #print(new_ticks.astype(str))
        #print((new_ticks))
        plt.xticks(new_ticks,new_ticks.astype(str))
        plt.yticks(new_ticks,new_ticks.astype(str))
        #restrict data region to calculate
        plt.xlim(-dis,dis)
        plt.ylim(-dis,dis)
        
        #print(count)

        if png_filename is None:
            plt.show()
        else:
            #example: "neighbor_cloud_1st_minima_index1369.png"
            plt.savefig(png_filename)
            plt.close()

    def draw_bonds_conditional_bond_cloud(self,check=[0.9,2.0],png_filename=None):
        R"""
        Parameters:
            png_filename: "prefix/bond_plot_index1513"
        Examples:
            import points_analysis_2D as pa
            s = "/home/tplab/Downloads/"
            index_num = 1387
            index_name = "index"+str(index_num)
            fname = s+index_name
            bb = pa.PointsAnalysis2D(filename=fname)
            oname1 = s+"bond_hist_"+index_name
            bb.draw_bond_length_distribution_and_first_minima(png_filename=oname1)
            oname2 = s+"bond_plot_"+index_name
            bb.draw_bonds_conditional_bond(check=[0.9, bb.bond_first_minima_left], png_filename=oname2)
        """
        
        #add lines for edges
        for i in range(np.shape(self.bond_length)[0]):
            if (self.bond_length[i,2] > check[0])&(self.bond_length[i,2] < check[1]) :
                edge = self.bond_length[i,0:2].astype(int)
                pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
                line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b')
                plt.gca().add_line(line)
        #plt.show()
        plt.title("bond_length:"+str(np.around(check,2))+"um")
        if not png_filename is None:
            plt.savefig(png_filename)
        #plt.close()

    def __cut_edge_of_positions(self,effective_lattice_constant = 3):
        R"""
        Variables:
            effective_lattice_constant： the distance between a particle and its neighbor
            edge_cut_positions_list: list the rows of particles' positions at given snapshot.
        """
        a = effective_lattice_constant
        #sz = len(self.init_positions)#np.size
        #xy = self.init_positions
        xmax = max(self.points[:,0])
        ymax = max(self.points[:,1])
        xmin = min(self.points[:,0]) 
        ymin = min(self.points[:,1])
        xmax = xmax - a
        ymax = ymax - a
        xmin = xmin + a
        ymin = ymin + a
        
        #That directly using 'and' has been banned, so 'logical_and' is necessary
        list_xmin = self.points[:,0] > xmin
        list_xmax = self.points[:,0] < xmax
        list_ymin = self.points[:,1] > ymin
        list_ymax = self.points[:,1] < ymax
        list_x = np.logical_and(list_xmin,list_xmax)
        list_y = np.logical_and(list_ymin,list_ymax)
        list_xy = np.logical_and(list_x,list_y)

        self.edge_cut_positions_list = np.where(list_xy)
        self.edge_cut_positions_bool = list_xy # T for body, F for edge.

    def get_bond_orientational_order(self,plot=False,png_filename=None,k_set=6):
        R"""
        Parameters:
            k_set: set the k-th disclination order parameter psi_k.
        """
        box=np.zeros(2)
        x1=min(self.points[:,0])
        x2=max(self.points[:,0])
        Lx=x2-x1#get box size in x-direction
        y1=min(self.points[:,1])
        y2=max(self.points[:,1])
        Ly=y2-y1#get box size in y-direction
        box[0]=Lx+1
        box[1]=Ly+1
        hex_order = freud.order.Hexatic(k=k_set)

        sp=np.shape(self.points)
        pts=np.zeros((sp[0],sp[1]+1))
        pts[:,0:2]=self.points[:]
        hex_order.compute(system=(box,pts[:]))

        self.Psi_k=hex_order.particle_order#local BOO
        self.Psi_k_abs=abs(hex_order.particle_order)
        self.Psi_k_global_with_edge=np.average(self.Psi_k_abs)#abs checked right
        self.Psi_k_global_cut_edge=np.average(self.Psi_k_abs[self.edge_cut_positions_list])#abs checked right

        if plot:
            plt.figure()
            plt.scatter(self.points[:,0],self.points[:,1],c=self.Psi_k_abs)
            plt.colorbar()
            if not png_filename is None:
                plt.savefig(png_filename)
            plt.close()

    def get_bond_orientational_order_selected(self,k_set=6,lower_limit=0.9):
        R"""
        Introduction:
            get psi_k ratio
        Parameters:
            k_set: set k-fold to calculate psi_k.
            lower_limit: those particles whose psi_k are larger than lower_limit 
                will be chosen as numerator.
        """
        #check whether Psi_k_abs exist
        try:
            pk=self.Psi_k_abs
        except AttributeError:
            self.get_bond_orientational_order(k_set=k_set)
            pk=self.Psi_k_abs

        place=np.where(pk[:]>lower_limit)
        sp_3=np.shape(pk[place])
        sp_all=np.shape(pk)
        self.Psi_k_rate=sp_3[0]/sp_all[0]

    def draw_bonds(self,fignum=1,show=False):
        #draw a figure with bonds
        plt.figure(fignum)
        plt.scatter(self.points[:,0],self.points[:,1],color='k')
        #add lines for edges
        for edge in self.bond_length[:,0:2].astype(int):
            #print(edge)
            pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
            #plt.gca().add_line(plt.Line2D(pt1,pt2))
            line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b')
            plt.gca().add_line(line)
        if show:
            plt.show()

    def draw_bonds_conditional_bond(self,check=[0.9,2.0],png_filename=None,nb_change=None,x_unit='(um)',
                                    show_traps=False,LinearCompressionRatio=0.79,
                                    trap_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12-part1"):
        R"""
        Parameters:
            png_filename: "prefix/bond_plot_index1513"
            nb_change: particle ids which change neighbors.
        Examples:
            import points_analysis_2D as pa
            s = "/home/tplab/Downloads/"
            index_num = 1387
            index_name = "index"+str(index_num)
            fname = s+index_name
            bb = pa.PointsAnalysis2D(filename=fname)
            oname1 = s+"bond_hist_"+index_name
            bb.draw_bond_length_distribution_and_first_minima(png_filename=oname1)
            oname2 = s+"bond_plot_"+index_name
            bb.draw_bonds_conditional_bond(check=[0.9, bb.bond_first_minima_left], png_filename=oname2)
        """
        bond_check= [check[0],check[1]]
        
        #if dis is None:
        xmax = max(self.points[:,0]) #- 3
        ymax = max(self.points[:,1]) #- 3
        xmin = min(self.points[:,0]) #+ 3
        ymin = min(self.points[:,1]) #+ 3
        dis = min(xmax-xmin,ymax-ymin)/2.0#half of the length of the system.
        dis = dis - 0 #cut the edge if necessary(eg. size & scale of images not match)

        center = [(xmax+xmin)*0.5,(ymax+ymin)*0.5]
        center_origin_distance = np.abs(np.dot(center,center))
        if  center_origin_distance < 1.0:# center is really close to origin
            center = [0,0]

        #draw a figure with edges
        plt.figure()
        plt.axis('equal')
        plt.xlabel('x'+x_unit)
        plt.ylabel('y'+x_unit)
        plt.title("bond_length:"+str(np.around(check,2))+x_unit)
        """
        the displayed image size will be the smaller one 
        between axis limitation for xlim/ylim or data itself. 
        Hence the effective xlim/ylim should be set smaller than data.

        To ensure that xlim/ylim < data region, 
        we add physically meaningless points.
        """
        #restrict ticks to show
        """
        #(let all the images share the same size)
        new_ticks = np.linspace(-dis,dis,int(2*dis+1))
        new_ticks = new_ticks.astype(int)
        #print(new_ticks.astype(str))
        #print((new_ticks))
        plt.xticks(new_ticks,new_ticks.astype(str))
        plt.yticks(new_ticks,new_ticks.astype(str))
        """
        #restrict data region to show
        plt.xlim(-dis+center[0],dis+center[0])
        plt.ylim(-dis+center[1],dis+center[1])
        
        #add lines for edges
        for i in range(np.shape(self.bond_length)[0]):
            if (self.bond_length[i,2] > bond_check[0])&(self.bond_length[i,2] < bond_check[1]) :
                edge = self.bond_length[i,0:2].astype(int)
                pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
                line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b')
                plt.gca().add_line(line)

        plt.scatter(self.points[:,0],self.points[:,1],color='k')
        """
        particle_ids = np.linspace(0,self.points.shape[0]-1,self.points.shape[0],dtype=int)
        particle_ids_str = particle_ids.astype(str)
        for j in particle_ids:
            plt.annotate(particle_ids_str[j],self.points[j])
        """
        if not nb_change is None:
            plt.scatter(self.points[nb_change,0],self.points[nb_change,1],color='orange')
        
        if show_traps:
            traps=np.loadtxt(trap_filename)
            traps=np.multiply(traps,LinearCompressionRatio)
            plt.scatter(traps[:,0], 
                   traps[:,1],
                   c='r',marker = 'x')

        if not png_filename is None:
            plt.savefig(png_filename)
        plt.close()
    
    def draw_bonds_conditional_bond_oop(self,check,png_filename=None,xy_stable=None,nb_change=None,x_unit='(um)',
                                    LinearCompressionRatio=None,trap_filename=None,
                                    axis_limit=None):
        R"""
        xy(self.points) must be a frame of txyz( not txyz_stable)!
        xy_stable: a frame of txyz_stable
        LinearCompressionRatio=0.79
        trap_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12-part1"
        """
        self.draw_bonds = bond_plot_module()
        self.draw_bonds.restrict_axis_property_relative(self.points,x_unit=x_unit)
        if not axis_limit is None:
            xlim = [0,axis_limit[0]]
            ylim = [0,axis_limit[1]]
            self.draw_bonds.restrict_axis_limitation(xlim,ylim)
        self.draw_bonds.draw_points_with_conditional_bond(self.points,self.bond_length,bond_length_limmit=check)
        #self.draw_bonds.draw_bonds_conditional_bond(self.points,self.bond_length,bond_length_limmit=check,x_unit=x_unit)
        if not trap_filename is None:
            self.draw_bonds.plot_traps(trap_filename=trap_filename,LinearCompressionRatio=LinearCompressionRatio)
        if not nb_change is None:
            self.draw_bonds.plot_neighbor_change(xy_stable,nb_change)
        if not png_filename is None:
            self.draw_bonds.save_figure(png_filename)

    def draw_bonds_conditional_bond_for_image_oop(self,image,check,png_filename=None,xy_stable=None,nb_change=None,x_unit='(um)',
                                    LinearCompressionRatio=None,trap_filename=None,
                                    axis_limit=None):
        R"""
        xy(self.points) must be a frame of txyz( not txyz_stable)!
        xy_stable: a frame of txyz_stable
        LinearCompressionRatio=0.79
        trap_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12-part1"
        axis_limit=[xmin,xmax,ymin,ymax]
        """
        self.draw_bonds = bond_plot_module_for_image(image)
        self.draw_bonds.restrict_axis_property_relative(self.points,x_unit=x_unit)
        if not axis_limit is None:
            xlim = [axis_limit[0],axis_limit[1]]
            ylim = [axis_limit[2],axis_limit[3]]
            self.draw_bonds.restrict_axis_limitation(xlim,ylim)
        line = self.draw_bonds.draw_points_with_conditional_bond(self.points,self.bond_length,bond_length_limmit=check)
        #self.draw_bonds.draw_bonds_conditional_bond(self.points,self.bond_length,bond_length_limmit=check,x_unit=x_unit)
        if not trap_filename is None:
            self.draw_bonds.plot_traps(trap_filename=trap_filename,LinearCompressionRatio=LinearCompressionRatio)
        if not nb_change is None:
            self.draw_bonds.plot_neighbor_change(xy_stable,nb_change)
        if not png_filename is None:
            self.draw_bonds.save_figure(png_filename)
        return line
    
    def print_benchmark(self,bond=True,ridge=True):
        if bond :
            print('bond')
            print('min',min(self.bond_length[:,2]))
            print('median',np.median(self.bond_length[:,2]))
            print('max',max(self.bond_length[:,2]))
            print('\n')
        if ridge :
            print('ridge')
            print('min',min(self.voronoi.ridge_length))
            print('median',np.median(self.voronoi.ridge_length))
            print('max',max(self.voronoi.ridge_length))
            print('\n')

    def print_basic_property(self):
        R"""
        but this class can only get points and hence doesnt know energy parameters.
        [x]
        """
        #interaction
        epsilon=300
        kappa=0.25
        #substrate
        k=800.0
        rcut=1.0

        energy_interaction=epsilon*np.exp(-kappa)# when rcut=12.5, Uparticle = kBT = 1; rcut=15, Up=0.47
        energy_substrate=0.5*k*np.multiply(rcut,rcut)
        #pt=np.exp(1)
        print('energy_thermal     1')
        print('energy_interaction '+str(energy_interaction.astype(int)))
        print('energy_substrate   '+str(energy_substrate.astype(int)))
    
    def get_cluster_vertices(self):
        #check whether cutted_ridges exist 
        try:
            v_index_temp=self.voronoi.ridge_vertices[(self.voronoi.cutted_ridges).astype(int)]
        except AttributeError:
            self.draw_bonds_conditional_ridge()
            v_index_temp=self.voronoi.ridge_vertices[(self.voronoi.cutted_ridges).astype(int)]
            
        v_index_temp=v_index_temp[:,0]#v_index_temp is [[array]]
        #get index of vertices
        c=np.concatenate([v_index_temp[:,0],v_index_temp[:,1]])
        c=np.sort(c)
        shp=np.shape(c)
        for j in range(shp[0]-1):
            if c[j]==c[j+1]:
                c[j+1]=-1
        c=np.sort(c)
        for j in range(shp[0]-1):
            if not c[j]==-1:
                start=j
                break
        c=c[j:shp[0]]
        self.voronoi.cluster_vertices=c
        #plt.scatter(self.voronoi.vertices[c,0],self.voronoi.vertices[c,1],color='g')
        #plt.show()
        del c,v_index_temp

    def count_n_polygon_fast(self,n_polygon=6):
        R"""
        arguments:
            n_polygon: annotate n-polygon
        """
        #check whether cluster_vertices exist
        try:
            cv_temp=self.voronoi.vertices[self.voronoi.cluster_vertices]
        except AttributeError:
            self.get_cluster_vertices()
            cv_temp=self.voronoi.vertices[self.voronoi.cluster_vertices]
        
        shp=np.shape(cv_temp)
        cluster_count=np.zeros((shp[0]))
        check=self.voronoi.cut_ridge_length
        cv_dis=distance.cdist(cv_temp,cv_temp,'euclidean')
        
        for i in range(shp[0]):
            '''
            if cv_dis[i,0]==255:
                continue
            else:
                list=np.where(cv_dis[i,:]<check)
                print(list[0])
                #cv_dis[:,list]=255#delete found vertices 
                cv_dis[list,:]=255#delete found vertices 
                count=np.shape(list[0])
                print(cv_dis[0:5,0:5])
                count=count[0]#how many vertices in a cluster
                cluster_count[count]=cluster_count[count]+1
            '''
            list=np.where(cv_dis[i,:]<1.0*check)
            cv_dis[:,list]=255#delete found vertices 
            cv_dis[list,:]=255#delete found vertices 
            count=np.shape(list[0])
            count=count[0]#how many vertices in a cluster
            cluster_count[count]=cluster_count[count]+1
            if count==(n_polygon-2):
                plt.annotate(cluster_count[count].astype(int), (cv_temp[i,0], cv_temp[i,1]))
        """
        #some great polygons whose vertices linked closely(small dis),
        # while twisted polygons whose vertices are farther from each other 
        # than from different vertices clusters.
        # so an algorithm with more precision is needed 
        """
        print('Tetragon,Pentagon,Hexagon...')
        print(cluster_count[2:10])
        
    def count_polygon_n(self,n_edges):   
        R"""
        Introduction:
            a polygon with n edges must be merged by n-2 triangles, hence n-3 bonds must be removed. 
            That is, a quadrilateral contains 2 triangles with 1 bond removed.
        parameters:
            points: n rows of (x,y)
            --------------------------------
            delaunay: scipy.spatial.Delaunay
            bond_length: n rows of (start_point_row, end_point_row, bond_length)
            cutted_bonds: n rows of (start_point_row, end_point_row)
            simplices: n rows of (point_row_1, point_row_2, point_row_3) vertices of triangles
            --------------------------------
            linked_triangles: 
        """
        #check whether cutted_bonds exist
        if not ('self.cutted_bonds' in dir()):#if hasattr(self,'box')
            self.draw_bonds_conditional_bond_oop()
        #get linked triangles pair   
        shp=np.shape(self.cutted_bonds)
        self.linked_triangles = np.zeros(shp)
        for i in range(shp[0]):
            #1st round search
            index1=np.where(self.delaunay.simplices[:]==self.cutted_bonds[i,0])
            #2nd round search
            index2=np.where(self.delaunay.simplices[index1[0]]==self.cutted_bonds[i,1])
            temp=index1[0]
            self.linked_triangles[i]=temp[index2[0]]
        
        #get linked triangles chain
        chain=np.zeros((shp[0],8))#up to 10-polygon
        chain[:]=-1
        all_chain_count=np.zeros((shp[0],1))
        chain_num=1#now we are proceeding the chain_num-th chain [x]unfinished
        tri_count_max=shp[0]#at most the triangle pairs should be linked
        lt_temp=self.linked_triangles[:]#np.zeros((shp[0],2))

        i=0#if all the triangle pairs are read, break the loop
        while sum(all_chain_count)<tri_count_max:
            if lt_temp[i,0]==-1:
                #this triangle pair has been linked
                i=i+1
                continue
            else:
                #init a new chain
                this_chain_count=1
                lt_temp[i,:]=-1#delete the head triangle pair to avoid being searched
                chain[chain_num,this_chain_count-1]=self.linked_triangles[i,0]
                chain[chain_num,this_chain_count]=self.linked_triangles[i,1]
                #search linked triangles
                index1=np.where((lt_temp[:,0]==self.linked_triangles[i,1])|(lt_temp[:,1]==self.linked_triangles[i,1]))#rightward
                index_check1=np.shape(index1[0])
                if index_check1[0]==0:
                    print('error,the triangle pair rightward not found') 
                elif index_check1[0]==1:
                    #find triangles
                    if index1[1]==0:
                        self.linked_triangles[index1[0],1]
                        lt_temp[index1[0],:]=-1
                        this_chain_count=this_chain_count+1

                    elif index1[1]==1:#
                        self.linked_triangles[index1[0],0]
                        lt_temp[index1[0],:]=-1
                        this_chain_count=this_chain_count+1

                    else:
                        print('error')
                            
                elif index_check1[0]==2:
                    #two braches
                    for j in index1[0]:
                        if lt_temp[j,0]== self.linked_triangles[i,1]:
                            pass 

                else:
                    print('error')
                    
                #leftward
            if i==tri_count_max:
                #if all the triangle pairs are read, break the loop
                break
    
        '''
        for i in range(shp[0]):
            if i==0:
                chain[0,0:2]=lt_temp[i]
                lt_temp[i]=[-1,-1]#clear the i-th triangle to avoid being found
                chain_count[count]=chain_count[count]+1#the count of chain-0 is added one
                count=count+1#now there is one more chain
            else:
                index1=np.where((chain[:]==lt_temp[i,0])|(chain[:]==lt_temp[i,1]))
                temp=index1[0]#row-index of linked triangles
                for j in temp:
                    if chain_count[j,0]:
                    lt_temp[j,1]
                #that chain[0] and chain[5] should be linked may exist

        del lt_temp,temp
        
        #https://blog.csdn.net/u013378642/article/details/81775131
        self.linked_triangles_sorted=sorted(self.linked_triangles,key=lambda x:x[0])
        print('t1\n',self.linked_triangles,'\nt2\n',self.linked_triangles_sorted)
        #a chain of triangle could be like [1,162] [162,4] [4,88.]
        '''
        pass

    def count_polygon_n_new(self,n_edges):
        R"""
        Introduction:
            a polygon with n edges must be merged by n-2 triangles, hence n-3 bonds must be removed. 
            That is, a quadrilateral contains 2 triangles with 1 bond removed.
        parameters:
            points: n rows of (x,y)
            --------------------------------
            delaunay: scipy.spatial.Delaunay
            bond_length: n rows of (start_point_row, end_point_row, bond_length)
            cutted_bonds: n rows of (start_point_row, end_point_row)
            simplices: n rows of (point_row_1, point_row_2, point_row_3) vertices of triangles
            --------------------------------
            linked_triangles:
            list_polygon_n: 8
        """
        list_polygon_n = range(10)

        return list_polygon_n
    
    def plot_polygon(self,x,y):
         fig,ax = plt.subplots()
         ax.fill(x,y,color='g')

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

    def plot_trajectory(self,length_cut_edge=0):#checked right
        R"""
        Example:
            import getDataAndScatter as scatt
            simu_index = 1369
            prefix_gsd = '/home/tplab/hoomd-examples_0/trajectory_auto'
            postfix_gsd = '.gsd'
            filename_gsd = prefix_gsd+str(simu_index)+postfix_gsd
            scatt.save_image_stack(gsd_file=filename_gsd)
        """
        #from celluloid import camera
        matplotlib.use(backend="agg")#Backend agg is non-interactive backend. Turning interactive mode off.
        #Backend Qt5Agg is interactive backend. Turning interactive mode on.
        self.trajectory#self.trajectory=gsd.hoomd.open(self.filename_gsd)#open a gsd file
        self.num_of_frames#self.num_of_frames=len(self.trajectory)

        frame_num=0
        dis = None
        while frame_num < self.num_of_frames:
            snap=self.trajectory.read_frame(frame_num)#take a snapshot of the N-th frame
            positions=snap.particles.position

            plt.figure()
            if dis is None:
                xmax = max(positions[:,0]) #- 3
                ymax = max(positions[:,1]) #- 3
                xmin = min(positions[:,0]) #+ 3
                ymin = min(positions[:,1]) #+ 3
                dis = min(xmax,ymax,-xmin,-ymin)
                dis = dis - length_cut_edge #cut the edge if necessary(eg. size & scale of images not match)
            plt.scatter(positions[:,0],positions[:,1])

            plt.axis('equal')
            plt.xlabel('x(sigma)')
            plt.ylabel('y(sigma)')
            """
            the displayed image size will be the smaller one 
            between axis limitation for xlim/ylim or data itself. 
            Hence the effective xlim/ylim should be set smaller than data.

            To ensure that xlim/ylim < data region, 
            we add physically meaningless points.
            """
            #restrict ticks to show
            """
            #(let all the images share the same size)
            new_ticks = np.linspace(-dis,dis,int(2*dis+1))
            new_ticks = new_ticks.astype(int)
            #print(new_ticks.astype(str))
            #print((new_ticks))
            plt.xticks(new_ticks,new_ticks.astype(str))
            plt.yticks(new_ticks,new_ticks.astype(str))
            """
            #restrict data region to show
            plt.xlim(-dis,dis)
            plt.ylim(-dis,dis)
            
            #plt.show()

            prefix_old="/home/tplab/hoomd-examples_0"
            folder_name=self.filename_gsd.strip(prefix_old)
            folder_name=folder_name.strip(".gsd")#folder_name=prefix+png_filename_as_folder
            folder_name="t"+folder_name#"t" is necessary in case of "/t" deleted before
            prefix='/home/tplab/Downloads/'
            png_filename=prefix+folder_name+"/"+folder_name+"_"+str(frame_num)
            #check if the folder exists
            isExists=os.path.exists(prefix+folder_name)
            if isExists:
                plt.savefig(png_filename)
            else:
                os.makedirs(prefix+folder_name)
                plt.savefig(png_filename)
            plt.close()

            frame_num=frame_num+1

        return folder_name
    
    def get_displacement_field_xy(self,frame_index,plot=False,png_filename=None):
        R"""
        Introduction:
            The function draws a displacement vector field from init state to final state 
            with positions at edge removed  to clean abnormal displacement vector. 
        Example:
            import points_analysis_2D as pa
            gsd = pa.proceed_gsd_file(simu_index=1382)
            gsd.get_displacement_field(plot=True)
        """
        self.gdf_xy = dynamic_points_analysis_2d(self.txyz_stable)
        self.gdf_xy.displacement_field_module()
        self.gdf_xy.displacemnt_field.get_displacement_field_xy(frame_index,plot,png_filename)

    def get_displacement_field_distribution(self,frame_index,log_mode=False,png_filename=None):
        R"""
        Introduction:
            The function draws a displacement vector field from init state to final state 
            with positions at edge removed  to clean abnormal displacement vector. 
        Example:
            import points_analysis_2D as pa
            gsd = pa.proceed_gsd_file(simu_index=1382)
            gsd.get_displacement_field(plot=True)
        """
        self.gdf_xy.displacemnt_field.get_displacement_field_distribution(frame_index,log_mode,png_filename)

    def get_gr(self,frame_num):
        rdf = freud.density.RDF(bins=200, r_max=20.0,)#

        rdf.compute(system=self.trajectory.read_frame(frame_num))
        """
        print(rdf.bin_centers)
        print(rdf.bin_counts)
        rdf.plot()
        data_filename=prefix+'index_'+str_index+'gr.png'
        plt.savefig(data_filename)
        """
        plt.close()

class proceed_exp_file:
    R"""
    see particle_tracking.py to get trajectories of particles from a video.
    """
    pass

class dynamic_points_analysis_2d:#old_class_name: msd
    R"""
    Introduction:
        this class is designed to analyze dynamic properties of a group of paritcle trajectories.
    
    Methods:
        msd_module: mean square displacements.
        displacement_field_module: displacement field.
    
    example:
        #final cut bond plot of a simulation
        import data_analysis_cycle as dac
        import points_analysis_2D as pa
        import numpy as np
        daw = dac.data_analysis_workflow()
        directory,dataname=daw.gsd_to_txyz(simu_index=5387)
        txyz = np.load('/home/remote/Downloads/5387_9/txyz.npy')
        dpa = pa.dynamic_points_analysis_2d(txyz)
        dpa.plot_bond_neighbor_change_oop(data_name=dataname,prefix=directory,final_cut=True,
                        trap_filename='/home/remote/hoomd-examples_0/testhoneycomb3-8-12-part1',
                        trap_lcr=0.816)
    """
    def __init__(self,txyz_stable,mode='simu'):
        R"""
        parameters:
            txyz_stable:array[Nframes,Nparticles,xyz],
                for simu data, 
                    ensure that particles never move across the boundary(box)!
                for exp data, 
                    unit of xyz must be um! 
                    ensure that particles are always in the vision field!
            #BOX: [lx,ly,lz,xy,xz,yz]
            account: 'remote' or 'tplab'
            mode: 'simu' or 'exp'. 'simu' to select particles in box; 'exp' to direct compute msd
        """            
        if mode == 'simu':
            #self.box = box
            #self.__select_stable_trajectories_simu(plot_trajectory)
            self.txyz_stable = txyz_stable
            self.x_unit = '($\sigma$)'
            self.t_unit = '(step)'
        elif mode == 'exp':
            self.txyz_stable = txyz_stable
            self.x_unit = '(um)'
            self.t_unit = '(s)'

        #self.account = account
        self.mode = mode

    def plot_trajectory_all_in_one(self,png_prefix=''):
        frames,particles,dimensions=self.txyz_stable.shape
        list_stable_id = range(particles)#txyz_stable.shape[1]
        plt.figure()
        for index_particle in list_stable_id:
            txyz_ith = self.txyz_stable[:,index_particle,:]
            plt.plot(txyz_ith[:,0],txyz_ith[:,1])
        plt.xlabel("$x$ "+self.x_unit )
        plt.ylabel("$y$ "+self.x_unit )
        #png_filename = '/home/'+self.account+'/Downloads/'+'traj_stable.png'
        png_filename=png_prefix+'trajectory_stable.png'
        plt.savefig(png_filename)
        plt.close() 

    def plot_trajectory_single_particle_loop(self,png_prefix=''):
        R"""
        EXP:
            path_to_results = '/home/remote/xiaotian_file/data/20221129/video_5'
            txyz_npy_filename = path_to_results+'/'+'txyz_stable.npy'
            txyz_stable = np.load(txyz_npy_filename)
            msds = pa.dynamic_points_analysis_2d(txyz_stable,mode='exp')
            msds.plot_trajectory_single_particle(path_to_results+'/')
        """
        frames,particles,dimensions=self.txyz_stable.shape
        list_stable_id = range(particles)#txyz_stable.shape[1]
        for particle_id in list_stable_id:
            txyz_ith = self.txyz_stable[:,particle_id,:]
            plt.figure()
            plt.plot(txyz_ith[:,0],txyz_ith[:,1])
            plt.xlabel("$x$ "+self.x_unit )#'(sigma)'
            plt.ylabel("$y$ "+self.x_unit )
            png_filename = 'traj_stable_'+str(int(particle_id))+'.png'
            plt.savefig(png_prefix+png_filename)
            plt.close()   

    def plot_trajectory_single_particle(self,particle_id,png_prefix=''):
        R"""
        """
        txyz_ith = self.txyz_stable[:,particle_id,:]
        plt.figure()
        plt.plot(txyz_ith[:,0],txyz_ith[:,1])
        plt.xlabel("$x$ "+self.x_unit )#'(sigma)'
        plt.ylabel("$y$ "+self.x_unit )
        png_filename = 'traj_stable_'+str(int(particle_id))+'.png'
        plt.savefig(png_prefix+png_filename)
        plt.close()   

    def plot_a_frame_of_points(self,frame_index,png_filename):
        points = self.txyz_stable[frame_index]
        fig,ax = plt.subplots()
        ax.scatter(points[:,0],points[:,1],c='k')
        plt.savefig(png_filename)

    def msd_module(self):
        R"""
        return:
            (class)self.msd
        """
        self.msd = mean_square_displacement(self.txyz_stable)
    
    def displacement_field_module(self):
        R"""
        return:
            (class)self.displacemnt_field
        """
        self.displacemnt_field = displacemnt_field_2D(self.txyz_stable)
    
    def compute_nearest_neighbor_displacements(self,csv_prefix='',bond_cut_off=6):
        R"""
        parameters:
            input: 
                trajectory_stable (unit = sigma or um)
                'stable' means all the particles in the 
                field of view from the 1st frame to the end.
                if the trajectory is not stable, a particle id 
                in different frames would not represent the same particle.

                csv_prefix: (string)the directory for csv file to save.

            return: 
                ts_id_dxy:
                    a pandas.dataframe contains ["particle_id","num_neighbors","sum_id_neighbors", "x","y", "frame"]
                    which contains all the nn displacements over time.
                    Caution: hist_cutoff should be set nearly 2 times of bond length, 
                    or the 1st minima bond length would be incorrectly recognized.
                dict_c_nb:
                    a dict() contains {'particle id': [list id of neighbors]} at given frame.
                    
        introduction:
            DOI: 10.1103/PhysRevE.95.022602
            
            nearest-neighbor relative mean square displacement(nnmsd): <[r(i,t)-<r(j,t)>]^2>
            where j is the neighbors of particle i.

            The dynamic Lindemann parameter: gamma_l_tau = nnmsd/2/a^2
        question:
            if I should calculate the particle averaged first neighbor bond length( A(t) ), 
            and calculate the averaged time bond length(<A>)? 
            dynamic lindemann parameter suits only for the time-scale when the system is stable solid.
        """
        id_dxy_pd_columns = ["particle_id","num_neighbors","sum_id_neighbors", "x","y"]
        #get frame-wise 1st neighbor bond length
        #list_1st_bond_length_framewise
        self.frames,self.particles,self.dimensions=self.txyz_stable.shape
        list_frames = range(self.frames)
        list_1st_bond_length_framewise = np.zeros((self.frames,))
        for frame in list_frames:
            result = static_points_analysis_2d(self.txyz_stable[frame,:,:])
            result.get_first_minima_bond_length_distribution(lattice_constant=1,hist_cutoff=bond_cut_off)#here lattice_constant is just used to normalize figure, hence set 2.0 is ok
            #print('recognized bond length: '+str(result.bond_length_median*lc)+'+-'+str(result.bond_length_std*lc)+' '+unit)
            list_1st_bond_length_framewise[frame] = result.bond_length_median#unit = sigma!
            #displacement neighbor
            id_dxy,dict_c_nb = result.get_nearest_neighbor_dispalcement()
            id_dxy_pd = pd.DataFrame(id_dxy)
            """
            list_id_dxy_nb: 
                list[particle_id,num_neighbors, sum_id_neighbors,dx,dy], 
                id of center particles and their displacements relative to neighbors.
                if the particle has new neighbors, unless inter-frame time is too long to
                catch the event of breaking old bond and forming new bond, number of neighbors
                will change together with the event.  
            """
            #print(id_dxy_pd.head())
            if frame == 0:
                #ts_id_dxy = t_id_dxy
                id_dxy_pd.columns = id_dxy_pd_columns
                id_dxy_pd['frame'] = frame
                ts_id_dxy = id_dxy_pd
                #print(ts_id_dxy.head())
                #print(ts_id_dxy.tail())
            else:
                #ts_id_dxy = np.concatenate((ts_id_dxy,t_id_dxy), axis=0) #[frame,id,dx,dy]
                id_dxy_pd.columns = id_dxy_pd_columns#[id,dx,dy,frame]
                id_dxy_pd['frame'] = frame
                ts_id_dxy = pd.concat([ts_id_dxy,id_dxy_pd])
                #print(ts_id_dxy.tail())
                #print(ts_id_dxy.shape)
        #select stable trajectories
        """
        print(ts_id_dxy[:,1].max)
        list_ids = np.unique(ts_id_dxy[:,1])
        id_check = np.zeros((list_ids.shape[0],),dtype=bool)
        for id in list_ids:
            ts_dxy = ts_id_dxy[ts_id_dxy[:,1]==id]

            if ts_dxy.shape[0] == self.frames:
                id_check[id] = True
        ts_id_dxy_stable = ts_id_dxy[:,id_check,:,:]
        """
        average_1st_bond_length = np.average(list_1st_bond_length_framewise) 
        #particle_tracking.select_stable_trajectory()is great!
        pd.DataFrame.to_csv(ts_id_dxy,csv_prefix+'ts_id_dxy.csv')
        return ts_id_dxy,average_1st_bond_length
    
    def monitor_neighbor_change_event(self,ts_id_dxy,csv_prefix=''):
        R"""
        input:
            ts_id_dxy: 
                a pandas.dataframe contains ["particle_id","num_neighbors","sum_id_neighbors", "x","y", "frame"]
                which contains all the nn displacements over time. 
                get from self.compute_nearest_neighbor_displacements()
            csv_prefix:
                (string)the directory for csv file to save.
        return:
            list_sum_id_nb_stable:
                (csv file)["frame", "particle_id", "sum_id_neighbors", 'if_nb_change'].
            if_nb_change_int:
                (ndarray)[Nframe_nb_stable,Nparticle_nb_stable]. 0 nb stable, 1 nb change.
            n_particle_nb_stable:(int)
        consideration:
            given edge cut in each frame, some of the particles are removed,
            check each particle_id if contains Nframes.
            count id = Nframes, then record; or delete
        """
        #scan particle_id
        #list_sum_id_nb_stable_np [Nframes, Nparticles]: sum_id_neighbors (int)
        list_sum_id_nb = ts_id_dxy[["frame", "particle_id", "sum_id_neighbors"]]
        """
        file_list_sum_id_nb = prefix + 'list_sum_id_nb.csv'
        pd.DataFrame.to_csv(list_sum_id_nb,file_list_sum_id_nb)
        """
        list_id = list_sum_id_nb['particle_id'].values
        list_id_counts = np.unique(list_id,return_counts=True) 
        list_id_counts = np.array(list_id_counts)

        list_id_nb_not_cut_edge= list_id_counts[0,list_id_counts[1,:] == max(list_id_counts[1])]#list_id_counts[1].max]
        n_frame_nb_stable = max(list_id_counts[1])
        n_particle_nb_stable = len(list_id_nb_not_cut_edge)
        
        list_sum_id_nb_stable = list_sum_id_nb[list_sum_id_nb["particle_id"].isin(list_id_nb_not_cut_edge)]
        list_sum_id_nb_stable_np = list_sum_id_nb_stable["sum_id_neighbors"].values
        list_sum_id_nb_stable_np = np.reshape(list_sum_id_nb_stable_np,(n_frame_nb_stable,n_particle_nb_stable))
        #sum_id [Nframes-1, Nparticles]: sum_id_neighbors (int)
        #sum_id_0 [0 ~ Nframes-1, Nparticles] where Nparticles are never cut by edges.
        #sum_id_1 [1 ~ Nframes, Nparticles]
        sum_id_0 = list_sum_id_nb_stable_np[:-1]
        sum_id_1 = list_sum_id_nb_stable_np[1:]
        #sum_id_diff = sum_id_1 - sum_id_0
        sum_id_diff = np.absolute(sum_id_1 - sum_id_0)
        # if sum_id_diff[:,:] != 0 then = 1
        #compare sum_id_neighbors between frame+1 and frame
        #different -> change_neighbor = true
        sum_id_diff = np.array(sum_id_diff,dtype=bool)
        if_nb_change = np.array(list_sum_id_nb_stable_np,dtype=bool)
        if_nb_change[0] = False
        if_nb_change[1:] = sum_id_diff
        if_nb_change_int = np.array(if_nb_change,dtype=int)
        #reshape if_nb_change and plug back into list_sum_id_nb_stable
        list_sum_id_nb_stable.loc[:,'if_nb_change'] = tuple(np.reshape(if_nb_change,n_frame_nb_stable*n_particle_nb_stable)) 
        file_list_sum_id_nb_stable = csv_prefix + 'list_sum_id_nb_stable.csv'
        pd.DataFrame.to_csv(list_sum_id_nb_stable,file_list_sum_id_nb_stable)
        return if_nb_change_int,n_particle_nb_stable
        
    def get_hist_neighbor_change_event(self,if_nb_change_int,n_particle_nb_stable,prefix=None):
        # histogram: frame VS count change_neighbor_events
        count_nb_change_event = np.sum(if_nb_change_int,axis=1)
        count_nb_change_event_rate = count_nb_change_event/n_particle_nb_stable
        np.save(prefix+'count_nb_change_event_rate',count_nb_change_event_rate)
        self.plot_hist_neighbor_change_event(count_nb_change_event_rate,prefix)
    
    def plot_hist_neighbor_change_event(self,count_nb_change_event_rate,prefix=None):
        # histogram: frame VS count change_neighbor_events
        fig,ax = plt.subplots()
        #ax.plot(count_nb_change_event_rate)
        #ax.semilogx(count_nb_change_event_rate)
        ax.loglog(count_nb_change_event_rate)
        ax.set_xlabel('time steps(1000)')
        ax.set_ylabel('nb_change_particles(%)')
        if not prefix is None:
            png_filename = 'hist_nb_change_loglog.png'#_semilog
            fig.savefig(prefix+png_filename)

    def plot_bond_neighbor_change(self,data_name='default_exp',prefix='',
                                    final_cut=False,nb_change=None,bond_cut_off=6,
                                    show_traps=False,trap_filename=None,trap_lcr=None):
        R"""
        Introduction:
            final_cut: true to proceed the last frame only.
            bond_cut_off: should always be set the same condition in a data_set!
        """
        #prefix='/home/'+account+'/Downloads/'#'/home/tplab/Downloads/'
        #log_prefix='/home/'+account+'/hoomd-examples_0/'#'/home/tplab/hoomd-examples_0/'
        #load time steps
        str_index=data_name
        num_of_frames = self.txyz_stable.shape[0]
        
        for i in range(num_of_frames):
            if final_cut:
                i = num_of_frames-1#i=9#!!! 23
            
            a_frame = static_points_analysis_2d(points=self.txyz_stable[i])#hide_figure=False

            if not nb_change is None:
                #num_frames = list_sum_id_nb_stable['frame'].values.max()+1
                #frames = range(num_frames)
                #for frame in frames:
                list_sum_id_nb_stable = nb_change
                snap = list_sum_id_nb_stable[list_sum_id_nb_stable['frame']==i]
                snap_part = snap[snap['if_nb_change'] == True]
                ids = snap_part["particle_id"].values
                #points_nb_chg = txyz_stable[frame,ids,:2]

            if final_cut:
                #bond_plot+trap_plot
                png_filename1 = prefix +'bond_hist_index'+str_index+'_'+str(int(i))+'.png'
                png_filename2 = prefix +'bond_plot_1st_minima_index'+str_index+'_'+str(int(i))+'.png'
            else:
                folder_name=prefix+"record_"+str_index#+"/"
                #check if the folder exists
                isExists=os.path.exists(folder_name)
                if not isExists:
                    os.makedirs(folder_name)
                #bond_plot+trap_plot
                png_filename1 = folder_name+"/" +'bond_hist_index'+str_index+'_'+str(int(i))+'.png'
                png_filename2 = folder_name+"/" +'bond_plot_1st_minima_index'+str_index+'_'+str(int(i))+'.png'
            
            a_frame.get_first_minima_bond_length_distribution(lattice_constant=1,hist_cutoff=bond_cut_off,png_filename=png_filename1)#,png_filename=png_filename1
            a_frame.draw_bonds_conditional_bond(check=[0.4, a_frame.bond_first_minima_left], png_filename=png_filename2,
                                            show_traps=show_traps,LinearCompressionRatio=trap_lcr,trap_filename=trap_filename,
                                            nb_change=ids,x_unit=self.x_unit)
        
            if final_cut:
                break
    
    def plot_bond_neighbor_change_oop(self,data_name='default_exp',prefix='',
                                    final_cut=False,init_cut=False,nb_change=None,bond_cut_off=6,
                                    trap_filename=None,trap_lcr=None):#show_traps=False,
        R"""
        Introduction:
            final_cut: true to proceed the last frame only.
            bond_cut_off: should always be set the same condition in a data_set!
        return:
            
        """
        #prefix='/home/'+account+'/Downloads/'#'/home/tplab/Downloads/'
        #log_prefix='/home/'+account+'/hoomd-examples_0/'#'/home/tplab/hoomd-examples_0/'
        #load time steps
        str_index=data_name
        num_of_frames = self.txyz_stable.shape[0]
        if self.mode == 'simu':
            txyz = np.load(prefix+'txyz.npy')
        elif self.mode == 'exp':
            txyz = self.txyz_stable
        for i in range(num_of_frames):
            if final_cut:
                i = num_of_frames-1#i=9#!!! 23
            if init_cut:
                i = 0
            a_frame = static_points_analysis_2d(points=txyz[i])#hide_figure=False

            if nb_change is None:
                ids = None
            else:
                #num_frames = list_sum_id_nb_stable['frame'].values.max()+1
                #frames = range(num_frames)
                #for frame in frames:
                list_sum_id_nb_stable = nb_change
                snap = list_sum_id_nb_stable[list_sum_id_nb_stable['frame']==i]
                snap_part = snap[snap['if_nb_change'] == True]
                ids = snap_part["particle_id"].values
                #points_nb_chg = txyz_stable[frame,ids,:2]

            if final_cut or init_cut:
                #bond_plot+trap_plot
                png_filename1 = prefix +'bond_hist_index'+str_index+'_'+str(int(i))+'.png'
                png_filename2 = prefix +'bond_plot_1st_minima_index'+str_index+'_'+str(int(i))+'.png'
                ids = None
            else:
                folder_name=prefix+"record_"+str_index#+"/"
                #check if the folder exists
                isExists=os.path.exists(folder_name)
                if not isExists:
                    os.makedirs(folder_name)
                #bond_plot+trap_plot
                png_filename1 = folder_name+"/" +'bond_hist_index'+str_index+'_'+str(int(i))+'.png'
                png_filename2 = folder_name+"/" +'bond_plot_1st_minima_index'+str_index+'_'+str(int(i))+'.png'
            
            a_frame.get_first_minima_bond_length_distribution(lattice_constant=1,hist_cutoff=bond_cut_off,png_filename=png_filename1)#,png_filename=png_filename1
            #a_frame.draw_bonds_conditional_bond(check=[0.4, a_frame.bond_first_minima_left], png_filename=png_filename2,
            #                               show_traps=show_traps,LinearCompressionRatio=trap_lcr,trap_filename=trap_filename,
            #                                nb_change=ids,x_unit=self.x_unit)
            a_frame.draw_bonds_conditional_bond_oop(check=[0.4, a_frame.bond_first_minima_left], png_filename=png_filename2,
                                                    xy_stable=self.txyz_stable[i],nb_change=ids,x_unit=self.x_unit,
                            LinearCompressionRatio=trap_lcr, trap_filename=trap_filename)
            if final_cut or init_cut:
                break

    def plot_neighbor_change_evolution(self,frame_init,frame_final,prefix='',nb_change=None,arrow=False,
                                    data_name='default_exp',ids=None,
                                    bond_cut_off=None,
                                    trap_filename=None,trap_lcr=None):#show_traps=False,):
        R"""
        input:
            arrow: 'direct' or 'annotate'
            bond_cut_off: 6 recommended(2 times lattice constant)
        return:
            ids
        plots:    
            black dots: particles of initial state.
            orange circles: particles of final state.
            green arrows: displacements between final nad initial states.
        """
        str_index=data_name
        txyz = np.load(prefix+'txyz.npy')

        xy_init = self.txyz_stable[frame_init,:,:2]
        xy_final = self.txyz_stable[frame_final,:,:2]
        particle_size=50
        circle_size = particle_size
        lw = 2#circle linewidths
        circle_color = 'limegreen'#'orange'
        arrow_color = 'limegreen'

        bpm = bond_plot_module()
        bpm.restrict_axis_property_relative(xy_init,'($\sigma$)')
        #bpm.ax.set_xlim(3,16)#plt.xlim(-dis+center[0],dis+center[0])
        #bpm.ax.set_ylim(-12,1)#plt.ylim(-dis+center[1],dis+center[1])

        #draw bonds
        if not bond_cut_off is None:
            a_frame = static_points_analysis_2d(points=xy_init)
            a_frame.get_first_minima_bond_length_distribution(lattice_constant=1,hist_cutoff=bond_cut_off)#,png_filename=png_filename1
            check=[0.4, a_frame.bond_first_minima_left]
            bpm.draw_points_with_conditional_bond(xy_init,a_frame.bond_length,check,particle_size=particle_size)
        else:
            bpm.draw_points_with_conditional_bond(xy_init)
        #bpm.ax.set_xlim(3,16)#plt.xlim(-dis+center[0],dis+center[0])
        #bpm.ax.set_ylim(-12,1)#plt.ylim(-dis+center[1],dis+center[1])

        #draw circles & arrows of the next frame
        if ids is None:
            if nb_change is None:
                ids = None
                bpm.ax.scatter(xy_final[:,0],xy_final[:,1],facecolors='none',edgecolors=circle_color,marker='o',s=circle_size,linewidths=lw)
            else:
                list_sum_id_nb_stable = nb_change
                snap = list_sum_id_nb_stable[list_sum_id_nb_stable['frame']==frame_final]
                snap_part = snap[snap['if_nb_change'] == True]
                ids = snap_part["particle_id"].values
                
                if arrow == 'direct':#direct_arrow
                    bpm.ax.quiver(xy_init[ids,0],xy_init[ids,1],xy_final[ids,0]-xy_init[ids,0],xy_final[ids,1]-xy_init[ids,1],color=arrow_color,angles='xy', scale_units='xy', scale=1)
                elif arrow == 'annotate':#annotate_arrow
                    dxy,uv,ids = self.__annotate_arrow(xy_init,xy_final,ids)
                    # move arrow a little away from in situ points.
                    bpm.ax.quiver(xy_init[ids,0]+dxy[:,0],xy_init[ids,1]+dxy[:,1],uv[:,0],uv[:,1],color=arrow_color,angles='xy', scale_units='xy', scale=1)

                bpm.ax.scatter(xy_final[ids,0],xy_final[ids,1],facecolors='none',edgecolors=circle_color,marker='o',s=circle_size,linewidths=lw)#circle, facecolors='none',fillstyle='none'[x]
                """
                ax.annotate("", xy=(0.5, 0.5), xytext=(0, 0),
                ...             arrowprops=dict(arrowstyle="->"))
                https://matplotlib.org/stable/tutorials/text/annotations.html#sphx-glr-tutorials-text-annotations-py
                """
        else:# force ax to draw scatter using given ids
            bpm.ax.scatter(xy_final[ids,0],xy_final[ids,1],facecolors='none',edgecolors=circle_color,marker='o',s=circle_size,linewidths=lw)

        #draw traps
        if not trap_filename is None:
            bpm.plot_traps(trap_filename,trap_lcr,'map')
            

        png_filename2 = prefix +'string_like_motion'+str_index+'_'+str(int(frame_init))+'.png'#'.pdf'
        bpm.save_figure(png_filename2)

        return ids

    def plot_string_like_motion(self,frame_init,frame_final,prefix='',
                                    data_name='default_exp',ids=None,
                                    bond_cut_off=None,
                                    trap_filename=None,trap_lcr=None):
        R"""
        return 
            ids: particles whose displacemnt is large
        """
        str_index=data_name
        txyz = np.load(prefix+'txyz.npy')

        xy_init = self.txyz_stable[frame_init,:,:2]
        xy_final = self.txyz_stable[frame_final,:,:2]
        particle_size=30
        circle_size = particle_size
        lw = 2#circle linewidths
        circle_color = 'limegreen'#'orange'
        arrow_color = 'limegreen'

        
        bpm = bond_plot_module()
        bpm.restrict_axis_property_relative(xy_init,'($\sigma$)')  
        #draw bonds
        a_frame = static_points_analysis_2d(points=xy_init)
        a_frame.get_first_minima_bond_length_distribution(lattice_constant=1,hist_cutoff=bond_cut_off)#,png_filename=png_filename1
        check=[0.4, a_frame.bond_first_minima_left]
        bpm.draw_points_with_conditional_bond(xy_init,a_frame.bond_length,check,particle_size=particle_size)  
        #draw arrows
        df2 = displacemnt_field_2D(self.txyz_stable,bpm.ax,bpm.fig)  
        ids_in = df2.get_displacement_field_xy_id(frame_init,frame_final)
        #draw circles
        if ids is None:
            ids = ids_in
            bpm.ax.scatter(xy_final[ids,0],xy_final[ids,1],facecolors='none',
                edgecolors=circle_color,marker='o',s=circle_size,linewidths=lw,zorder =2)
        else:
            bpm.ax.scatter(xy_final[ids,0],xy_final[ids,1],facecolors='none',
                edgecolors=circle_color,marker='o',s=circle_size,linewidths=lw,zorder =2)#circle, facecolors='none',fillstyle='none'[x]
        #draw traps
        if not trap_filename is None:
            bpm.plot_traps(trap_filename,trap_lcr,'map')
            
        #save figure
        png_filename2 = prefix +'string_like_motion'+str_index+'_'+str(int(frame_init))+'_'+str(int(frame_final))+'.png'#'.pdf'
        bpm.save_figure(png_filename2)     
        return ids_in    

    def plot_string_like_motion_rank(self,frame_init,frame_final,prefix='',
                                    data_name='default_exp',nb_change=None,ids=None,
                                    bond_cut_off=None,
                                    trap_filename=None,trap_lcr=None):
        R"""
        return :
            ids: particles whose displacemnt is large
        warning:
            ids would lost partly for edge_cut ins static_points_analysis_2D
        """
        str_index=data_name
        txyz = np.load(prefix+'txyz.npy')

        xy_init = self.txyz_stable[frame_init,:,:2]
        xy_final = self.txyz_stable[frame_final,:,:2]
        particle_size=30
        circle_size = particle_size
        lw = 2#circle linewidths
        circle_color = 'limegreen'#'orange'
        arrow_color = 'limegreen'
        
        
        bpm = bond_plot_module()
        bpm.restrict_axis_property_relative(xy_init,'($\sigma$)')  
        #draw bonds
        a_frame = static_points_analysis_2d(points=xy_init)
        a_frame.get_first_minima_bond_length_distribution(lattice_constant=1,hist_cutoff=bond_cut_off)#,png_filename=png_filename1
        check=[0.4, a_frame.bond_first_minima_left]
        bpm.draw_points_with_conditional_bond(xy_init,a_frame.bond_length,check,particle_size=particle_size)  
        #draw arrows
        df2 = displacemnt_field_2D(self.txyz_stable,bpm.ax,bpm.fig)  
        df2.get_displacements(frame_final,frame_init)
        uv,ids_active = df2.select_long_displacement_arrow()
        #df2.get_displacement_field_xy_id(frame_init)

        rank = df2.compare_displacement_method(uv,ids_active,frame_init,frame_final)
        #not so good, nan
        #rank = self.neighbor_change_id_method(nb_change,frame_init,frame_final,ids_active)
        
        df2.get_displacement_field_xy_rank(rank,frame_init,frame_final)

        #draw circles
        if ids is None:
            ids = ids_active
            bpm.ax.scatter(xy_final[ids,0],xy_final[ids,1],c=rank,
                marker='o',s=circle_size,linewidths=lw,zorder =2)#edgecolors=,
        else:
            bpm.ax.scatter(xy_final[ids,0],xy_final[ids,1],facecolors='none',
                edgecolors=circle_color,marker='o',s=circle_size,linewidths=lw,zorder =2)#circle, facecolors='none',fillstyle='none'[x]
        #draw traps
        if not trap_filename is None:
            bpm.plot_traps(trap_filename,trap_lcr,'map')
            
        #save figure
        png_filename2 = prefix +'string_like_motion_rank'+str_index+'_'+str(int(frame_init))+'_'+str(int(frame_final))+'.png'#'.pdf'
        bpm.save_figure(png_filename2)     
        return ids_active      
    
    def neighbor_change_id_method(self,nb_change,frame_init,frame_final,ids_active):
        list_sum_id_nb_stable = nb_change
        list_frames = (list_sum_id_nb_stable['frame']>=frame_init)&(list_sum_id_nb_stable['frame']<=frame_final)
        snap = list_sum_id_nb_stable[list_frames]
        snap_true = snap[snap['if_nb_change'] == True]
        #print(snap_true.head(10))
        
        ids_active = np.array(ids_active)
        num_ids_active = np.sum(ids_active.astype(int))
        rank = np.zeros((num_ids_active,))
        list_ids_active = np.zeros((num_ids_active,))
        #get list_ids_active
        j=0
        for i in range(len(ids_active)):
            if ids_active[i]:
                list_ids_active[j]=i
                j=j+1
                if j==num_ids_active:
                    break
        #print(list_ids_active)
        j=0
        for i in list_ids_active.astype(int):
            id_trajectory = snap_true[snap_true["particle_id"] == i]
            rank[j] = id_trajectory["frame"].max()
            j=j+1
        print(rank) #some ids are nan for edge_cut
        #or xy is in need to check displacement.
        return rank

    def __annotate_arrow(self,xy_init,xy_final,ids):
        R"""
        input:
            xy_init: from txyz_stable
            xy_final: from txyz_stable
            ids: defined by neighbor change events.
        output:
            dxy: vectors from points to places to annotate arrows, selected by ids only.
            uv: vectors of arrow, selected by ids only.
            ids: selected by displacements(>0.2)
        """
        #normalize arrow length part1
        uv = xy_final[ids]-xy_init[ids]
        dr2 = uv*uv
        dr = np.sqrt(dr2[:,0]+dr2[:,1])
        
        #select particles with long displacement
        dr_long_list = dr[:]>0.2
        ids=ids[dr_long_list]
        uv = xy_final[ids]-xy_init[ids]
        dr2 = uv*uv
        dr = np.sqrt(dr2[:,0]+dr2[:,1])

        ##normalize arrow length part2
        uv[:,0] = uv[:,0]/dr
        uv[:,1] = uv[:,1]/dr
        rotate = np.zeros((2,2))
        rotate[0,0] = np.cos(np.pi/2)
        rotate[1,1] = np.cos(np.pi/2)
        rotate[1,0] = np.sin(-np.pi/2)
        rotate[0,1] = np.sin(np.pi/2)
        # move arrow a little away from in situ points.
        scale_away = 0.5
        #rotate operator is on the right side, 
        # so rotational direction is inverted.
        #to let rotation match intuition, 
        # I inverted rotation operator too.
        dxy = np.matmul(uv,rotate)*scale_away
        #tune the length of arrow
        scale_arrow =1
        uv = uv*scale_arrow
        
        return dxy,uv,ids


class mean_square_displacement:
    R"""
    Introduction:
        defined by wikipedia(https://en.wikipedia.org/wiki/Mean_squared_displacement)
        this function is designed to calculate mean squared displacements.
        particle index: i
        displacement: dr(i,t-0) = r(i,t) - r(i,0)
        squared displacement: [ dr(i,t-0) ]^2
        mean squared displacement: MSD(t) = < [ dr(i,t-0) ]^2 >_i
        averaged time MSD: MSD(dt) where dt = t2 -t1 to replace t -0, which is used to improve statistical performance.
        caution: averaged time MSD would hide dynamic events like collapse.

    example_experiment:
        import numpy
        import particle_tracking as ptt
        ana = ptt.particle_track()
        dir_path='/home/remote/xiaotian_file/data/20221129/DefaultVideo_5'
        traj = ana.select_stable_trajectory()
        traj_um = traj/32*3
        traj_sigma = traj_um/2
        pa = msd(traj_sigma,mode='exp')
        pa.compute_atmsd_t_chips(interval_max=0.9)
        path = '/home/remote/xiaotian_file/data/20221129/video5/'
        time_log = np.loadtxt(path+'DefaultVideo_5.txt')

        pa.plot_msd(time_log,path+'msd_chips_loglog_'+'20221129_video5_90%'+'.png')
        pa.plot_trajectory(path+'trajectory_stable_'+'20221129_video5'+'.png')

    EXAMPLE:
        print(msds.result_msd)
        print(np.shape(msds.result_msd))
    example:
        import points_analysis_2D as pa
        gsd_data = pa.proceed_gsd_file(simu_index=5208,seed=9)
        gsd_data.get_trajectory_data()

        msd_class = pa.msd(gsd_data.txyz,gsd_data.box)

        import freud
        msds = freud.msd.MSD(gsd_data.box)#the class is fault,,'direct'
        msds.compute(positions=msd_class.txyz_stable)
        import matplotlib.pyplot as plt 
        plt.figure()
        plt.plot(msds.msd)
        plt.title("Mean Squared Displacement")
        plt.xlabel("$t$")
        plt.ylabel("MSD$(t)$")
        png_filename = 'msd_'+'index5208_9'+'.png'
        plt.savefig(png_filename)#png_filename
        plt.close()
        #print(tr.shape)
    exp-get nn msd:
        ana = ptt.particle_track()
        traj = ana.select_stable_trajectory()
        traj_um = traj/32*3
        traj_sigma = traj_um/2
        pa = dynamic_points_analysis_2d(traj_sigma,mode='exp')
        ts_id_dxy = pa.compute_nearest_neighbor_displacements(unit='sigma')
        import pandas as pd
        ts_id_dxy = pd.read_csv('ts_id_dxy.csv')
        ts_id_dxy['z'] = 0
        pa = dynamic_points_analysis_2d(ts_id_dxy,mode='exp')
        txyz_ids_stable = pa.compute_nearest_neighbor_displacements_stable(pa.txyz_stable)
        pa = dynamic_points_analysis_2d(txyz_ids_stable,mode='exp')
        pa.compute_atmsd_t_chips(0.9)
        path_to_folder = '/home/remote/xiaotian_file/data/20221129/video5/'
        time_log = path_to_folder+'DefaultVideo_5.txt'
        time_log = np.loadtxt(time_log)
        pa.plot_msd_t_chips(time_log)
    """
    def __init__(self,txyz_stable,mode):
        self.txyz_stable = txyz_stable
        self.frames,self.particles,self.dimensions=self.txyz_stable.shape
        self.mode = mode

    def compute_msd_normal(self):
        m_max = int(0.5*self.frames)#to ensure the robustness of statistics
        #k_max = 0
        start = 1
        self.record_msd=np.zeros((m_max-start+1,2))#[interval m, msd_m]
        
        list_m=np.arange(start,m_max+1)
        self.record_msd[:,0]=list_m
        for m in list_m:
            self.dm_xyz = self.txyz_stable[m:,:,:] - self.txyz_stable[0,:,:]
            sq_dm_xyz2 = np.square(self.dm_xyz)
            sum_sq_dm_r2 = np.sum(sq_dm_xyz2)
            self.record_msd[m-start,1]=sum_sq_dm_r2/self.dm_xyz.shape[0]/self.particles

    def compute_atmsd_scan_t(self):
        R"""
        method:
            scan t axis, 1 frame for each interval. here may exist some overlaps between two intervals who share the same length of time interval m.
        record_msd:[interval m, msd_m]
        k: the frame to start scan 
        m: time interval
        m_max: max time interval
        """
        m_max = 9000#int(0.5*self.frames)#to ensure the robustness of statistics
        #k_max = 0
        start = 1000
        self.record_msd=np.zeros((m_max-start+1,2))#[interval m, msd_m]
        
        list_m=np.arange(start,m_max+1)
        self.record_msd[:,0]=list_m
        for m in list_m:
            self.dm_xyz = self.txyz_stable[m:,:,:] - self.txyz_stable[:-m,:,:]
            sq_dm_xyz2 = np.square(self.dm_xyz)
            sum_sq_dm_r2 = np.sum(sq_dm_xyz2)
            self.record_msd[m-start,1]=sum_sq_dm_r2/self.dm_xyz.shape[0]/self.particles
    
    def compute_atmsd_scan_t_log_record(self):
        R"""
        method:
            scan t axis, 1 frame for each interval. here may exist some overlaps between two intervals who share the same length of time interval m.
        record_msd:[interval m, msd_m]
        k: the frame to start scan 
        m: time interval
        m_max: max time interval
        """
        frames = np.log10(self.frames)
        print(self.frames)
        bins = 20
        alpha = np.linspace(0,frames,bins)
        list_m = np.int32(np.power(10,alpha)) 
        print(list_m)
        self.record_msd=np.zeros((np.shape(list_m)[0],2))#[interval m, msd_m]
        
        snap_to_step = 1000
        self.record_msd[:,0]=list_m*snap_to_step#+500000-100
        i=0
        for m in list_m:
            self.dm_xyz = self.txyz_stable[m:,:,:] - self.txyz_stable[:-m,:,:]
            sq_dm_xyz2 = np.square(self.dm_xyz)
            sum_sq_dm_r2 = np.sum(sq_dm_xyz2)
            self.record_msd[i,1]=sum_sq_dm_r2/self.dm_xyz.shape[0]/self.particles
            i=i+1
        
        return self.record_msd

    def compute_atmsd_t_chips(self,interval_max=0.1,msd_per_particle=False):
        R"""
        parameters:
            interval_max: (0,1) to which the longest interval msd is computed.
            k: the frame to start scan 
            m: time interval
            m_max: max time interval
            msd_per_particle: control the program to calculate msd for each particle
        
        return:
            self.record_msd: [interval m, atmsd(m)]
            or record_msd_id

        method:
            cut t axis into several chips for each interval. there are not any two 
            intervals who share the same length of time interval m would overlap.
            Hence we get the averaged time mean square displacement(atmsd).

        """
        if not hasattr(self,'frames'):
            self.frames,self.particles,self.dimensions=self.txyz_stable.shape
        m_max = int(interval_max*self.frames)#to ensure the robustness of statistics
        start = 1
        self.record_msd=np.zeros((m_max-start+1,2))#[interval m, msd_m]
        #you can select points [log 1, 1/100*log m_max,log m_max]
        list_m=np.arange(start,m_max+1)
        self.record_msd[:,0]=list_m
        
        for m in list_m:
            chips = int((self.txyz_stable.shape[0]-1)/m)# how many chips the time axis is divided into
            list_frames = m*np.arange(0,chips+1)
            m_xyz = self.txyz_stable[list_frames,:,:]
            self.dm_xyz = m_xyz[1:,:,:] - m_xyz[:-1,:,:]
            sq_dm_xyz2 = np.square(self.dm_xyz)
            sum_sq_dm_r2 = np.sum(sq_dm_xyz2)
            self.record_msd[m-start,1]=sum_sq_dm_r2/self.dm_xyz.shape[0]/self.particles

            if msd_per_particle:
                id_sum_sq_dm_r2_chips = np.sum(sq_dm_xyz2,axis=2)
                print(id_sum_sq_dm_r2_chips.shape)
                id_sum_sq_dm_r2 = np.sum(id_sum_sq_dm_r2_chips,axis=0)
                print(id_sum_sq_dm_r2.shape)
                if 'record_msd_id' in locals():
                    record_msd_id[m-start,:]=id_sum_sq_dm_r2/self.dm_xyz.shape[0]
                else:
                    record_msd_id=np.zeros((m_max-start+1,self.particles))
                    record_msd_id[m-start,:]=id_sum_sq_dm_r2/self.dm_xyz.shape[0]

        if msd_per_particle:
            return record_msd_id

    def compute_nn_msd_X(self,unit='sigma'):
        R"""
        parameters:
            input: trajectory (unit = sigma or um)

        introduction:
            DOI: 10.1103/PhysRevE.95.022602
            
            nearest-neighbor relative mean square displacement(nnmsd): <[r(i,t)-<r(j,t)>]^2>
            where j is the neighbors of particle i.

            The dynamic Lindemann parameter: gamma_l_tau = nnmsd/2/a^2
        """
        #get frame-wise 1st neighbor bond length
        #list_1st_bond_length_framewise
        list_frames = range(self.frames)
        list_1st_bond_length_framewise = np.zeros((self.frames,))
        for frame in list_frames:
            result = static_points_analysis_2d(self.txyz_stable[frame,:,:])
            if unit == 'sigma':
                lc = 1.0#lattice constant, particle diameter(2 um) as default
            elif unit == 'um':
                lc = 2.0
            result.get_first_minima_bond_length_distribution(lattice_constant=lc,hist_cutoff=5)#here lattice_constant is just used to normalize figure, hence set 2.0 is ok
            #print('recognized bond length: '+str(result.bond_length_median*lc)+'+-'+str(result.bond_length_std*lc)+' '+unit)
            list_1st_bond_length_framewise[frame] = result.bond_length_median#unit = sigma!
            #displacement neighbor
            #result.bond_first_neighbour
            id_dxy = result.get_nearest_neighbor_dispalcement()
            t_id_dxy = np.zeros((id_dxy.shape[0],id_dxy.shape[1]+1))  
            t_id_dxy[:,0] = frame
            t_id_dxy[:,1:] = id_dxy  
            if frame == 0:
                ts_id_dxy = t_id_dxy
            else:
                ts_id_dxy = np.concatenate((ts_id_dxy,t_id_dxy), axis=0) #[frame,id,dx,dy]
            #print(ts_id_dxy.shape)
        #select stable trajectories
        print(ts_id_dxy[:,1].max)
        list_ids = np.unique(ts_id_dxy[:,1])
        id_check = np.zeros((list_ids.shape[0],),dtype=bool)
        for id in list_ids:
            ts_dxy = ts_id_dxy[ts_id_dxy[:,1]==id]

            if ts_dxy.shape[0] == self.frames:
                id_check[id] = True
        ts_id_dxy_stable = ts_id_dxy[:,id_check,:,:]
        #particle_tracking.select_stable_trajectory()is great!

    def plot_msd_t_chips(self,time_log=None,png_filename='msd_chips_long_loglog.png',sigma=False,um_sec=True,lindemann=False):
        R"""
        introduction:
            input: 'txyz_stable.csv' 
            output: msd plot
        """
        import matplotlib.pyplot as plt 
        plt.figure()
        if time_log is None:
            plt.semilogx(self.record_msd[:,0],self.record_msd[:,1])#loglog semilogx
        else:
            time_msd=time_log[:self.record_msd.shape[0]]
            plt.loglog(time_msd,self.record_msd[:,1])

        if self.mode =='simu':
            plt.title("Mean Squared Displacement")
            plt.xlabel("$t$ (steps)" )
            plt.ylabel("MSD$(t)(\sigma^2)$ ")
        elif self.mode == 'exp':
            if lindemann:
                plt.title("Mean Squared Displacement")
                plt.xlabel("$t$ (sec)")
                plt.ylabel("MSD$(t)$ (sigma^2/A^2)")
            elif sigma:
                plt.title("Mean Squared Displacement")
                plt.xlabel("$t$ (sec)")
                plt.ylabel("MSD$(t)$ (sigma^2)")
            elif um_sec:
                plt.title("Mean Squared Displacement")
                plt.xlabel("$t$ (sec)")
                plt.ylabel("MSD$(t)$ (um^2)")
        #png_filename = 'msd_chips_long_loglog_'+'index5208_9'+'.png'
        plt.savefig(png_filename)#png_filename
        plt.close()
        """
        list_dm_frames = np.arange(self.dm_xyz.shape[0])
        for i in list_dm_frames:
            #np.dot()
            self.dm_xyz_tr=np.transpose(self.dm_xyz)
            np.matmul(self.dm_xyz),self.dm_xyz_tr)
        """
        
        
        """
        m=0#frame interval
        record_disp2_1 = np.zeros((self.sp[0]-1,3))
        while m<self.sp[0]-1:#m < N frame
            k=0#frame scanning
            record_disp2_m = np.zeros((1,3))
            
            while k<self.sp[0]-1-m:
                disp_m = positions[k+m]-positions[k]#199-0 what is the problem?
                disp2_m = np.dot(np.transpose(disp_m),disp_m) 

                if disp2_m[0,0]>0:
                    print(m)

                record_disp2_m = record_disp2_m + [disp2_m[0,0],disp2_m[1,1],disp2_m[2,2]]
                k=k+1

            record_disp2_m = record_disp2_m/(self.sp[0]-1-m)#mean of all the chips with k frames interval
            record_disp2_1[m] = record_disp2_m[0]
            
            m=m+1
        self.result_msd_xyz = record_disp2_1/self.sp[1]#divided by N_particles
        self.result_msd = np.zeros((self.sp[0]-1,1))
        self.result_msd = self.result_msd_xyz[:,0]+self.result_msd_xyz[:,1]+self.result_msd_xyz[:,2]        
        """
    
    def plot_lindemann_msd(self,m_msd,average_1st_bond_length,time_log=None,png_filename='lindemann_msds_loglog.png'):
        R"""
        input: 
            m_msd: n rows of [interval m, msd(m)] . 
            
            time_log: at least n rows of [interval m], 
                here m is real seconds, frames or others 
                you want to set as time interval.

        output: msd plot

        introduction:
            
        """
        import matplotlib.pyplot as plt 
        a2 = average_1st_bond_length*average_1st_bond_length
        plt.figure()
        if time_log is None:
            plt.loglog(m_msd[:,1],m_msd[:,1]/a2/2.0)
        else:
            time_msd=time_log[:m_msd.shape[0]]
            plt.loglog(time_msd,m_msd[:,1]/a2/2.0)

        plt.title("dynamic Lindemann Mean Squared Displacement")
        plt.xlabel("$t$ (sec)")
        plt.ylabel("gamma_L$(t)$ (1)")
        #png_filename = 'msd_chips_long_loglog_'+'index5208_9'+'.png'
        plt.savefig(png_filename)#png_filename
        plt.close()

    def plot_msd_particle_wise_X(self,m_msd,time_log=None,png_filename='msds_loglog.png'):
        s = 0
        self.plot_lindemann_msd()

class displacemnt_field_2D:
    def __init__(self,txyz_stable,ax=None,fig=None):
        self.txyz_stable = txyz_stable
        if ax is None:
            fig,ax = plt.subplots()
        self.fig = fig
        self.ax = ax      

    def get_displacements(self,frame_index_end=-1,frame_index_start=0):
        init_positions = self.txyz_stable[frame_index_start]
        final_positions = self.txyz_stable[frame_index_end]
        self.displacements =  final_positions - init_positions

    def get_displacement_field_xy(self,frame_index_start=0,frame_index_end=-1,plot=False,png_filename=None,x_unit='($\sigma$)',limit=False):
        R"""
        Introduction:
            The function draws a displacement vector field from init state to final state 
            with positions at edge removed  to clean abnormal displacement vector. 
        input:
            frame_index: -1
            plot:True or False
            png_filename: 'displacement_field_xy.png'
            x_unit:'(sigma)', '(um)' or '(1)'
        Example:
            import points_analysis_2D as pa
            gsd = pa.proceed_gsd_file(simu_index=1382)
            gsd.get_displacement_field(plot=True)
        """
        self.get_displacements(frame_index_end,frame_index_start)

        xy = self.txyz_stable[frame_index_start]#init_positions
        xye = self.txyz_stable[frame_index_end]
        uv = self.displacements

        if plot:
            #self.ax.scatter(self.final_positions[:,0],self.final_positions[:,1])#final_state
            self.ax.quiver(xy[:,0],xy[:,1],uv[:,0],uv[:,1],color='orange',angles='xy', scale_units='xy', scale=1,width=0.01)
            self.ax.scatter(xye[:,0],xye[:,1],c='k')#init_state
            self.ax.set_title('displacement field ')#+'index:'+str(self.simu_index)
            self.ax.set_xlabel('x'+x_unit)
            self.ax.set_ylabel('y'+x_unit)
            self.ax.set_aspect('equal','box')

            if limit:#False
                self.ax.set_xlim(-4,12)
                self.ax.set_ylim(-5,16)
                pass
            if not png_filename is None:
                plt.savefig(png_filename)
            plt.close()
        return uv

    def get_displacement_field_xy_id(self,frame_index_start=0,frame_index_end=-1,x_unit='($\sigma$)',color='limegreen'):
        R"""
        Introduction:
            The function draws a displacement vector field from init state to final state 
            with positions at edge removed  to clean abnormal displacement vector. 
        input:
            frame_index: -1
            plot:True or False
            png_filename: 'displacement_field_xy.png'
            x_unit:'($\sigma$)', '(um)' or '(1)'
        return:
            ids: the ids of particles whose displacements are large enough.
        Example:
            import points_analysis_2D as pa
            gsd = pa.proceed_gsd_file(simu_index=1382)
            gsd.get_displacement_field(plot=True)
        """
        self.get_displacements(frame_index_end,frame_index_start)
        uv,ids = self.select_long_displacement_arrow()

        xy = self.txyz_stable[frame_index_start]#init_positions
        xye = self.txyz_stable[frame_index_end]


        #self.ax.scatter(self.final_positions[:,0],self.final_positions[:,1])#final_state
        self.ax.quiver(xy[ids,0],xy[ids,1],uv[:,0],uv[:,1],color=color,angles='xy', scale_units='xy', scale=1)
        """
        self.ax.scatter(xy[:,0],xy[:,1],c='k')#init_state
        self.ax.set_title('displacement field ')#+'index:'+str(self.simu_index)
        self.ax.set_xlabel('x'+x_unit)
        self.ax.set_ylabel('y'+x_unit)
        self.ax.set_aspect('equal','box')

        """
        
        if False:
            self.ax.xlim(5,20)
            self.ax.ylim(-15,0)
        return ids

    def get_displacement_field_xy_rank(self,rank,frame_index_start=0,frame_index_end=-1,x_unit='($\sigma$)'):#,color='limegreen'
        R"""
        Introduction:
            The function draws a displacement vector field from init state to final state 
            with positions at edge removed  to clean abnormal displacement vector. 
        input:
            frame_index: -1
            plot:True or False
            png_filename: 'displacement_field_xy.png'
            x_unit:'($\sigma$)', '(um)' or '(1)'
        return:
            ids: the ids of particles whose displacements are large enough.
        Example:
            import points_analysis_2D as pa
            gsd = pa.proceed_gsd_file(simu_index=1382)
            gsd.get_displacement_field(plot=True)
        """
        self.get_displacements(frame_index_end,frame_index_start)
        uv,ids = self.select_long_displacement_arrow()

        xy = self.txyz_stable[frame_index_start]#init_positions


        #self.ax.scatter(self.final_positions[:,0],self.final_positions[:,1])#final_state
        #arrow
        #plt.figure()
        #plt.quiver(xy[ids,0],xy[ids,1],uv[:,0],uv[:,1],c=rank, scale=1)
        #plt.show()
        cc = np.hypot(uv[:,0],uv[:,1])
        mapp = self.ax.quiver(xy[ids,0],xy[ids,1],uv[:,0],uv[:,1],rank,angles='xy', scale_units='xy', scale=1)#,rank
        self.fig.colorbar(mapp,ax=self.ax)#ax=self.ax

    def get_displacement_field_distribution(self,frame_index,log_mode=False,png_filename=None):
        R"""
        Introduction:
            The function draws a displacement vector field from init state to final state 
            with positions at edge removed  to clean abnormal displacement vector. 
        Example:
            import points_analysis_2D as pa
            gsd = pa.proceed_gsd_file(simu_index=1382)
            gsd.get_displacement_field(plot=True)
        """
        self.get_displacements(frame_index)

        xy = self.txyz_stable[0]#init_positions
        uv = self.displacements
        uv2 = uv*uv#let each element  i be i^2
        square_displacement = uv2[:,0]+uv2[:,1]
        rsd = np.sqrt(square_displacement)

        if log_mode:
            rsd = np.log10(square_displacement)
            plt.figure()
            count_bins=plt.hist(rsd,bins=20,range=[-2,2])#
            plt.xlabel('log(displacement)/um')
        else:
            plt.figure()
            count_bins=plt.hist(rsd,bins=20)#,range=[0,2]

        plt.savefig(png_filename)
        plt.close()

    def select_long_displacement_arrow(self):
        R"""
        output:
            uv: vectors of arrow, selected by ids only.
            ids: selected by displacements(>0.2)
        """
        #normalize arrow length part1
        uv = self.displacements
        dr2 = uv*uv
        dr = np.sqrt(dr2[:,0]+dr2[:,1])
        
        #select particles with long displacement
        dr_long_list = dr[:]>1.0#2.5(3*lcr),1.0(sigma),0.2(precise)
        ids=dr_long_list
        uv = uv[ids]

        #tune the length of arrow
        scale_arrow =1
        uv = uv*scale_arrow
        
        return uv,ids

    def compare_displacement_method(self,uv,ids,frame_index_start=0,frame_index_end=-1):
        R"""
        introduction:
            compare the displacements from start to set frame and 
            the displacements from start to end. hence a rank of frame for particles 
            who arrive destinations earlier or later
        return:
            rank_relative: an array of frame [j]. relative to the frame to start, the i-th particle arrive whose
            destination at the j-th frame.
        """
        ids_active = np.array(ids)
        num_ids_active = np.sum(ids_active.astype(int))
        rank_record = np.zeros((num_ids_active,frame_index_end-frame_index_start))
        rank = np.zeros((num_ids_active,))
        list_ids_active = np.zeros((num_ids_active,))
        #get list_ids_active
        j = 0
        for id in range(len(ids_active)):
            if ids_active[id]:
                list_ids_active[j]=id
                j=j+1
                if j==num_ids_active:
                    break
        #rank = np.zeros((np.shape(uv)[0],))

        j = 0
        for dframe in range(frame_index_end-frame_index_start):
            temp_frame_index_end = frame_index_start+dframe+1
            self.get_displacements(temp_frame_index_end,frame_index_start)
            duv = uv - self.displacements[ids,:]#
            duv2 = duv*duv
            druv = np.sqrt(duv2[:,0]+duv2[:,1])
            arrive_bool = druv[:]<1.0
            rank_record[:,j] = arrive_bool
            j=j+1
        #rank_frame = frame_index_end+1-np.sum(rank_record,1)#frame1,frame2,frame3
        rank_relative = frame_index_end-frame_index_start+1-np.sum(rank_record,1)#1,2,3
        rank_normalized = rank_relative/np.max(rank_relative)
        #print(rank_frame)
        return rank_relative#normalized

class bond_plot_module:
    def __init__(self,fig=None,ax=None):
        if ax is None:
            self.fig,self.ax = plt.subplots()
        else:
            self.fig = fig
            self.ax = ax
        
    def restrict_axis_property_relative(self,xy,x_unit='(um)'):
        R"""
        Parameters:
            txyz: all the particle positions, no one removed.
            bond_length: [particle_id1,particle_id2, bond_length] for txyz.
            check: limit the shortest and longest bond( in bond_length) to draw.
            png_filename: "prefix/bond_plot_index1513.png"
        weight of shapes:
            bond(blue line) < particles(black circle) < neighbor_change(orange circle) < traps(red cross)
            0   1   2   3   
        """
        
        #if dis is None:
        self.points = xy
        xmax = max(self.points[:,0]) #- 3
        ymax = max(self.points[:,1]) #- 3
        xmin = min(self.points[:,0]) #+ 3
        ymin = min(self.points[:,1]) #+ 3
        dis = min(xmax-xmin,ymax-ymin)/2.0#half of the length of the system.
        dis = dis - 0 #cut the edge if necessary(eg. size & scale of images not match)

        center = [(xmax+xmin)*0.5,(ymax+ymin)*0.5]
        center_origin_distance = np.abs(np.dot(center,center))
        if  center_origin_distance < 1.0:# center is really close to origin
            center = [0,0]

        #draw a figure with edges
        if x_unit == '(sigma)':
            """
            plt.rcParams.update({
            "text.usetex": True,
            "font.family": "Helvetica"
            })
            """
            x_unit = '($\sigma$)'#pip install latex is necessary for plt.savefig
        self.ax.set_aspect('equal','box')#plt.axis('equal')
        self.ax.set_xlabel('x'+x_unit)  # Add an x-label to the axes.
        self.ax.set_ylabel('y'+x_unit)  # Add a y-label to the axes.
        """
        the displayed image size will be the smaller one 
        between axis limitation for xlim/ylim or data itself. 
        Hence the effective xlim/ylim should be set smaller than data.

        To ensure that xlim/ylim < data region, 
        we add physically meaningless points.
        """
        #restrict ticks to show
        """
        #(let all the images share the same size)
        new_ticks = np.linspace(-dis,dis,int(2*dis+1))
        new_ticks = new_ticks.astype(int)
        #print(new_ticks.astype(str))
        #print((new_ticks))
        plt.xticks(new_ticks,new_ticks.astype(str))
        plt.yticks(new_ticks,new_ticks.astype(str))
        """
        #restrict data region to show
        self.ax.set_xlim(-dis+center[0],dis+center[0])#plt.xlim(-dis+center[0],dis+center[0])
        self.ax.set_ylim(-dis+center[1],dis+center[1])#plt.ylim(-dis+center[1],dis+center[1])
        self.x_unit = x_unit

    def restrict_axis_limitation(self,xlim,ylim):
        #restrict data region to show
        self.ax.set_xlim(xlim[0],xlim[1])
        self.ax.set_ylim(ylim[0],ylim[1])

    def draw_points_with_conditional_bond(self,xy,bond_length=None,bond_length_limmit=[0.9,2.0],particle_size=None):
        R"""
        Parameters:
            xy: particle positions of a frame, with no one removed.
            bond_length: [particle_id1,particle_id2, bond_length] for txyz.
            bond_length_limmit: limit the shortest and longest bond( in bond_length) to draw.
        weight of shapes:
            bond(blue line) < particles(black circle) < neighbor_change(orange circle) < traps(red cross)
            0   1   2   3   
        Examples:
            import points_analysis_2D as pa
            s = "/home/tplab/Downloads/"
            index_num = 1387
            index_name = "index"+str(index_num)
            fname = s+index_name
            bb = pa.PointsAnalysis2D(filename=fname)
            oname1 = s+"bond_hist_"+index_name
            bb.draw_bond_length_distribution_and_first_minima(png_filename=oname1)
            oname2 = s+"bond_plot_"+index_name
            bb.draw_bonds_conditional_bond(check=[0.9, bb.bond_first_minima_left], png_filename=oname2)
        """
        self.points = xy
        if not (bond_length is None):
            bond_check= tuple([bond_length_limmit[0],bond_length_limmit[1]])
            #add lines for edges
            for i in range(np.shape(bond_length)[0]):
                if (bond_length[i,2] > bond_check[0])&(bond_length[i,2] < bond_check[1]) :
                    edge = tuple(bond_length[i,0:2].astype(int))
                    pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
                    line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b',zorder=0)#,lineStyle='dashed'
                    self.ax.add_line(line)
            self.ax.set_title("bond_length:"+str(np.around(bond_length_limmit,2))+self.x_unit)  # Add a title to the axes

        if not particle_size is None:
            self.ax.scatter(xy[:,0],xy[:,1],color='k',zorder=1,s=particle_size)
        else:
            self.ax.scatter(xy[:,0],xy[:,1],color='k',zorder=1)

    def draw_points_with_conditional_vertices(self,xy,bond_length=None,vertex_bonds_index=None,particle_size=None):
        R"""
        Parameters:
            xy: particle positions of a frame, with no one removed.
            bond_length: [particle_id1,particle_id2, bond_length] for txyz.
            vertex_bonds_index: [particle_id1,particle_id2]
        weight of shapes:
            bond(blue line) < particles(black circle) < neighbor_change(orange circle) < traps(red cross)
            0   1   2   3   
        Examples:
        """
        self.points = xy
        if not (bond_length is None):
            #add lines for edges
            for i in range(np.shape(vertex_bonds_index)[0]):
                #vertex_bonds_index = vertex_bonds_index[i,0:2].astype(int)
                #bond_index=vertex_bonds_index[i]
                #pt1,pt2 = [self.points[int(bond_length[bond_index,0])],self.points[int(bond_length[bond_index,1])]]
                pt1,pt2 = [self.points[vertex_bonds_index[i,0]],self.points[vertex_bonds_index[i,1]]]
                line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b')
                plt.gca().add_line(line)
            self.ax.set_title("bond_length: vertices"+self.x_unit)  # Add a title to the axes

        if not particle_size is None:
            self.ax.scatter(xy[:,0],xy[:,1],color='k',zorder=1,s=particle_size)
        else:
            self.ax.scatter(xy[:,0],xy[:,1],color='k',zorder=1)

    def plot_neighbor_change(self,xy_stable,nb_change):
        R"""
        txyz_stable:array[Nframes,Nparticles,xyz]
                for simu data, 
                    ensure that particles never move across the boundary(box)!
                for exp data, 
                    unit of xyz must be um! 
                    ensure that particles are always in the vision field!
        nb_change: particle ids( in txyz_stable) which change neighbors.
        """
        self.ax.scatter(xy_stable[nb_change,0],xy_stable[nb_change,1],color='orange',zorder=2)
            
    def plot_traps(self,trap_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12-part1",LinearCompressionRatio=0.79,mode='array'):
        R"""
        trap_filename:
                '/home/remote/hoomd-examples_0/testhoneycomb3-8-12'
                '/home/remote/hoomd-examples_0/testhoneycomb3-8-12-part1'
                '/home/remote/hoomd-examples_0/testkagome3-11-6'
                '/home/remote/hoomd-examples_0/testkagome_part3-11-6'
        mode: 'array'(scatter) or 'map'(pcolormesh)
        """
        traps=np.loadtxt(trap_filename)
        traps=np.multiply(traps,LinearCompressionRatio)
        if mode=='array':
            #x_scale = 200
            self.ax.scatter(traps[:,0], traps[:,1],c='r',marker = 'x',zorder=3)#,s=x_scale
        elif mode=='map':
            """
            #get points
            N = 256
            vals = np.ones((N, 4))
            vals[:, 0] = np.linspace(1, 1, N)
            vals[:, 1] = np.linspace(1, 128/256, N)
            vals[:, 2] = np.linspace(1, 128/256, N)
            newcmp = ListedColormap(vals)#LinearSegmentedColormap(vals)#ListedColormap(vals)
            
            cmp = plt.get_cmap('autumn')
            cmp.reversed('autumn_r')
            """
            rcut=1.0
            cmap_name = 'Reds'#'autumn_r'#'autumn'#newcmp#'binary'#
            transparency = 0.5#0.3


            #set traps
            max = np.max(traps)
            min = np.min(traps)
            length = (max - min)
            steps = length/(rcut/10.0)
            #plt.style.use('_mpl-gallery-nogrid')

            # make data
            X, Y = np.meshgrid(np.linspace(min, max, steps.astype(int)), np.linspace(min, max, steps.astype(int)))
            HarmonicK = 100
            #origin = np.zeros((1,2))
            sz = np.shape(traps)
            i = 0
            Z = ( (0.50*HarmonicK*rcut*rcut-0.50*HarmonicK*((X-traps[i,0])**2 + (Y-traps[i,1])**2))\
                *(((X-traps[i,0])**2 + (Y-traps[i,1])**2) < rcut*rcut) )
            i = i+1
            while i<sz[0]:#sz[0]
                Zi = (0.50*HarmonicK*rcut*rcut-0.50*HarmonicK*((X-traps[i,0])**2 + (Y-traps[i,1])**2))\
                    *(((X-traps[i,0])**2 + (Y-traps[i,1])**2) < rcut*rcut)
                Z = Z + Zi
                i = i+1
            
            self.ax.pcolormesh(X, Y, Z,cmap=cmap_name,zorder = -1,alpha=transparency)#,zorder=1

    def save_figure(self,png_filename):
        R"""
        parameter:
            png_filename: "prefix/bond_plot_index1513.png"

        if latex is used in matplolib,
        'pip install latex' is necessary for plt.savefig()
        """
        self.fig.savefig(png_filename)#plt.savefig(png_filename)
        plt.close() # closes the current active figure
        #del self.ax,self.fig

    def draw_bonds_conditional_bond(self,xy,bond_length,bond_length_limmit=[0.9,2.0],x_unit='(um)'):#compatible module
        R"""
        Parameters:
            txyz: all the particle positions, no one removed.
            bond_length: [particle_id1,particle_id2, bond_length] for txyz.
            check: limit the shortest and longest bond( in bond_length) to draw.
            png_filename: "prefix/bond_plot_index1513.png"
        weight of shapes:
            bond(blue line) < particles(black circle) < neighbor_change(orange circle) < traps(red cross)
            0   1   2   3   
        Examples:
            import points_analysis_2D as pa
            s = "/home/tplab/Downloads/"
            index_num = 1387
            index_name = "index"+str(index_num)
            fname = s+index_name
            bb = pa.PointsAnalysis2D(filename=fname)
            oname1 = s+"bond_hist_"+index_name
            bb.draw_bond_length_distribution_and_first_minima(png_filename=oname1)
            oname2 = s+"bond_plot_"+index_name
            bb.draw_bonds_conditional_bond(check=[0.9, bb.bond_first_minima_left], png_filename=oname2)
        """
        bond_check= tuple([bond_length_limmit[0],bond_length_limmit[1]])
        
        #if dis is None:
        self.points = xy
        xmax = max(self.points[:,0]) #- 3
        ymax = max(self.points[:,1]) #- 3
        xmin = min(self.points[:,0]) #+ 3
        ymin = min(self.points[:,1]) #+ 3
        dis = min(xmax-xmin,ymax-ymin)/2.0#half of the length of the system.
        dis = dis - 0 #cut the edge if necessary(eg. size & scale of images not match)

        center = [(xmax+xmin)*0.5,(ymax+ymin)*0.5]
        center_origin_distance = np.abs(np.dot(center,center))
        if  center_origin_distance < 1.0:# center is really close to origin
            center = [0,0]

        #draw a figure with edges
        if x_unit == '(sigma)':
            """
            plt.rcParams.update({
            "text.usetex": True,
            "font.family": "Helvetica"
            })
            """
            x_unit = '($\sigma$)'#pip install latex is necessary for plt.savefig
        self.ax.set_aspect('equal','box')#plt.axis('equal')
        self.ax.set_xlabel('x'+x_unit)  # Add an x-label to the axes.
        self.ax.set_ylabel('y'+x_unit)  # Add a y-label to the axes.
        self.ax.set_title("bond_length:"+str(np.around(bond_length_limmit,2))+x_unit)  # Add a title to the axes
        """
        the displayed image size will be the smaller one 
        between axis limitation for xlim/ylim or data itself. 
        Hence the effective xlim/ylim should be set smaller than data.

        To ensure that xlim/ylim < data region, 
        we add physically meaningless points.
        """
        #restrict ticks to show
        """
        #(let all the images share the same size)
        new_ticks = np.linspace(-dis,dis,int(2*dis+1))
        new_ticks = new_ticks.astype(int)
        #print(new_ticks.astype(str))
        #print((new_ticks))
        plt.xticks(new_ticks,new_ticks.astype(str))
        plt.yticks(new_ticks,new_ticks.astype(str))
        """
        #restrict data region to show
        self.ax.set_xlim(-dis+center[0],dis+center[0])#plt.xlim(-dis+center[0],dis+center[0])
        self.ax.set_ylim(-dis+center[1],dis+center[1])#plt.ylim(-dis+center[1],dis+center[1])
        
        #add lines for edges
        for i in range(np.shape(bond_length)[0]):
            if (bond_length[i,2] > bond_check[0])&(bond_length[i,2] < bond_check[1]) :
                edge = tuple(bond_length[i,0:2].astype(int))
                pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
                line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b',zorder=0)
                self.ax.add_line(line)
        """
        particle_ids = np.linspace(0,self.points.shape[0]-1,self.points.shape[0],dtype=int)
        particle_ids_str = particle_ids.astype(str)
        for j in particle_ids:
            plt.annotate(particle_ids_str[j],self.points[j])
        """
        self.ax.scatter(xy[:,0],xy[:,1],color='k',zorder=1)

class bond_plot_module_for_image:
    def __init__(self,image):
        self.fig,self.ax = plt.subplots()

        self.ax.imshow(image,zorder=-1)
        
    def restrict_axis_property_relative(self,xy,x_unit='(um)'):
        R"""
        Parameters:
            txyz: all the particle positions, no one removed.
            bond_length: [particle_id1,particle_id2, bond_length] for txyz.
            check: limit the shortest and longest bond( in bond_length) to draw.
            png_filename: "prefix/bond_plot_index1513.png"
        weight of shapes:
            bond(blue line) < particles(black circle) < neighbor_change(orange circle) < traps(red cross)
            0   1   2   3   
        """
        
        #if dis is None:
        self.points = xy
        xmax = max(self.points[:,0]) #- 3
        ymax = max(self.points[:,1]) #- 3
        xmin = min(self.points[:,0]) #+ 3
        ymin = min(self.points[:,1]) #+ 3
        dis = min(xmax-xmin,ymax-ymin)/2.0#half of the length of the system.
        dis = dis - 0 #cut the edge if necessary(eg. size & scale of images not match)

        center = [(xmax+xmin)*0.5,(ymax+ymin)*0.5]
        center_origin_distance = np.abs(np.dot(center,center))
        if  center_origin_distance < 1.0:# center is really close to origin
            center = [0,0]

        #draw a figure with edges
        if x_unit == '(sigma)':
            """
            plt.rcParams.update({
            "text.usetex": True,
            "font.family": "Helvetica"
            })
            """
            x_unit = '($\sigma$)'#pip install latex is necessary for plt.savefig
        self.ax.set_aspect('equal','box')#plt.axis('equal')
        self.ax.set_xlabel('x'+x_unit)  # Add an x-label to the axes.
        self.ax.set_ylabel('y'+x_unit)  # Add a y-label to the axes.
        """
        the displayed image size will be the smaller one 
        between axis limitation for xlim/ylim or data itself. 
        Hence the effective xlim/ylim should be set smaller than data.

        To ensure that xlim/ylim < data region, 
        we add physically meaningless points.
        """
        #restrict ticks to show
        """
        #(let all the images share the same size)
        new_ticks = np.linspace(-dis,dis,int(2*dis+1))
        new_ticks = new_ticks.astype(int)
        #print(new_ticks.astype(str))
        #print((new_ticks))
        plt.xticks(new_ticks,new_ticks.astype(str))
        plt.yticks(new_ticks,new_ticks.astype(str))
        """
        #restrict data region to show
        self.ax.set_xlim(-dis+center[0],dis+center[0])#plt.xlim(-dis+center[0],dis+center[0])
        self.ax.set_ylim(-dis+center[1],dis+center[1])#plt.ylim(-dis+center[1],dis+center[1])
        self.x_unit = x_unit

    def restrict_axis_limitation(self,xlim,ylim):
        #restrict data region to show
        self.ax.set_xlim(xlim[0],xlim[1])
        self.ax.set_ylim(ylim[0],ylim[1])

    def draw_points_with_conditional_bond(self,xy,bond_length=None,bond_length_limmit=[0.9,2.0],particle_size=None):
        R"""
        Parameters:
            xy: particle positions of a frame, with no one removed.
            bond_length: [particle_id1,particle_id2, bond_length] for txyz.
            bond_length_limmit: limit the shortest and longest bond( in bond_length) to draw.
        weight of shapes:
            bond(blue line) < particles(black circle) < neighbor_change(orange circle) < traps(red cross)
            0   1   2   3   
        Examples:
            import points_analysis_2D as pa
            s = "/home/tplab/Downloads/"
            index_num = 1387
            index_name = "index"+str(index_num)
            fname = s+index_name
            bb = pa.PointsAnalysis2D(filename=fname)
            oname1 = s+"bond_hist_"+index_name
            bb.draw_bond_length_distribution_and_first_minima(png_filename=oname1)
            oname2 = s+"bond_plot_"+index_name
            bb.draw_bonds_conditional_bond(check=[0.9, bb.bond_first_minima_left], png_filename=oname2)
        """
        if not (bond_length is None):
            bond_check= tuple([bond_length_limmit[0],bond_length_limmit[1]])
            #add lines for edges
            for i in range(np.shape(bond_length)[0]):
                if (bond_length[i,2] > bond_check[0])&(bond_length[i,2] < bond_check[1]) :
                    edge = tuple(bond_length[i,0:2].astype(int))
                    pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
                    line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b',zorder=0,lw=1)#,lineStyle='dashed'
                    self.ax.add_line(line)
            self.ax.set_title("bond_length:"+str(np.around(bond_length_limmit,2))+self.x_unit)  # Add a title to the axes

        self.points = xy
        if not particle_size is None:
            self.ax.scatter(xy[:,0],xy[:,1],color='k',zorder=1,s=particle_size)
        else:
            pass#self.ax.scatter(xy[:,0],xy[:,1],color='k',zorder=1)
        return line

    def plot_neighbor_change(self,xy_stable,nb_change):
        R"""
        txyz_stable:array[Nframes,Nparticles,xyz]
                for simu data, 
                    ensure that particles never move across the boundary(box)!
                for exp data, 
                    unit of xyz must be um! 
                    ensure that particles are always in the vision field!
        nb_change: particle ids( in txyz_stable) which change neighbors.
        """
        self.ax.scatter(xy_stable[nb_change,0],xy_stable[nb_change,1],color='orange',zorder=2)
            
    def plot_traps(self,trap_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12-part1",LinearCompressionRatio=0.79,mode='array'):
        R"""
        trap_filename:
                '/home/remote/hoomd-examples_0/testhoneycomb3-8-12'
                '/home/remote/hoomd-examples_0/testhoneycomb3-8-12-part1'
                '/home/remote/hoomd-examples_0/testkagome3-11-6'
                '/home/remote/hoomd-examples_0/testkagome_part3-11-6'
        mode: 'array'(scatter) or 'map'(pcolormesh)
        """
        traps=np.loadtxt(trap_filename)
        traps=np.multiply(traps,LinearCompressionRatio)
        if mode=='array':
            #x_scale = 200
            self.ax.scatter(traps[:,0], traps[:,1],c='r',marker = 'x',zorder=3)#,s=x_scale
        elif mode=='map':
            """
            #get points
            N = 256
            vals = np.ones((N, 4))
            vals[:, 0] = np.linspace(1, 1, N)
            vals[:, 1] = np.linspace(1, 128/256, N)
            vals[:, 2] = np.linspace(1, 128/256, N)
            newcmp = ListedColormap(vals)#LinearSegmentedColormap(vals)#ListedColormap(vals)
            
            cmp = plt.get_cmap('autumn')
            cmp.reversed('autumn_r')
            """
            rcut=1.0
            cmap_name = 'Reds'#'autumn_r'#'autumn'#newcmp#'binary'#
            transparency = 0.5#0.3


            #set traps
            max = np.max(traps)
            min = np.min(traps)
            length = (max - min)
            steps = length/(rcut/10.0)
            #plt.style.use('_mpl-gallery-nogrid')

            # make data
            X, Y = np.meshgrid(np.linspace(min, max, steps.astype(int)), np.linspace(min, max, steps.astype(int)))
            HarmonicK = 100
            #origin = np.zeros((1,2))
            sz = np.shape(traps)
            i = 0
            Z = ( (0.50*HarmonicK*rcut*rcut-0.50*HarmonicK*((X-traps[i,0])**2 + (Y-traps[i,1])**2))\
                *(((X-traps[i,0])**2 + (Y-traps[i,1])**2) < rcut*rcut) )
            i = i+1
            while i<sz[0]:#sz[0]
                Zi = (0.50*HarmonicK*rcut*rcut-0.50*HarmonicK*((X-traps[i,0])**2 + (Y-traps[i,1])**2))\
                    *(((X-traps[i,0])**2 + (Y-traps[i,1])**2) < rcut*rcut)
                Z = Z + Zi
                i = i+1
            
            self.ax.pcolormesh(X, Y, Z,cmap=cmap_name,zorder = -1,alpha=transparency)#,zorder=1

    def save_figure(self,png_filename):
        R"""
        parameter:
            png_filename: "prefix/bond_plot_index1513.png"

        if latex is used in matplolib,
        'pip install latex' is necessary for plt.savefig()
        """
        self.fig.savefig(png_filename)#plt.savefig(png_filename)
        plt.close() # closes the current active figure
        #del self.ax,self.fig

    def draw_bonds_conditional_bond(self,xy,bond_length,bond_length_limmit=[0.9,2.0],x_unit='(um)'):#compatible module
        R"""
        Parameters:
            txyz: all the particle positions, no one removed.
            bond_length: [particle_id1,particle_id2, bond_length] for txyz.
            check: limit the shortest and longest bond( in bond_length) to draw.
            png_filename: "prefix/bond_plot_index1513.png"
        weight of shapes:
            bond(blue line) < particles(black circle) < neighbor_change(orange circle) < traps(red cross)
            0   1   2   3   
        Examples:
            import points_analysis_2D as pa
            s = "/home/tplab/Downloads/"
            index_num = 1387
            index_name = "index"+str(index_num)
            fname = s+index_name
            bb = pa.PointsAnalysis2D(filename=fname)
            oname1 = s+"bond_hist_"+index_name
            bb.draw_bond_length_distribution_and_first_minima(png_filename=oname1)
            oname2 = s+"bond_plot_"+index_name
            bb.draw_bonds_conditional_bond(check=[0.9, bb.bond_first_minima_left], png_filename=oname2)
        """
        bond_check= tuple([bond_length_limmit[0],bond_length_limmit[1]])
        
        #if dis is None:
        self.points = xy
        xmax = max(self.points[:,0]) #- 3
        ymax = max(self.points[:,1]) #- 3
        xmin = min(self.points[:,0]) #+ 3
        ymin = min(self.points[:,1]) #+ 3
        dis = min(xmax-xmin,ymax-ymin)/2.0#half of the length of the system.
        dis = dis - 0 #cut the edge if necessary(eg. size & scale of images not match)

        center = [(xmax+xmin)*0.5,(ymax+ymin)*0.5]
        center_origin_distance = np.abs(np.dot(center,center))
        if  center_origin_distance < 1.0:# center is really close to origin
            center = [0,0]

        #draw a figure with edges
        if x_unit == '(sigma)':
            """
            plt.rcParams.update({
            "text.usetex": True,
            "font.family": "Helvetica"
            })
            """
            x_unit = '($\sigma$)'#pip install latex is necessary for plt.savefig
        self.ax.set_aspect('equal','box')#plt.axis('equal')
        self.ax.set_xlabel('x'+x_unit)  # Add an x-label to the axes.
        self.ax.set_ylabel('y'+x_unit)  # Add a y-label to the axes.
        self.ax.set_title("bond_length:"+str(np.around(bond_length_limmit,2))+x_unit)  # Add a title to the axes
        """
        the displayed image size will be the smaller one 
        between axis limitation for xlim/ylim or data itself. 
        Hence the effective xlim/ylim should be set smaller than data.

        To ensure that xlim/ylim < data region, 
        we add physically meaningless points.
        """
        #restrict ticks to show
        """
        #(let all the images share the same size)
        new_ticks = np.linspace(-dis,dis,int(2*dis+1))
        new_ticks = new_ticks.astype(int)
        #print(new_ticks.astype(str))
        #print((new_ticks))
        plt.xticks(new_ticks,new_ticks.astype(str))
        plt.yticks(new_ticks,new_ticks.astype(str))
        """
        #restrict data region to show
        self.ax.set_xlim(-dis+center[0],dis+center[0])#plt.xlim(-dis+center[0],dis+center[0])
        self.ax.set_ylim(-dis+center[1],dis+center[1])#plt.ylim(-dis+center[1],dis+center[1])
        
        #add lines for edges
        for i in range(np.shape(bond_length)[0]):
            if (bond_length[i,2] > bond_check[0])&(bond_length[i,2] < bond_check[1]) :
                edge = tuple(bond_length[i,0:2].astype(int))
                pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
                line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b',zorder=0)
                self.ax.add_line(line)
        """
        particle_ids = np.linspace(0,self.points.shape[0]-1,self.points.shape[0],dtype=int)
        particle_ids_str = particle_ids.astype(str)
        for j in particle_ids:
            plt.annotate(particle_ids_str[j],self.points[j])
        """
        self.ax.scatter(xy[:,0],xy[:,1],color='k',zorder=1)

class dynamical_facilitation_module:
    R"""
    example:
        import points_analysis_2D as pa
        df = pa.dynamical_facilitation_module()
        df.plot_displacement_t(0)
    """
    def __init__(self):
        pass
    
    def plot_displacement_t(self,trajectory,ax=None):
        R"""
        input:    
            trajectory: Nframes rows of [x,y,z], a trajectory for a particle
        """
        uv = trajectory - trajectory[0]
        dr2 = uv*uv
        dr = np.sqrt(dr2[:,0]+dr2[:,1])
        if ax is None:
            fig,ax = plt.subplots()
            ax.set_xlabel('t (k steps)')
            ax.set_ylabel('$\delta$r ($\sigma$)')
        ax.plot(dr)
        #plt.show()
        return ax
        """
        #test right
        traj = np.ones((10,3))
        t0 = [1,2,3]
        t0 = np.array(t0)
        dt = traj - t0
        print(dt)
        """
    def plot_scan_displacement_t(self,trajectory,dt=1,ax=None):
        R"""
        input:    
            trajectory: Nframes rows of [x,y,z], a trajectory for a particle
            dt: time interval to calculate displacement.
        """
        dr = self.scan_displacement_t(trajectory,dt)
        if ax is None:
            fig,ax = plt.subplots()
            ax.set_xlabel('t (k steps)')
            ax.set_ylabel('$\delta$r ($\sigma$)')
        ax.plot(dr)
        #plt.show()
        return ax
    
    def scan_displacement_t(self,trajectory,dt=1):
        R"""
        input:    
            trajectory: Nframes rows of [x,y,z], a trajectory for a particle
            dt: time interval to calculate displacement.
        """
        uv = trajectory[dt:] - trajectory[:-dt]
        dr2 = uv*uv
        dr = np.sqrt(dr2[:,0]+dr2[:,1])

        return dr

    def compute_atmsd_scan_t(self):#reference
        R"""
        method:
            scan t axis, 1 frame for each interval. here may exist some overlaps between two intervals who share the same length of time interval m.
        record_msd:[interval m, msd_m]
        k: the frame to start scan 
        m: time interval
        m_max: max time interval
        """
        m_max = 9000#int(0.5*self.frames)#to ensure the robustness of statistics
        #k_max = 0
        start = 1000
        self.record_msd=np.zeros((m_max-start+1,2))#[interval m, msd_m]
        
        list_m=np.arange(start,m_max+1)
        self.record_msd[:,0]=list_m
        for m in list_m:
            self.dm_xyz = self.txyz_stable[m:,:,:] - self.txyz_stable[:-m,:,:]
            sq_dm_xyz2 = np.square(self.dm_xyz)
            sum_sq_dm_r2 = np.sum(sq_dm_xyz2)
            self.record_msd[m-start,1]=sum_sq_dm_r2/self.dm_xyz.shape[0]/self.particles

    def plot_displacement_t_all(self,txyz_stable,dt=1,ax=None):#[x]
        R"""
        input:    
            trajectory: Nframes rows of [x,y,z], a trajectory for a particle
        """
        
        uv = txyz_stable[dt:] - txyz_stable[:-dt]#[Nframe,Nparticle,xyz]
        dr2 = uv*uv
        sz = np.shape(dr2)
        print(sz)
        dr = np.sqrt(dr2[:,0]+dr2[:,1])
    
    def get_activation_event(self,trajectory,ta,a):#[x]
        R"""
        input:
            trajectory: 
        output:

        """
        dt = 0
        t = 0
        sz = np.size(trajectory)#(Nframe,Ndimension)
        tt = ta/2 - dt
        uv = trajectory[t+tt] - trajectory[t-tt]
        dr2 = uv*uv
        dr = np.sqrt(dr2[:,0]+dr2[:,1])
        abs(dr[t+tt] - dr[t-tt])

    def coarse_grain(self):#[x]
        pass

    def get_facilitation_waiting_time_distributions_1(self,reference_occupation,prefix=None,io_only=False):#,png_filename=None
        R"""
        input:
            reference_occupation:
                (Nframe,Nreference,)[frame:int, reference_point: int, isvacancy:bool]
                Check if a reference point is a vacancy; 
                If so, check if its neighbors are occupied in future frames 
            prefix: directory to save file.
            io_only: not proceeding data, just draw a figure as result.
        output:
            txt file: 'list_waiting_time_1.txt' e.g.: 11111222334
        example:
            import points_analysis_2D as pa
            import numpy as np
            #get reference positions
            dfm = pa.dynamical_facilitation_module()
            pt_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12"
            pt = np.loadtxt(pt_filename)
            pt = pt[:,:2]*0.81
            #get box from gsd
            gsd_filename=None#"/home/remote/hoomd-examples_0/trajectory_auto4302_9.gsd"
            pgsd=pa.proceed_gsd_file(gsd_filename,'remote',4302,9)
            box = pgsd.box[:2]
            #get reference positions inside box.
            pt = dfm.cut_edge_of_positions(pt,box)


            txyz_filename="/home/remote/Downloads/4302_9/txyz_stable.npy"
            result_filename="/home/remote/Downloads/4302_9/reference_occupation_cut.npy"
            txyz = np.load(txyz_filename)
            dfm.get_reference_occupation(pt,txyz,result_filename,1)#[:1000]
            reference_occupation=np.load(result_filename)
            prefix = "/home/remote/Downloads/4302_9/reference_occupation/"
            dfm.get_reference_neighbors(pt)
            dfm.get_facilitation_waiting_time_distributions_1(reference_occupation,prefix)#,io_only=True
            #non-square figure not restrained ax.set_xlim(-dis+center[0],dis+center[0]);ax.set_ylim(-dis+center[1],dis+center[1])
            #dfm.plot_reference_occupation(pt,reference_occupation[:106],prefix)
        """
        if not io_only:
            sz = np.shape(reference_occupation)
            list_waiting_time = np.zeros(sz)
            #list_waiting_time = np.ones(sz)
            #list_waiting_time[:] = -list_waiting_time[:]

            for frame in range(sz[0]):
                for id in range(sz[1]):
                    if reference_occupation[frame,id]:
                        nb = self.search_reference_neighbors(id)#search_nb.
                        count_nb = len(nb)
                        list_waiting_time_fr_id = np.zeros((count_nb,))
                        for i in range(count_nb):
                            bool0 = bool(reference_occupation[frame,nb[i]])
                            for frame_s in range(sz[0]):
                                if frame_s > frame:
                                    bool1 = bool(reference_occupation[frame_s,nb[i]])
                                    if (bool0 != bool1):#flop event
                                        list_waiting_time_fr_id[i] = frame_s-frame
                                        break
                        
                        if np.max(list_waiting_time_fr_id)==0:
                            waiting_time_fr_id = 0
                        elif np.min(list_waiting_time_fr_id)==0:
                            #choose the time > 0
                            #print('error! waiting time is 0!')
                            list_waiting_time_fr_id_sort = np.sort(list_waiting_time_fr_id)
                            #from small to large
                            for wt in list_waiting_time_fr_id_sort:
                                if wt > 0:
                                    waiting_time_fr_id = wt
                                    break
                        elif np.min(list_waiting_time_fr_id)>0:
                            waiting_time_fr_id = np.min(list_waiting_time_fr_id)
                        list_waiting_time[frame,id] = waiting_time_fr_id
            lwt_filename = prefix + 'list_waiting_time_frame_id_1.npy'
            np.save(lwt_filename,list_waiting_time)
            sz = np.shape(list_waiting_time)
            list_waiting_time_1d = np.reshape(list_waiting_time,sz[0]*sz[1])
            list_waiting_time_sort = np.sort(list_waiting_time_1d)
            count_wt_max = sz[0]*sz[1]
            print(count_wt_max)
            for i in range(count_wt_max):
                if list_waiting_time_sort[i]>0:#The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
                    #count_wt_1 = i
                    list_waiting_time_eff = list_waiting_time_sort[i:]
                    break
            #count_wt = count_wt_max - count_wt_1
            txt_filename = prefix + 'list_waiting_time_1.txt'
            np.savetxt(txt_filename,list_waiting_time_eff)
        else:
            txt_filename = prefix + 'list_waiting_time_1.txt'
            list_waiting_time_eff = np.loadtxt(txt_filename)

        if not prefix is None:
            fig,ax = plt.subplots()
            count_bins = ax.hist (list_waiting_time_eff,bins=200,log=True)#bins=20,range=[0,100]
            #plt.show()
            _count=count_bins[0]
            _bins=count_bins[1]
            fig2,ax2 = plt.subplots()
            ax2.loglog(_bins[1:],_count)#semilogy,loglog
            ax2.set_xlabel('waiting time (k steps)')
            ax2.set_ylabel('count (1)')
            png_filename = prefix + 'list_waiting_time_1_bin.png'
            plt.savefig(png_filename)
            txt_filename = prefix + 'list_waiting_time_1_bin.txt'
            cb = np.zeros((_count.size,2))
            cb[:,0] = _bins[1:]
            cb[:,1] = _count
            np.savetxt(txt_filename,cb)


    def get_facilitation_waiting_time_distributions_0(self,search_nb,reference_occupation,prefix=None):#[x]
        R"""
        flopping simultaneously (coarse grained step) is counted in waiting time
        """
        pass

    def cut_facilitation_waiting_time(self,list_waiting_time_frame_id,wait_max,prefix=None):#[x]
        R"""
        flopping simultaneously (coarse grained step) is counted in waiting time
        """
        sz = np.shape(list_waiting_time_frame_id)

        for frame in range(sz[0]):
            for id in range(sz[1]):
                if list_waiting_time_frame_id[frame,id]>=wait_max:
                    list_waiting_time_frame_id[frame,id]=0
        return list_waiting_time_frame_id

    def get_reference_neighbors(self,reference_positions):
        R"""
        Input: 
            reference_positions(Nreference,2)[x,y]
        Output: 
            activate search_reference_neighbors()
            List_reference_neighbors (Nreference,)[list_neighbors_id]

        example:
            import numpy as np
            import data_analysis_cycle as dac
            trap_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12-part1"
            trap = np.loadtxt(trap_filename)
            pt_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12"
            pt = np.loadtxt(pt_filename)
        """
        self.get_nb = static_points_analysis_2d(points=reference_positions)
        self.get_nb.get_first_minima_bond_length_distribution(lattice_constant=3)
        #get_nb.draw_bonds_conditional_bond()
        self.get_nb.get_coordination_number_conditional()#__coordination_bond recorded
        
    def search_reference_neighbors(self,id):
        R"""
        input:
            id
        output:
            neighobrs
        caution: function get_reference_neighbors has been runned.
        """
        neighbors = self.get_nb.search_neighbor_id(id)
        return neighbors
        
    def get_reference_occupation(self,reference_positions,txyz,result_filename=None,r_trap=1):#,txyz
        R"""
        Input: 
            txyz(Nframe,Nparticle,2)[x,y]
        Output: 
            reference_occupation: save .npy file recording list_reference_occupation
        Variables:
            list_reference_occupation (Nframe,Nreference,)[frame:int, reference_point: int, occupied:bool]        
            bond_length: (Nref,Nxy)[Length_ij]
            ref_dis: (Nreference,)
        example:
            import numpy as np
            import points_analysis_2D as pa
            pt_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12"
            pt = np.loadtxt(pt_filename)
            dfm = pa.dynamical_facilitation_module()
            dfm.get_reference_neighbors(pt)`
        """
        txy = txyz[:,:,:2]
        reference_positions = reference_positions[:,:2]
        sz_txy = np.shape(txy)
        sz_ref = np.shape(reference_positions)
        list_reference_occupation = np.zeros((sz_txy[0],sz_ref[0]))
        for frame in range(sz_txy[0]):
            xy = txy[frame]
            self.bond_length_temp=distance.cdist(reference_positions,xy,'euclidean')
            #print(self.bond_length_temp)
            
            ref_dis = np.zeros((sz_ref[0],))
            for irow in range(sz_ref[0]):
                min_irow = np.min(self.bond_length_temp[irow])
                ref_dis[irow] = min_irow
            #print(ref_dis)
            ref_dis[:] = ref_dis[:]>r_trap# 0 for committment; 1 for active
            #print(ref_dis)
            list_reference_occupation[frame] = ref_dis
        print('size of list_reference_occupation:')
        print(np.shape(list_reference_occupation))
        if not result_filename is None:
            np.save(result_filename,list_reference_occupation)
        #np.min()
        #print()
        """
        shp=np.shape(self.voronoi.ridge_points)
        self.bond_length = np.zeros((shp[0],shp[1]+1))#(start_point,end_point,length)
        self.bond_length_temp=distance.cdist(self.voronoi.points[self.voronoi.ridge_points[:,0],:], self.voronoi.points[self.voronoi.ridge_points[:,1],:], 'euclidean')
        self.bond_length[:,0:2]=self.voronoi.ridge_points
        self.bond_length[:,2]=self.bond_length_temp.diagonal()
        """
    
    def plot_reference_occupation(self,reference_positions,reference_occupation,prefix=None):
        sz = np.shape(reference_occupation)
        for frame in range(sz[0]):
            fig,ax = plt.subplots()
            circle_size=50
            vac = reference_occupation[frame,:]==1
            occu = np.logical_not(vac) 
            lw=2
            ax.scatter(reference_positions[vac,0],reference_positions[vac,1],facecolors='none',edgecolors='k',marker='o',s=circle_size,linewidths=lw)
            ax.scatter(reference_positions[occu,0],reference_positions[occu,1],c='k',marker='o',s=circle_size)
            
            png_filename = prefix+'reference_occupation'+str(int(frame))+'.png'
            #scatter size
            #axis_equal
            x_unit = '($\sigma$)'#pip install latex is necessary for plt.savefig
            ax.set_aspect('equal','box')#plt.axis('equal')
            ax.set_xlabel('x'+x_unit)  # Add an x-label to the axes.
            ax.set_ylabel('y'+x_unit)  # Add a y-label to the axes.
            plt.savefig(png_filename)
            plt.close()
    
    def plot_reference_positions_waiting_time_colored(self,reference_positions,reference_occupation,prefix=None):
        R"""
        EXAMPLE:
            #get box from gsd
            gsd_filename=None#"/home/remote/hoomd-examples_0/trajectory_auto4302_9.gsd"
            pgsd=pa.proceed_gsd_file(gsd_filename,'remote',4302,9)
            box = pgsd.box[:2]
            pt_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12"
            pt = np.loadtxt(pt_filename)
            pt = pt[:,:2]*0.81
            dfm = pa.dynamical_facilitation_module()
            pt = dfm.cut_edge_of_positions(pt,box)
            print(np.shape(pt))
            txyz_filename="/home/remote/Downloads/4302_9/txyz_stable.npy"
            result_filename="/home/remote/Downloads/4302_9/reference_occupation_cut.npy"
            txyz = np.load(txyz_filename)
            #dfm.get_reference_occupation(pt,txyz,result_filename,1)#[:1000]
            reference_occupation=np.load(result_filename)
            prefix = "/home/remote/Downloads/4302_9/reference_occupation/"
            #dfm.get_reference_neighbors(pt)
            #dfm.get_facilitation_waiting_time_distributions_1(reference_occupation,prefix)#,io_only=True
            #non-square figure not restrained ax.set_xlim(-dis+center[0],dis+center[0]);ax.set_ylim(-dis+center[1],dis+center[1])
            dfm.plot_reference_occupation(pt,reference_occupation[:110],prefix)

            #draw scatter reference positions waiting time colored
            lwt_filename = prefix + 'list_waiting_time_frame_id_1.npy'
            list_waiting_time_frame_id = np.load(lwt_filename)
            list_waiting_time_frame_id = dfm.cut_facilitation_waiting_time(list_waiting_time_frame_id,5,prefix)
            dfm.plot_reference_positions_waiting_time_colored(pt,list_waiting_time_frame_id,prefix)
        """
        sz = np.shape(reference_occupation)
        for frame in range(sz[0]):
            if frame<=200:
                fig,ax = plt.subplots()
                circle_size=50
                mapp = ax.scatter(reference_positions[:,0],reference_positions[:,1],c=reference_occupation[frame],marker='o',s=circle_size)
                png_filename = prefix+'reference_waiting_time_color'+str(int(frame))+'.png'
                #scatter size
                #axis_equal
                x_unit = '($\sigma$)'#pip install latex is necessary for plt.savefig
                ax.set_aspect('equal','box')#plt.axis('equal')
                ax.set_xlabel('x'+x_unit)  # Add an x-label to the axes.
                ax.set_ylabel('y'+x_unit)  # Add a y-label to the axes.
                fig.colorbar(mapp,ax=ax)#plt.cm.ScalarMappable(norm=plt.autoscale()),
                #plt.colorbar()
                plt.savefig(png_filename)
                plt.close()
        
    def cut_edge_of_positions(self,points,box):
        R"""
        input:
            points: [x,y]
            box: [Lx,Ly] get from gsd file
        return:
            points_selected
        Variables:
            effective_lattice_constant： the distance between a particle and its neighbor
            edge_cut_positions_list: list the rows of particles' positions at given snapshot.
        example:
            gsd_filename="/home/remote/hoomd-examples_0/trajectory_auto4302_9.gsd"
            pgsd=pa.proceed_gsd_file(gsd_filename,'remote',4302,9)
            box = pgsd.box[:2]
        """
        xmax = box[0]/2
        ymax = box[1]/2
        xmin = -box[0]/2
        ymin = -box[1]/2
        
        #That directly using 'and' has been banned, so 'logical_and' is necessary
        list_xmin = points[:,0] > xmin
        list_xmax = points[:,0] < xmax
        list_ymin = points[:,1] > ymin
        list_ymax = points[:,1] < ymax
        list_x = np.logical_and(list_xmin,list_xmax)
        list_y = np.logical_and(list_ymin,list_ymax)
        list_xy = np.logical_and(list_x,list_y)

        #edge_cut_positions_list = np.where(list_xy)
        #edge_cut_positions_bool = list_xy 
        points_selected = points[list_xy]
        print('size of traps:')
        print(np.shape(points_selected))
        return points_selected

class polygon_analyzer:
    def __init__(self,points,bonds,faces):
        R"""
        
        """
        pass