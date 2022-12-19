import numpy 
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import  Delaunay
from scipy.spatial import  distance
import matplotlib 
import os
import freud

from particle_tracking import particle_track

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
                self.points = numpy.loadtxt(filename)
                self.points = self.points[:,0:2]
                self.points_analysis()
        else :
            self.points = points[:,:2]
            self.points_analysis()

        #not show figures
        if hide_figure:
            matplotlib.use(backend="agg")#Backend agg is non-interactive backend. Turning interactive mode off. 'QtAgg' is interactive mode
        #set basic parameters
        #self.prefix='/home/tplab/Downloads/'
            
    def points_analysis(self):
        self.voronoi = Voronoi(self.points)
        self.delaunay = Delaunay(self.points)
        self.get_ridge_length()
        self.get_bond_length()
        #cut the edge of positions
        self.__cut_edge_of_positions()#effective lattice constant is given 3 defaultly
    #def bonds_analysis(self):

    def get_ridge_length(self):
        #where the ridge length of vertex-i & vertex-j is ridge_length[i,j]
        self.voronoi.ridge_vertices=numpy.array(self.voronoi.ridge_vertices)#选不了ridge_vertices,不是数组格式
        self.voronoi.ridge_length=distance.cdist(self.voronoi.vertices[self.voronoi.ridge_vertices[:,0],:],self.voronoi.vertices[self.voronoi.ridge_vertices[:,1],:],'euclidean')
        self.voronoi.ridge_length=self.voronoi.ridge_length.diagonal()

    def get_bond_length(self):
        #self.bond_length: (start_point_index,end_point_index,bond_length)
        shp=numpy.shape(self.voronoi.ridge_points)
        self.bond_length = numpy.zeros((shp[0],shp[1]+1))#(start_point,end_point,length)
        self.bond_length_temp=distance.cdist(self.voronoi.points[self.voronoi.ridge_points[:,0],:], self.voronoi.points[self.voronoi.ridge_points[:,1],:], 'euclidean')
        self.bond_length[:,0:2]=self.voronoi.ridge_points
        self.bond_length[:,2]=self.bond_length_temp.diagonal()
        #should I list point pair like small-to-large？
        del self.bond_length_temp #delete temp data

    def get_first_minima_bond_length_distribution(self,lattice_constant=3.0,method="global_compare_method",png_filename=None,hist_cutoff=2):
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
        self.bond_sorted=numpy.sort(self.bond_length[:,2]/self.lattice_constant)
        #find the minimum bin, then set the left side as bond_first_neighbour
        i_max=self._bins.size
        i=0
        if method=="local_compare_method":
            while i < i_max-3:#i start from 0, which have to -1;compares i,i+1 and i+2, which have to -2, hence -3
                #numpy.where( self.bond_sorted[:]>res.bins[i] & self.bond_sorted[:]< res.bins[i+1] ) 
                #print(self._count[i])
                if self._count[i] > self._count[i+1]:
                    if self._count[i+1] <= self._count[i+2]:
                        if self._bins[i+1] > 1/self.lattice_constant:
                            # bond_length should be larger than sigma of particle
                            i+=1
                            break
                i+=1

        elif method=="global_compare_method":
            self._count = numpy.array(self._count)
            self.count_sorted=numpy.sort(self._count)
            """
            print(self.count_sorted)
            print(self.count_sorted[-1])#max1
            print(self.count_sorted[-2])#max2
            """
            i_bins_for_count_peak_1 = numpy.where(self._count[:]==self.count_sorted[-1])
            i_bins_for_count_peak_2 = numpy.where(self._count[:]==self.count_sorted[-2])
            i_bins_for_count_peak_1 = i_bins_for_count_peak_1[0]
            i_bins_for_count_peak_2 = i_bins_for_count_peak_2[0]
            #if there are more than one bins share the same count, select the smallest bin number.
            if numpy.shape(i_bins_for_count_peak_1)[0]>1:
                i_bins_for_count_peak_1 = min(i_bins_for_count_peak_1)
            if numpy.shape(i_bins_for_count_peak_2)[0]>1:
                i_bins_for_count_peak_2 = min(i_bins_for_count_peak_2)
            i_bins_for_count_1st_peak = min(i_bins_for_count_peak_1,i_bins_for_count_peak_2)
            i_bins_for_count_1st_peak = int(i_bins_for_count_1st_peak)#force i be a simple data_type(int)
            """
            print(i_bins_for_count_peak_1) 
            print(i_bins_for_count_peak_2)
            print(i_bins_for_count_1st_peak) 
            """
            #list_xy = numpy.logical_and(list_x,list_y)
            #self.edge_cut_positions_list = numpy.where(==)
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
        self.bond_first_neighbour=self.bond_sorted[numpy.where(self.bond_sorted[:]<self.bond_first_minima_left)]
        
        self.bond_length_median = numpy.around( numpy.median(self.bond_first_neighbour),2)
        self.bond_length_std = numpy.around(numpy.std(self.bond_first_neighbour),2)

        if not png_filename  is None:
            #plot divide line
            plt.plot([self.bond_first_minima_left,self.bond_first_minima_left],[0,self._count[i]],linestyle='--')
            plt.title("1st peak:"+str( self.bond_length_median )+"+-"+str( self.bond_length_std ))
            plt.xlabel("bond length (unit=A)")
            plt.ylabel("count (1)")
            #plt.show()
            #png_filename=prefix+"bond_hist_index"+str(index_num)
            plt.savefig(png_filename)
        plt.close()
        #let bond_first_minima_left be a absolute length, not normalized one
        self.bond_first_minima_left=self.bond_first_minima_left*self.lattice_constant

    def get_first_minima_radial_distribution_function(self,rdf,lattice_constant=3.0,png_filename=None):#[x]
        #print(rdf.bin_centers) print(rdf.bin_counts)
        self.lattice_constant=lattice_constant
        #locate the 1st minima of bond-length distribution
        self._count=rdf.bin_counts
        self._bins=rdf.bin_centers
        self.bond_sorted=numpy.sort(self.bond_length[:,2]/self.lattice_constant)
        #find the minimum bin, then set the left side as bond_first_neighbour
        
        "global_compare_method"
        self._count = numpy.array(self._count)
        self._bins = numpy.array(self._bins)
        self.count_sorted=numpy.sort(self._count)

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
        i_bins_for_count_peak_1 = numpy.where(self._count[:]==self.count_sorted[-1])
        i_bins_for_count_peak_2 = numpy.where(self._count[:]==self.count_sorted[-2])
        i_bins_for_count_peak_1 = i_bins_for_count_peak_1[0]
        i_bins_for_count_peak_2 = i_bins_for_count_peak_2[0]
        #if there are more than one bins share the same count, select the smallest bin number.
        if numpy.shape(i_bins_for_count_peak_1)[0]>1:
            i_bins_for_count_peak_1 = min(i_bins_for_count_peak_1)
        if numpy.shape(i_bins_for_count_peak_2)[0]>1:
            i_bins_for_count_peak_2 = min(i_bins_for_count_peak_2)
        i_bins_for_count_1st_peak = min(i_bins_for_count_peak_1,i_bins_for_count_peak_2)
        i_bins_for_count_1st_peak = int(i_bins_for_count_1st_peak)#force i be a simple data_type(int)
        """
        print(i_bins_for_count_peak_1) 
        print(i_bins_for_count_peak_2)
        print(i_bins_for_count_1st_peak) 
        """
        #list_xy = numpy.logical_and(list_x,list_y)
        #self.edge_cut_positions_list = numpy.where(==)
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
        self.bond_first_neighbour=self.bond_sorted[numpy.where(self.bond_sorted[:]<self.bond_first_minima_left)]
        
        if not png_filename  is None:
            #plot divide line
            plt.plot([self.bond_first_minima_left,self.bond_first_minima_left],[0,0.5],linestyle='--')#self._count[i]
            plt.title("1st peak:"+str( numpy.around( numpy.median(self.bond_first_neighbour),2) )+"+-"+str( numpy.around(numpy.std(self.bond_first_neighbour),2) ))
            plt.xlabel("interparticle distance (unit=A)")
            plt.ylabel("count (1)")
            #plt.show()
            #png_filename=prefix+"bond_hist_index"+str(index_num)
            plt.savefig(png_filename)
        plt.close()
        #let bond_first_minima_left be a absolute length, not normalized one
        self.bond_first_minima_left=self.bond_first_minima_left*self.lattice_constant
        
    def get_coordination_number_conditional(self):
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
            self.get_first_minima_bond_length_distribution()
            bond_length_limit=self.bond_first_minima_left

        #select the bonds within the limit(bond length of 1st neighbour)
        self._bond_min=0.0
        self._bond_max=bond_length_limit
        place=numpy.where((self.bond_length[:,2]>self._bond_min) & (self.bond_length[:,2]<self._bond_max))
        self.__coordination_bond=self.bond_length[place,0:2]
        self.__coordination_bond=self.__coordination_bond[0]#remove a [] in [[]] structure
        
        
        #self._particle_id_max=numpy.max(self.__coordination_bond[:])#get the largest id of particle to restrain while loop
        #print(self.__coordination_bond)

        cn_max = 12#for honeycomb, cn10 truly appears in Index4340_seed9_frame2!
        self.count_coordination=numpy.zeros([cn_max,1])
        
        for id in self.edge_cut_positions_list[0]:
            place=numpy.where(self.__coordination_bond[:]==id)
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
        
    def get_nearest_neighbor_dispalcement(self):
        R"""
        Introduction:
            lindemann_dispalcement: r(i,t)-<r(j,t)>
        Return: 
            list_id_dxy_nb: 
                list[particle_id,num_neighbors,dx,dy], 
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
        place=numpy.where((self.bond_length[:,2]>self._bond_min) & (self.bond_length[:,2]<self._bond_max))
        self.__coordination_bond=self.bond_length[place,0:2]
        self.__coordination_bond=self.__coordination_bond[0]#remove a [] in [[]] structure
        
        
        for id in self.edge_cut_positions_list[0]:
            place=numpy.where(self.__coordination_bond[:]==id)
            n_nb=len(place[0])#number of neighbors
            place_nb=numpy.array(place) 
            #print(place_nb)
            place_nb[1]=1-place_nb[1]
            #print(place_nb)
            id_nb=self.__coordination_bond[place_nb[0],place_nb[1]].astype(int)
            #print(id_nb)#id s 
            
            #print(self.points[id_nb])
            center_nb = numpy.average(self.points[id_nb],axis = 0)
            dxy_id = self.points[id]-center_nb
            #print(center_nb)
            if 'list_id_dxy_nb' in locals():
                list_id_dxy_nb.append([id,n_nb,dxy_id[0],dxy_id[1]]) 
                dict_c_nb[str(id)] = id_nb
            else:
                list_id_dxy_nb = []
                list_id_dxy_nb.append([id,n_nb,dxy_id[0],dxy_id[1]])
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
        new_ticks = numpy.linspace(-dis,dis,int(2*dis+1))
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
        #count=numpy.zeros((self.points.shape[0],1))
        for i,j in pairs:
            if self.edge_cut_positions_bool[i]:
                #numpy.np.concatenate([array1,array2]) 
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
        new_ticks = numpy.linspace(-dis,dis,int(2*dis+1))
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
        for i in range(numpy.shape(self.bond_length)[0]):
            if (self.bond_length[i,2] > check[0])&(self.bond_length[i,2] < check[1]) :
                edge = self.bond_length[i,0:2].astype(int)
                pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
                line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b')
                plt.gca().add_line(line)
        #plt.show()
        plt.title("bond_length:"+str(numpy.around(check,2))+"um")
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
        #sz = len(self.init_positions)#numpy.size
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
        list_x = numpy.logical_and(list_xmin,list_xmax)
        list_y = numpy.logical_and(list_ymin,list_ymax)
        list_xy = numpy.logical_and(list_x,list_y)

        self.edge_cut_positions_list = numpy.where(list_xy)
        self.edge_cut_positions_bool = list_xy # T for body, F for edge.


    def get_bond_orientational_order(self,plot=False,png_filename=None,k_set=6):
        R"""
        Parameters:
            k_set: set the k-th disclination order parameter psi_k.
        """
        box=numpy.zeros(2)
        x1=min(self.points[:,0])
        x2=max(self.points[:,0])
        Lx=x2-x1#get box size in x-direction
        y1=min(self.points[:,1])
        y2=max(self.points[:,1])
        Ly=y2-y1#get box size in y-direction
        box[0]=Lx+1
        box[1]=Ly+1
        hex_order = freud.order.Hexatic(k=k_set)

        sp=numpy.shape(self.points)
        pts=numpy.zeros((sp[0],sp[1]+1))
        pts[:,0:2]=self.points[:]
        hex_order.compute(system=(box,pts[:]))

        self.Psi_k=hex_order.particle_order#local BOO
        self.Psi_k_abs=abs(hex_order.particle_order)
        self.Psi_k_global_with_edge=numpy.average(self.Psi_k_abs)#abs checked right
        self.Psi_k_global_cut_edge=numpy.average(self.Psi_k_abs[self.edge_cut_positions_list])#abs checked right

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

        place=numpy.where(pk[:]>lower_limit)
        sp_3=numpy.shape(pk[place])
        sp_all=numpy.shape(pk)
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

    def draw_bonds_conditional_bond(self,check=[0.9,2.0],png_filename=None,show_traps=False,LinearCompressionRatio=0.79,trap_filename="/home/tplab/hoomd-examples_0/testhoneycomb3-8-12-part1"):
        #bond_plot+trap_plot[X]
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
        bond_check= [check[0],check[1]]
        
        #if dis is None:
        xmax = max(self.points[:,0]) #- 3
        ymax = max(self.points[:,1]) #- 3
        xmin = min(self.points[:,0]) #+ 3
        ymin = min(self.points[:,1]) #+ 3
        dis = min(xmax-xmin,ymax-ymin)/2.0#half of the length of the system.
        dis = dis - 0 #cut the edge if necessary(eg. size & scale of images not match)

        center = [(xmax+xmin)*0.5,(ymax+ymin)*0.5]
        center_origin_distance = numpy.abs(numpy.dot(center,center))
        if  center_origin_distance < 1.0:# center is really close to origin
            center = [0,0]

        #draw a figure with edges
        plt.figure()
        plt.scatter(self.points[:,0],self.points[:,1],color='k')
        if show_traps:
            traps=numpy.loadtxt(trap_filename)
            plt.scatter(traps[:,0]*LinearCompressionRatio, 
                   traps[:,1]*LinearCompressionRatio,
                   c='r',marker = 'x')
            #scatter_trap_plot[X] red cross
            pass
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
        new_ticks = numpy.linspace(-dis,dis,int(2*dis+1))
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
        for i in range(numpy.shape(self.bond_length)[0]):
            if (self.bond_length[i,2] > bond_check[0])&(self.bond_length[i,2] < bond_check[1]) :
                edge = self.bond_length[i,0:2].astype(int)
                pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
                line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b')
                plt.gca().add_line(line)
        #plt.show()
        plt.title("bond_length:"+str(numpy.around(check,2))+"um")
        if not png_filename is None:
            plt.savefig(png_filename)
        plt.close()

    def draw_bonds_conditional_ridge(self,fignum=1,check=[0.2,5.0]):
        #record cutted bond
        count=0
        shp=numpy.shape(self.voronoi.ridge_length)
        cutted_temp = numpy.zeros((shp[0],2))
        cutted_ridge_temp = numpy.zeros((shp[0],1))
        #[x]
        median_temp = numpy.median(self.voronoi.ridge_length)
        ridge_check= [check[0]*median_temp,check[1]*median_temp]
        self.voronoi.cut_ridge_length=ridge_check[0]
        #check_min can be too large, but can not be too small
        #for amorphous system, check is hard to set
        del median_temp

        #draw a figure with edges
        plt.figure(fignum)
        plt.scatter(self.points[:,0],self.points[:,1],color='k')
        plt.title(str(check)+'median')
        #add lines for edges
        for i in range(numpy.shape(self.voronoi.ridge_length)[0]):
            if (self.voronoi.ridge_length[i] > ridge_check[0])&(self.voronoi.ridge_length[i] < ridge_check[1]) :
                edge = self.voronoi.ridge_points[i]
                pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
                line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b')
                plt.gca().add_line(line)
            elif (self.voronoi.ridge_length[i] < ridge_check[0])&(self.voronoi.ridge_length[i] > 0):
                cutted_temp[count]=self.voronoi.ridge_points[i]
                cutted_ridge_temp[count]=i
                count=count+1
                
        self.cutted_bonds = cutted_temp[0:count,:]
        self.voronoi.cutted_ridges = cutted_ridge_temp[0:count]
        del cutted_temp,cutted_ridge_temp
    
    def print_benchmark(self,bond=True,ridge=True):
        if bond :
            print('bond')
            print('min',min(self.bond_length[:,2]))
            print('median',numpy.median(self.bond_length[:,2]))
            print('max',max(self.bond_length[:,2]))
            print('\n')
        if ridge :
            print('ridge')
            print('min',min(self.voronoi.ridge_length))
            print('median',numpy.median(self.voronoi.ridge_length))
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

        energy_interaction=epsilon*numpy.exp(-kappa)# when rcut=12.5, Uparticle = kBT = 1; rcut=15, Up=0.47
        energy_substrate=0.5*k*numpy.multiply(rcut,rcut)
        #pt=numpy.exp(1)
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
        c=numpy.concatenate([v_index_temp[:,0],v_index_temp[:,1]])
        c=numpy.sort(c)
        shp=numpy.shape(c)
        for j in range(shp[0]-1):
            if c[j]==c[j+1]:
                c[j+1]=-1
        c=numpy.sort(c)
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
        
        shp=numpy.shape(cv_temp)
        cluster_count=numpy.zeros((shp[0]))
        check=self.voronoi.cut_ridge_length
        cv_dis=distance.cdist(cv_temp,cv_temp,'euclidean')
        
        for i in range(shp[0]):
            '''
            if cv_dis[i,0]==255:
                continue
            else:
                list=numpy.where(cv_dis[i,:]<check)
                print(list[0])
                #cv_dis[:,list]=255#delete found vertices 
                cv_dis[list,:]=255#delete found vertices 
                count=numpy.shape(list[0])
                print(cv_dis[0:5,0:5])
                count=count[0]#how many vertices in a cluster
                cluster_count[count]=cluster_count[count]+1
            '''
            list=numpy.where(cv_dis[i,:]<1.0*check)
            cv_dis[:,list]=255#delete found vertices 
            cv_dis[list,:]=255#delete found vertices 
            count=numpy.shape(list[0])
            count=count[0]#how many vertices in a cluster
            cluster_count[count]=cluster_count[count]+1
            if count==(n_polygon-2):
                plt.annotate(cluster_count[count].astype(int), (cv_temp[i,0], cv_temp[i,1]))
        print('Tetragon,Pentagon,Hexagon...')
        print(cluster_count[2:10])
        #some great polygons whose vertices linked closely(small dis),
        # while twisted polygons whose vertices are farther from each other 
        # than from different vertices clusters.
        # so an algorithm with more precision is needed 
            
'''
def count_honeycomb(self):
    #check whether cutted_bonds exist
    if not ('self.cutted_bonds' in dir()):
        self.draw_bonds_conditional_ridge()
    #get linked triangles pair   
    shp=numpy.shape(self.cutted_bonds)
    self.linked_triangles = numpy.zeros(shp)
    for i in range(shp[0]):
        #1st round search
        index1=numpy.where(self.delaunay.simplices[:]==self.cutted_bonds[i,0])
        #2nd round search
        index2=numpy.where(self.delaunay.simplices[index1[0]]==self.cutted_bonds[i,1])
        temp=index1[0]
        self.linked_triangles[i]=temp[index2[0]]
    del temp
    #get linked triangles chain
    chain=numpy.zeros((shp[0],8))#up to 10-polygon
    chain[:]=-1
    all_chain_count=numpy.zeros((shp[0],1))
    chain_num=1#now we are proceeding the chain_num-th chain [x]unfinished
    tri_count_max=shp[0]#at most the triangle pairs should be linked
    lt_temp=self.linked_triangles[:]#numpy.zeros((shp[0],2))

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
            index1=numpy.where((lt_temp[:,0]==self.linked_triangles[i,1])|(lt_temp[:,1]==self.linked_triangles[i,1]))#rightward
            index_check1=numpy.shape(index1[0])
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
                      
'''
for i in range(shp[0]):
    if i==0:
        chain[0,0:2]=lt_temp[i]
        lt_temp[i]=[-1,-1]#clear the i-th triangle to avoid being found
        chain_count[count]=chain_count[count]+1#the count of chain-0 is added one
        count=count+1#now there is one more chain
    else:
        index1=numpy.where((chain[:]==lt_temp[i,0])|(chain[:]==lt_temp[i,1]))
        temp=index1[0]#row-index of linked triangles
        for j in temp:
            if chain_count[j,0]:
            lt_temp[j,1]
        #that chain[0] and chain[5] should be linked may exist

del lt_temp,temp
'''

'''
#https://blog.csdn.net/u013378642/article/details/81775131
self.linked_triangles_sorted=sorted(self.linked_triangles,key=lambda x:x[0])
print('t1\n',self.linked_triangles,'\nt2\n',self.linked_triangles_sorted)
#a chain of triangle could be like [1,162] [162,4] [4,88.]
'''

'''
#code for checking results
i=int(0)
...
print(self.linked_triangles[i,0].astype(int))
print(self.cutted_bonds[i])
print(self.delaunay.simplices[self.linked_triangles[i,0].astype(int)])
print(self.delaunay.simplices[self.linked_triangles[i,1].astype(int)])
'''

'''
for n in [i for i in self.delaunay.simplices]:
self.cutted_bonds[i,0]
self.cutted_bonds[i,1]
self.linked_triangles[i]
'''
'''
voronoi_plot_2d(vor)

plt.figure(1)
plt.scatter(newData[:,0],newData[:,1])


plt.figure(2)
plt.scatter(newData[:,0],newData[:,1],)
plt.show()
'''

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
                #self.data_gsd = numpy.loadtxt(self.filename_gsd)
                self.__open_gsd()

                prefix_gsd = '/home/'+account+'/hoomd-examples_0/trajectory_auto'
                simu_index = filename_gsd_seed.strip(prefix_gsd)
                id=simu_index.index('_')
                self.simu_index = simu_index[0:id]

        else :
            self.simu_index = simu_index
            if not seed is None:
                simu_index = str(simu_index)+'_'+str(seed)
            prefix_gsd = '/home/'+account+'/hoomd-examples_0/trajectory_auto'
            postfix_gsd = '.gsd'
            self.filename_gsd = prefix_gsd+str(simu_index)+postfix_gsd
            self.__open_gsd()
        
        self.box = self.trajectory.read_frame(-1).configuration.box
        self.__cut_edge_of_positions()

    def __open_gsd(self):
        import gsd.hoomd
        self.trajectory=gsd.hoomd.open(self.filename_gsd)#open a gsd file
        self.num_of_frames=len(self.trajectory)
        
    def read_a_frame(self,frame_num):
        snap=self.trajectory.read_frame(frame_num)#take a snapshot of the N-th frame
        positions=snap.particles.position[:,0:2]#just record [x,y] ignoring z
        return positions
        #self.N = self.snap.particles.N
    
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
        self.get_displacements(frame_index)

        #cut the edge of positions
        self.__cut_edge_of_positions()#effective lattice constant is given 3 defaultly
        list = self.edge_cut_positions_list
        xy = self.init_positions[list]
        uv = self.displacements[list]

        if plot:
            plt.figure()
            #plt.scatter(self.init_positions[:,0],self.init_positions[:,1])#init_state
            #plt.scatter(self.final_positions[:,0],self.final_positions[:,1])#final_state
            plt.quiver(xy[:,0],xy[:,1],uv[:,0],uv[:,1],angles='xy', scale_units='xy', scale=1)
            plt.title('displacement field '+'index:'+str(self.simu_index))
            plt.xlabel('x(sigma)')
            plt.ylabel('y(sigma)')
            if not png_filename is None:
                plt.savefig(png_filename)
            plt.close()

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

        #cut the edge of positions
        self.__cut_edge_of_positions()#effective lattice constant is given 3 defaultly
        list = self.edge_cut_positions_list
        xy = self.init_positions[list]
        uv = self.displacements[list]
        uv2 = uv*uv#let each element  i be i^2
        square_displacement = uv2[:,0]+uv2[:,1]
        rsd = numpy.sqrt(square_displacement)

        if log_mode:
            rsd = numpy.log10(square_displacement)
            plt.figure()
            count_bins=plt.hist(rsd,bins=20,range=[-2,2])#
            plt.xlabel('log(displacement)/um')
        else:
            plt.figure()
            count_bins=plt.hist(rsd,bins=20)#,range=[0,2]
        
        """
        prefix_old="/home/tplab/hoomd-examples_0"
        folder_name=self.filename_gsd.strip(prefix_old)
        folder_name=folder_name.strip(".gsd")#folder_name=prefix+png_filename_as_folder
        folder_name="t"+folder_name#"t" is necessary in case of "/t" deleted before
        prefix='/home/tplab/Downloads/'
        png_filename=prefix+folder_name+"/"+folder_name+"_"+str(frame_index)+'_'+str(int(i))+'Displacement_Field_hist'+'.png'
        #check if the folder exists
        isExists=os.path.exists(prefix+folder_name)
        if isExists:
            plt.savefig(png_filename)
        else:
            os.makedirs(prefix+folder_name)
            plt.savefig(png_filename)
        """
        plt.savefig(png_filename)
        plt.close()
        #plt.show()
        """
        #check code
        import numpy as np
        uv = np.zeros((4,2))
        i=0
        while i<4:
            uv[i]=[i,i+1]
            i+=1
        print(uv)
        uv2=uv*uv
        print(uv2)
        square_displacement = uv2[:,0]+uv2[:,1]
        print(square_displacement)
        rsd = np.sqrt(square_displacement)
        print(rsd)
        """ 
        #square_displacement = uv2[:,0]+uv2[:,1]
    def get_displacements(self,frame_index):
        self.init_positions = numpy.array(self.read_a_frame(0)) 
        final_positions = numpy.array(self.read_a_frame(frame_index)) 
        self.displacements =  final_positions - self.init_positions
        #return displacements

    def __cut_edge_of_positions(self,effective_lattice_constant = 3):
        R"""
        Variables:
            effective_lattice_constant： the distance between a particle and its neighbor
            edge_cut_positions_list: list the rows of particles' positions at given snapshot.
        """
        self.init_positions = numpy.array(self.read_a_frame(0)) 
        a = effective_lattice_constant
        #sz = len(self.init_positions)#numpy.size
        #xy = self.init_positions
        xmax = max(self.init_positions[:,0])
        ymax = max(self.init_positions[:,1])
        xmin = min(self.init_positions[:,0]) 
        ymin = min(self.init_positions[:,1])
        xmax = xmax - a
        ymax = ymax - a
        xmin = xmin + a
        ymin = ymin + a
        
        #That directly using 'and' has been banned, so 'logical_and' is necessary
        list_xmin = self.init_positions[:,0] > xmin
        list_xmax = self.init_positions[:,0] < xmax
        list_ymin = self.init_positions[:,1] > ymin
        list_ymax = self.init_positions[:,1] < ymax
        list_x = numpy.logical_and(list_xmin,list_xmax)
        list_y = numpy.logical_and(list_ymin,list_ymax)
        list_xy = numpy.logical_and(list_x,list_y)

        self.edge_cut_positions_list = numpy.where(list_xy)

    def get_coordination_number(self):#[x]
        """
        import points_analysis_2D
        obj_of_simu_index = points_analysis_2D.PointsAnalysis2D(filename=data_filename)
        obj_of_simu_index.get_coordination_number_conditional()
        ccn=obj_of_simu_index.count_coordination_ratio
        """
        pass
    
    def get_trajectory_data(self):
        R"""
        introduction:
            transform gsd file into an array [N frames,N particles,3],
            recording the trajectory of particles.
        """
        iframe = 0
        snapi = self.trajectory.read_frame(iframe)
        pos_list = numpy.zeros([self.num_of_frames,snapi.particles.N,3])#gsd_data.trajectory[0].particles.N,
        while iframe < self.num_of_frames:
            pos_list[iframe] = self.trajectory.read_frame(iframe).particles.position
            #print(self.trajectory.read_frame(iframe).configuration.box)
            iframe = iframe + 1
        
        self.txyz = pos_list
    
    def get_trajectory(self,length_cut_edge=0):#checked right
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
            new_ticks = numpy.linspace(-dis,dis,int(2*dis+1))
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
    pass

class dynamic_points_analysis_2d:#old_class_name: msd
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
        time_log = numpy.loadtxt(path+'DefaultVideo_5.txt')

        pa.plot_msd(time_log,path+'msd_chips_loglog_'+'20221129_video5_90%'+'.png')
        pa.plot_trajectory(path+'trajectory_stable_'+'20221129_video5'+'.png')

    EXAMPLE:
        print(msds.result_msd)
        print(numpy.shape(msds.result_msd))
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
    def __init__(self,txyz,box=None,plot_trajectory=False,mode='simu'):#,account='tplab'
        R"""
        parameters:
            TXYZ:array[frames,particles,xyz]
            BOX: [lx,ly,lz,xy,xz,yz]
            account: 'remote' or 'tplab'
            mode: 'simu' or 'exp'. 'simu' to select particles in box; 'exp' to direct compute msd
            chips_max: the largest number of chips the frames can be cut into.
        """
        self.txyz = txyz
            
        if mode == 'simu':
            self.box = box
            self.select_stable_trajectories_simu(plot_trajectory)
        elif mode == 'exp':
            self.txyz_stable = txyz
        self.frames,self.particles,self.dimensions=self.txyz.shape
        self.chips_max = self.frames-1
        #self.account = account
        self.mode = mode

    def select_stable_trajectories_simu(self):
        R"""
        introduction:
            transform trajectory data from simulation with periodic boundary condition 
            into trajectories of which particles never move across boundary(box).
        parameters:
            positions: (N_frames,N_particles,3)

        exp:
            import points_analysis_2D as pa
            msd_class = pa.msd(gsd_data.txyz,gsd_data.box,account='remote')#pa.msd(txyz,box,account='remote')
            msd_class.compute_t_chips()
            msd_class.plot()

        """
        #dedrift?
        self.frames,self.particles,self.dimensions=self.txyz.shape
        if hasattr(self,'box'):
            #print(locals())#local variable not of class
            self.dtxyz = self.txyz[1:,:,:] - self.txyz[:self.frames-1,:,:]
            #cross is true
            list_crossleft = self.dtxyz[:,:,0] > 0.9*self.box[0]
            list_crossbottom = self.dtxyz[:,:,1] > 0.9*self.box[1]
            list_crossright = self.dtxyz[:,:,0] < -0.9*self.box[0]
            list_crosstop = self.dtxyz[:,:,1] < -0.9*self.box[1]
            #mark all the frames where cross event occur as True
            list_crossx = numpy.logical_or(list_crossleft,list_crossright)
            list_crossy = numpy.logical_or(list_crossbottom,list_crosstop)
            list_cross = numpy.logical_or(list_crossx,list_crossy)
            #mark all the particles who have experienced cross event as True
            list_cross_true = list_cross[0,:]
            #list_cross_true = list_cross_true[0]#remove empty extra dimension
            #print(list_cross_true.shape)
            i=0
            while i<self.particles:
                list_cross_true[i] = numpy.max(list_cross[:,i])
                i = i + 1
            list_stable_id = numpy.where(list_cross_true[:]==False)
            list_stable_id = list_stable_id[0]#remove empty extra dimension
            #print(list_stable_id.shape)
            
            self.txyz_stable = self.txyz[:,list_stable_id,:]
            self.particles = list_stable_id.shape[0]

    def compute_msd_normal(self,interval_max=0.8):
        R"""
        EXP:
            import points_analysis_2D as pa
            import numpy as np
            path_to_results = '/home/remote/xiaotian_file/data/20221129/video_5'
            txyz_npy_filename = path_to_results+'/'+'txyz_stable.npy'
            txyz_stable = np.load(txyz_npy_filename)
            msds = pa.dynamic_points_analysis_2d(txyz_stable,mode='exp')
            msds.compute_msd_normal(0.95)
            time_log_file = path_to_results+'/'+'DefaultVideo_5.txt'
            time_log = np.loadtxt(time_log_file)
            png_filename=path_to_results+'/'+'msd_normal_long_loglog.png'
            msds.plot_msd(time_log,png_filename)
        """
        dt_max = int(interval_max*self.chips_max)#to ensure the robustness of statistics
        #k_max = 0
        dt_start = 1
        self.record_msd=numpy.zeros((dt_max-dt_start+1,2))#[interval m, msd_m]
        
        list_dt=numpy.arange(dt_start,dt_max+1)
        self.record_msd[:,0]=list_dt
        """
        print(self.txyz_stable[0,:5,:2])
        print(self.txyz_stable[0,:5,:2]-1.0)
        """
        """
        ff = 1
        dt_xyz = self.txyz_stable[dt_start:dt_max+1,:,:]
        while ff < dt_max+1:
            dt_xyz[ff-1] = dt_xyz[ff-1] - self.txyz_stable[0]
            ff = ff+1
        """
        #print(dt_xyz[1,1])
        dt_xyz = self.txyz_stable[dt_start:dt_max+1,:,:] - self.txyz_stable[0,:,:]
        #print(dt_xyz_test[1,1])
        dt_xyz2 = numpy.square(dt_xyz)
        sum_dt_r2 = numpy.sum(dt_xyz2,axis=2)
        sum_dt_r2 = numpy.sum(sum_dt_r2,axis=1)
        self.record_msd[:,1]=sum_dt_r2[:]/self.particles

    def compute_atmsd_scan_t(self,interval_max=0.66):
        R"""
        method:
            scan t axis, 1 frame for each interval. here may exist some overlaps between two intervals who share the same length of time interval m.
        record_msd:[interval m, msd_m]
        dt_start: the frame to start scan 
        dt: time interval
        dt_max: max time interval
        """
        dt_max = int(interval_max*self.chips_max)#to ensure the robustness of statistics
        #k_max = 0
        dt_start = 1
        self.record_msd=numpy.zeros((dt_max-dt_start+1,2))#[interval m, msd_m]
        
        list_dt=numpy.arange(dt_start,dt_max+1)
        self.record_msd[:,0]=list_dt
        for dt in list_dt:
            dt_xyz = self.txyz_stable[dt:,:,:] - self.txyz_stable[:-dt,:,:]
            chips = dt_xyz.shape[0]
            dt_xyz2 = numpy.square(dt_xyz)
            sum_dt_r2 = numpy.sum(dt_xyz2)
            self.record_msd[dt-dt_start,1]=sum_dt_r2/chips/self.particles
            
    def compute_atmsd_t_chips(self,interval_max=0.1,msd_per_particle=False):
        R"""
        parameters:
            interval_max: (0,1) to which the longest interval dt to which msd is computed.
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
        dt_max = int(interval_max*self.chips_max)#to ensure the robustness of statistics
        dt_start = 1
        self.record_msd=numpy.zeros((dt_max-dt_start+1,2))#[interval m, msd_m]
        #you can select points [log 1, 1/100*log m_max,log m_max]
        list_dt=numpy.arange(dt_start,dt_max+1)
        self.record_msd[:,0]=list_dt
        
        for dt in list_dt:
            chips = int((self.chips_max)/dt)# how many chips the time axis is divided into
            list_frames = dt*numpy.arange(0,chips+1)
            dt_xyz = self.txyz_stable[list_frames,:,:]
            self.dt_xyz = dt_xyz[1:,:,:] - dt_xyz[:-1,:,:]
            sq_dm_xyz2 = numpy.square(self.dt_xyz)
            sum_sq_dm_r2 = numpy.sum(sq_dm_xyz2)
            self.record_msd[dt-dt_start,1]=sum_sq_dm_r2/chips/self.particles

            if msd_per_particle:
                id_sum_sq_dm_r2_chips = numpy.sum(sq_dm_xyz2,axis=2)
                print(id_sum_sq_dm_r2_chips.shape)
                id_sum_sq_dm_r2 = numpy.sum(id_sum_sq_dm_r2_chips,axis=0)
                print(id_sum_sq_dm_r2.shape)
                if 'record_msd_id' in locals():
                    record_msd_id[dt-dt_start,:]=id_sum_sq_dm_r2/self.dt_xyz.shape[0]
                else:
                    record_msd_id=numpy.zeros((dt_max-dt_start+1,self.particles))
                    record_msd_id[dt-dt_start,:]=id_sum_sq_dm_r2/self.dt_xyz.shape[0]

        if msd_per_particle:
            return record_msd_id

    def compute_nearest_neighbor_displacements(self,unit='sigma'):
        R"""
        parameters:
            input: 
                trajectory_stable (unit = sigma or um)
                'stable' means all the particles in the 
                field of view from the 1st frame to the end.
                if the trajectory is not stable, a particle id 
                in different frames would not represent the same particle.

            return: 
                ts_id_dxy:
                    a pandas.dataframe contains ["particle", "x","y"."frame"]
                    which contains all the nn displacements over time.
                #dict_c_nb:
                #    a dict() contains {'particle id': [list id of neighbors]} at given frame.

                    #shape = [N_frames, N_ids, dxy], dtype = float64
                    #shape = [N_frames, N_ids, if_nb_change], dtype = bool
                    
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
        id_dxy_pd_columns = ["particle","neighbors", "x","y"]
        import pandas as pd
        #get frame-wise 1st neighbor bond length
        #list_1st_bond_length_framewise
        self.frames,self.particles,self.dimensions=self.txyz.shape
        list_frames = range(self.frames)
        list_1st_bond_length_framewise = numpy.zeros((self.frames,))
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
            id_dxy,dict_c_nb = result.get_nearest_neighbor_dispalcement()
            id_dxy_pd = pd.DataFrame(id_dxy)
            #print(id_dxy_pd.head())
            if frame == 0:
                #ts_id_dxy = t_id_dxy
                id_dxy_pd.columns = id_dxy_pd_columns
                id_dxy_pd['frame'] = frame
                ts_id_dxy = id_dxy_pd
                #print(ts_id_dxy.head())
                #print(ts_id_dxy.tail())
            else:
                #ts_id_dxy = numpy.concatenate((ts_id_dxy,t_id_dxy), axis=0) #[frame,id,dx,dy]
                id_dxy_pd.columns = id_dxy_pd_columns#[id,dx,dy,frame]
                id_dxy_pd['frame'] = frame
                ts_id_dxy = pd.concat([ts_id_dxy,id_dxy_pd])
                #print(ts_id_dxy.tail())
                #print(ts_id_dxy.shape)
        #select stable trajectories
        """
        print(ts_id_dxy[:,1].max)
        list_ids = numpy.unique(ts_id_dxy[:,1])
        id_check = numpy.zeros((list_ids.shape[0],),dtype=bool)
        for id in list_ids:
            ts_dxy = ts_id_dxy[ts_id_dxy[:,1]==id]

            if ts_dxy.shape[0] == self.frames:
                id_check[id] = True
        ts_id_dxy_stable = ts_id_dxy[:,id_check,:,:]
        """
        average_1st_bond_length = numpy.average(list_1st_bond_length_framewise) 
        #particle_tracking.select_stable_trajectory()is great!
        pd.DataFrame.to_csv(ts_id_dxy,'ts_id_dxy.csv')
        return ts_id_dxy,average_1st_bond_length

    def compute_nearest_neighbor_displacements_stable(self,ts_id_dxy,unit='sigma'):#compute_nearest_neighbor_displacements_stable
        R"""
        input:
            ts_id_dxy:
                a pandas.dataframe contains ["particle", "x","y"."frame"]
                which contains all the nn displacements over time.
        return:
            x txy_ids_stable:
                it is (Nframes,Nparticles,xyz) 3-dimension array,
                which contains ['particle'] which are always in the field of view in video.

            x msd_c:
                a pandas.dataframe contains ["particle", "x","y"."frame"]
        """
        pt = particle_track()
        txy_ids_stable = pt.select_stable_trajectory(tpxyz=ts_id_dxy)
        """
        list_center_id= ts_id_dxy['particle'].values
        frame_max = ts_id_dxy['frame'].values.max()
        for id in list_center_id:
           ts_dxy = ts_id_dxy[ts_id_dxy['particle'] == id]
           ts_dxy = numpy.array(ts_dxy.values)
           print(ts_dxy.shape)
        """
        
            #self.compute_atmsd_t_chips()
        return txy_ids_stable

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
        list_1st_bond_length_framewise = numpy.zeros((self.frames,))
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
            t_id_dxy = numpy.zeros((id_dxy.shape[0],id_dxy.shape[1]+1))  
            t_id_dxy[:,0] = frame
            t_id_dxy[:,1:] = id_dxy  
            if frame == 0:
                ts_id_dxy = t_id_dxy
            else:
                ts_id_dxy = numpy.concatenate((ts_id_dxy,t_id_dxy), axis=0) #[frame,id,dx,dy]
            #print(ts_id_dxy.shape)
        #select stable trajectories
        print(ts_id_dxy[:,1].max)
        list_ids = numpy.unique(ts_id_dxy[:,1])
        id_check = numpy.zeros((list_ids.shape[0],),dtype=bool)
        for id in list_ids:
            ts_dxy = ts_id_dxy[ts_id_dxy[:,1]==id]

            if ts_dxy.shape[0] == self.frames:
                id_check[id] = True
        ts_id_dxy_stable = ts_id_dxy[:,id_check,:,:]
        #particle_tracking.select_stable_trajectory()is great!

    def plot_msd(self,time_log=None,png_filename='msd_loglog.png',um_sec=True,sigma=False,lindemann=False):
        R"""
        introduction:
            input: 'txyz_stable.csv' 
            output: msd plot
        """
        import matplotlib.pyplot as plt 
        plt.figure()
        if time_log is None:
            plt.loglog(self.record_msd[:,0],self.record_msd[:,1])
        else:
            nth_time = self.record_msd.shape[0]+1#start from 0
            time_msd=time_log[1:nth_time]
            plt.loglog(time_msd,self.record_msd[:,1])

        if self.mode =='simu':
            plt.title("Mean Squared Displacement")
            plt.xlabel("$t$ (steps)" )
            plt.ylabel("MSD$(t)$ (sigma^2)")
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
        list_dm_frames = numpy.arange(self.dm_xyz.shape[0])
        for i in list_dm_frames:
            #numpy.dot()
            self.dm_xyz_tr=numpy.transpose(self.dm_xyz)
            numpy.matmul(self.dm_xyz),self.dm_xyz_tr)
        """
        
        
        """
        m=0#frame interval
        record_disp2_1 = numpy.zeros((self.sp[0]-1,3))
        while m<self.sp[0]-1:#m < N frame
            k=0#frame scanning
            record_disp2_m = numpy.zeros((1,3))
            
            while k<self.sp[0]-1-m:
                disp_m = positions[k+m]-positions[k]#199-0 what is the problem?
                disp2_m = numpy.dot(numpy.transpose(disp_m),disp_m) 

                if disp2_m[0,0]>0:
                    print(m)

                record_disp2_m = record_disp2_m + [disp2_m[0,0],disp2_m[1,1],disp2_m[2,2]]
                k=k+1

            record_disp2_m = record_disp2_m/(self.sp[0]-1-m)#mean of all the chips with k frames interval
            record_disp2_1[m] = record_disp2_m[0]
            
            m=m+1
        self.result_msd_xyz = record_disp2_1/self.sp[1]#divided by N_particles
        self.result_msd = numpy.zeros((self.sp[0]-1,1))
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

    def plot_trajectory(self,png_filename='trajectory_stable.png'):
        list_stable_id = range(self.particles)#txyz_stable.shape[1]
        plt.figure()
        for index_particle in list_stable_id:
            txyz_ith = self.txyz_stable[:,index_particle,:]
            plt.plot(txyz_ith[:,0],txyz_ith[:,1])

        unit = '(um)' #'(sigma)'
        plt.xlabel("$x$ "+unit )
        plt.ylabel("$y$ "+unit )
        #png_filename = '/home/'+self.account+'/Downloads/'+'traj_stable.png'
        plt.savefig(png_filename)
        plt.close() 

    def plot_trajectory_single_particle(self,prefix='/home/',unit = '(um)'):
        R"""
        EXP:
            path_to_results = '/home/remote/xiaotian_file/data/20221129/video_5'
            txyz_npy_filename = path_to_results+'/'+'txyz_stable.npy'
            txyz_stable = np.load(txyz_npy_filename)
            msds = pa.dynamic_points_analysis_2d(txyz_stable,mode='exp')
            msds.plot_trajectory_single_particle(path_to_results+'/')
        """
        list_stable_id = range(self.particles)#txyz_stable.shape[1]
        
        for particle_id in list_stable_id:
            txyz_ith = self.txyz_stable[:,particle_id,:]
            plt.figure()
            plt.plot(txyz_ith[:,0],txyz_ith[:,1])
            plt.xlabel("$x$ "+unit )#'(sigma)'
            plt.ylabel("$y$ "+unit )
            png_filename = 'traj_stable_'+str(int(particle_id))+'.png'
            plt.savefig(prefix+png_filename)
            plt.close()    
