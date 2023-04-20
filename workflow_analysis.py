import numpy as np
import matplotlib.pyplot as plt
import numpy
import os

import scipy

class get_msd_from_gsd:
    def __init__(self):
        pass
    def get_msd(self,simu_index=5208,seed=9,account='remote'):
        import points_analysis_2D as pa
        self.prefix = "/home/"+account+"/Downloads/"
        self.str_index = str(int(simu_index))+'_'+str(int(seed))
        self.create_folder()
        self.prefix = self.prefix+self.str_index+"/"
        gsd_data = pa.proceed_gsd_file(account=account,simu_index=simu_index,seed=seed)
        
        gsd_data.get_trajectory_data(self.prefix)
        gsd_data.get_trajectory_stable_data(self.prefix)

        file_txyz_npy = self.prefix+'txyz_stable.npy'#_stable
        txyz_stable = np.load(file_txyz_npy)
        """
        import points_analysis_2D as pa
        dpa = pa.dynamic_points_analysis_2d(txyz_stable)
        dpa.plot_trajectory_single_particle(10)
        """
        msdm = pa.mean_square_displacement(txyz_stable,'simu')
        record=msdm.compute_atmsd_scan_t_log_record()
        png_filename = self.prefix+'msd_chips_long_loglog.png'
        msdm.plot_msd_t_chips(png_filename=png_filename)

        txt_filename = self.prefix+'msd_chips_long_loglog.txt'
        np.savetxt(txt_filename,record)
        return self.prefix

    def create_folder(self):
        folder_name=self.prefix+self.str_index#+"/"
        #check if the folder exists
        isExists=os.path.exists(folder_name)
        if isExists:
            pass
        else:
            os.makedirs(folder_name)
    
class get_displacement_field:
    def __init__(self):
        pass

    def get_displacement_field_1():
        import points_analysis_2D as pa
        #gsd_data = pa.proceed_gsd_file(account='remote',simu_index=5208,seed=9)
        
        save_prefix = "/home/remote/Downloads/5208_9/"
        file_txyz_npy = save_prefix+'txyz.npy'#_stable
        
        file_txyz_npy = save_prefix+'txyz_stable.npy'#_stable
        txyz_stable = np.load(file_txyz_npy)
        txyz_stable = txyz_stable[:,:,0:2]
        #dpa = pa.dynamic_points_analysis_2d(txyz_stable)
        #dpa.plot_trajectory_single_particle(save_prefix)
        df = pa.displacemnt_field_2D(txyz_stable)
        png_filename=save_prefix+'displacement_field.png'
        uv = df.get_displacement_field_xy(frame_index_start=0,plot=True,png_filename=png_filename)
        vector_avg = np.average(uv,0)
        vector_avg_dt = vector_avg/20000
        txyz_stable_tuned = np.zeros(np.shape(txyz_stable))
        for i in range(20001):
            txyz_stable_tuned[i,:] = txyz_stable[i,:]-vector_avg_dt*i
        png_filename=save_prefix+'displacement_field_tuned.png'
        df.get_displacement_field_xy(plot=True,png_filename=png_filename)
        filename_npy = save_prefix + 'txyz_stable_tuned'
        np.save(filename_npy,txyz_stable_tuned)
    
    def get_displacement_field_normal(self,save_prefix):
        import points_analysis_2D as pa
        #gsd_data = pa.proceed_gsd_file(account='remote',simu_index=5208,seed=9)
        
        #save_prefix = "/home/remote/Downloads/5208_9/"
        file_txyz_npy = save_prefix+'txyz.npy'#_stable
        
        file_txyz_npy = save_prefix+'txyz_stable.npy'#_stable
        txyz_stable = np.load(file_txyz_npy)
        txyz_stable = txyz_stable[:,:,0:2]
        #dpa = pa.dynamic_points_analysis_2d(txyz_stable)
        #dpa.plot_trajectory_single_particle(save_prefix)
        df = pa.displacemnt_field_2D(txyz_stable)
        png_filename=save_prefix+'displacement_field.png'
        uv = df.get_displacement_field_xy(plot=True,png_filename=png_filename)

    def get_displacement_field_dedrift():
        import points_analysis_2D as pa
        #gsd_data = pa.proceed_gsd_file(account='remote',simu_index=5208,seed=9)
        
        save_prefix = "/home/remote/Downloads/5208_9/"
        file_txyz_npy = save_prefix+'txyz.npy'#_stable
        
        file_txyz_npy = save_prefix+'txyz_stable.npy'#_stable
        txyz_stable = np.load(file_txyz_npy)
        txyz_stable = txyz_stable[:,:,0:2]

        np.shape(txyz_stable)
        txyz_stable_tuned = np.zeros(np.shape(txyz_stable))
        for i in range(20001):
            txyz_stable_tuned[i,:] = txyz_stable[i,:]#-vector_avg_dt*i
        #dpa = pa.dynamic_points_analysis_2d(txyz_stable)
        #dpa.plot_trajectory_single_particle(save_prefix)
        df = pa.displacemnt_field_2D(txyz_stable)
        png_filename=save_prefix+'displacement_field.png'
        uv = df.get_displacement_field_xy(frame_index_start=0,plot=True,png_filename=png_filename)
        vector_avg = np.average(uv,0)
        vector_avg_dt = vector_avg/20000
        
        png_filename=save_prefix+'displacement_field_tuned.png'
        df.get_displacement_field_xy(plot=True,png_filename=png_filename)
        filename_npy = save_prefix + 'txyz_stable_tuned'
        np.save(filename_npy,txyz_stable_tuned)

class check_box:
    def __init__(self):
        import points_analysis_2D as pa
        gsd_data = pa.proceed_gsd_file(account='remote',simu_index=5208,seed=9)
        
        lx = gsd_data.box[0]
        ly = gsd_data.box[1]
        fence = [[-lx/2,ly/2],[lx/2,ly/2],[lx/2,-ly/2],[-lx/2,-ly/2],[-lx/2,ly/2]]
        fence = np.array(fence)
        for i in range(10):
            pos = gsd_data.read_a_frame(i)
            posp = pos[:]+[0,gsd_data.box[1]]
            posn = pos[:]-[0,gsd_data.box[1]]
            posr = pos[:]+[gsd_data.box[0],0]
            posl = pos[:]-[gsd_data.box[0],0]
            pos_all = np.concatenate([pos,posp,posn,posr,posl],axis=0)
            fig,ax = plt.subplots()
            #ax.scatter(pos[:,0],pos[:,1],c='k')
            ax.scatter(pos_all[:,0],pos_all[:,1],c='k')
            ax.plot(fence[:,0],fence[:,1])
            ax.plot()
            ax.set_aspect('equal','box')

            plt.show()
    #def square_fence(self):

class merge_data:
    R"""
    import workflow_analysis as wa
    for seed in range(9):

        mg = wa.merge_data()
        d3 = mg.read_txt(5208,seed)
        d2 = mg.read_txt(5431,seed)
        d1 = mg.read_txt(5430,seed)
        d12 = mg.merge_txts(d1,d2)
        d123 = mg.merge_txts(d12,d3)
        png_filename = 'msd_chips_long_loglog_'+str(int(seed))+'.png'
        mg.plot_msd_t_chips(d123,png_filename)
    """
    def __init__(self):
        pass

    def read_txt(self,simu_index,seed):
        str_index = str(int(simu_index))+'_'+str(int(seed))
        prefix = "/home/remote/Downloads/"
        folder = str_index+"/"
        txt_filename = "msd_chips_long_loglog.txt"
        filename = prefix+folder+txt_filename
        record = np.loadtxt(filename)
        return record
    
    def merge_txts(self,record1,record2):
        #np.unique(record1,axis=0)
        record = np.concatenate((record1,record2))
        return record
    
    def plot_msd_t_chips(self,record_msd,png_filename='msd_chips_long_loglog.png'):
        R"""
        introduction:
            input: 'txyz_stable.csv' 
            output: msd plot
        """
        import matplotlib.pyplot as plt 
        plt.figure()

        plt.loglog(record_msd[:,0],record_msd[:,1])


        plt.title("Mean Squared Displacement")
        plt.xlabel("$t$ (steps)" )
        plt.ylabel("MSD$(t)(\sigma^2)$ ")
        
        plt.savefig(png_filename)
        plt.close()
    
    def generate_power_law_mix(self):
        R"""
        introduction: 
            ln(y) = k*ln(x), linear in loglog,
            when -k > 1 , the cumulation is convergent series
        return: 
            x,y
        """
        x0 = np.linspace(0,2,10)
        x = np.power(10,x0)
        y0 = np.linspace(-5,-1,10)
        y = np.power(10,y0)
        #xy = np.concatenate((x,y),axis=1)
        xy = np.zeros((10,2))
        xy[:,0]=x
        xy[:,1]=y

        x0 = np.linspace(3,5,10)
        x2 = np.power(10,x0)
        y0 = np.linspace(-0.5,0.5,10)
        y2 = np.power(10,y0)
        xx = np.concatenate((x,x2))
        yy = np.concatenate((y,y2))
        xy = np.zeros((20,2))
        xy[:,0]=xx
        xy[:,1]=yy

        x0 = np.linspace(5,6,10)
        x3 = np.power(10,x0)
        y0 = np.linspace(0.5,1.8,10)
        y3 = np.power(10,y0)
        xx = np.concatenate((x,x2,x3))
        yy = np.concatenate((y,y2,y3))
        xy = np.zeros((30,2))
        xy[:,0]=xx
        xy[:,1]=yy
        return xy
    
class optimize_polyfit:
    R"""
    intro:
        to fit the index of power law
    exp:
        import workflow_analysis as wa
        wa.optimize_polyfit()
    """
    def __init__(self):
        prefix = "/home/remote/Downloads/4302_9/reference_occupation/waiting_time_split_2_halves/"
        txt_filename = prefix+"list_waiting_time_1_bin.txt"
        xy = np.loadtxt(txt_filename)
        nums = 20
        lnxy = np.log10(xy[:nums])
        lnx = lnxy[:,0]
        lny = lnxy[:,1]
        z = np.polyfit(lnx,lny,1)
        k_str = str(np.around(z[0],3))
        print(z)

        fig,ax = plt.subplots()
        self.plot(ax,xy[:,0],xy[:,1])
        lny_s = z[0]*lnx+z[1]
        #y_s = np.exp(lny_s)
        y_s = np.power(10,lny_s)
        self.plot_dash(ax,xy[:nums,0],y_s)
        ax.annotate('k='+k_str,[10,100],c='orange')
        png_filename = prefix+"list_waiting_time_1_polyfit.png"
        plt.savefig(png_filename)
        #plt.show()
        #scipy.optimize.curve_fit()
        return z
    
    def plot_dash(self,ax,x,y):
        ax.loglog(x,y,dashes=[6, 2],c='orange')
    
    def plot(self,ax,x,y):
        ax.loglog(x,y)

R"""
Idea: the existence and advantages(speed and ratio); mechanism and universality of transformation
Show a single experiment/simulation:
Show the transition from hex to honeycomb;√
Show the class of phase transition( diffusionless transformation)(direct pin + interstitial filling)
Show the kinetics during the transition using the evolution of order parameters (psi3, the fraction of neighbor-changing events, smooth?);√
Show the kinetics during the transition using the evolution of configuration;？
Show the constrained diffusion mechanism(chain activation), and clips of displacement field or trajectory. arrow_circle color bar: the rank of nb change events among ids with long displacements. 4302_9_89-106(4chain√) 110-113(2chain) 1193-1195(diffusion√) 1361-1365(chain√)
Show the constrained diffusion mechanism(Wating time)(Apart from instanton, are vacancies in honeycomb transition activations?)method1: displacement & instanton;method2:vacancies located at traps or interstitials are activations, similar to spin-ups in the east model.
Show the constrained diffusion mechanism(dynamic heterogeneity)√
Show the pin mechanism√

Show a list of simulation:
Diagram of psi3 tuned by lcr and k;√
Show the difference between different conditions:
Direct trap VS dual lattice;√
Higher transformation ratio; faster transformation dynamics.
"""
class show_transition_from_hex_to_honeycomb:
    def __init__(self):
        import data_analysis_cycle as dac
        daw = dac.data_analysis_workflow()
        daw.get_bond_plot()
    
    def get_bond_plot(self,directory,data_name=None,trap_filename=None,trap_lcr=None,io_only=False):
        R"""
        input:
            directory from self.gsd_to_txyz
            trap_filename:
                '/home/remote/hoomd-examples_0/testhoneycomb3-8-12'
                '/home/remote/hoomd-examples_0/testhoneycomb3-8-12-part1'
                '/home/remote/hoomd-examples_0/testkagome3-11-6'
                '/home/remote/hoomd-examples_0/testkagome_part3-11-6'
            io_only: just return results, not proceeding data.
        return:
            a series of figures with particles(mark neighbor changes), bonds, traps
        example:
            import data_analysis_cycle as da
            get_traj = da.data_analysis()
            directory,data_name = get_traj.gsd_to_txyz('remote',4448,9,io_only=True)
            get_traj.txyz_to_bond_plot(directory,data_name,
                trap_filename='/home/remote/hoomd-examples_0/testkagome_part3-11-6',trap_lcr=0.89,
                    io_only=True)
        """
        #write a routine class
        import pandas as pd
        import points_analysis_2D as pa
        file_txyz_stable = directory + 'txyz_stable.npy'
        txyz_stable = numpy.load(file_txyz_stable)
        dpa = pa.dynamic_points_analysis_2d(txyz_stable,mode='simu')
        #particle id should be set as what in txyz_stable!
        bond_cut_off = 6
        if not io_only:
            dpa.compute_nearest_neighbor_displacements(csv_prefix=directory,bond_cut_off=bond_cut_off)
            file_ts_id_dxy = directory + 'ts_id_dxy.csv'
            ts_id_dxy = pd.read_csv(file_ts_id_dxy)
            if_nb_change_int,n_particle_nb_stable = dpa.monitor_neighbor_change_event(ts_id_dxy=ts_id_dxy,csv_prefix=directory)
            dpa.get_hist_neighbor_change_event(if_nb_change_int,n_particle_nb_stable,directory)
        count_nb_change_event_rate = numpy.load(directory+'count_nb_change_event_rate.npy')
        dpa.plot_hist_neighbor_change_event(count_nb_change_event_rate,directory)
        """
        if_nb_change_int, n_particle_nb_stable, png_filename==dpa.monitor_neighbor_change_event(ts_id_dxy=ts_id_dxy,csv_prefix=directory)
        dpa.plot_hist_neighbor_change_event(if_nb_change_int, n_particle_nb_stable, png_filename=)
        """
        if not data_name is None:
            file_list_sum_id_nb_stable = directory + 'list_sum_id_nb_stable.csv'
            list_sum_id_nb_stable = pd.read_csv(file_list_sum_id_nb_stable)
            #dpa.plot_bond_neighbor_change(nb_change=list_sum_id_nb_stable,data_name=data_name,prefix=directory,bond_cut_off=bond_cut_off,
            #                                    show_traps=True,trap_filename='/home/remote/hoomd-examples_0/testhoneycomb3-8-12',trap_lcr=0.79)
            dpa.plot_bond_neighbor_change_oop(data_name=data_name,prefix=directory,nb_change=list_sum_id_nb_stable,bond_cut_off=bond_cut_off,
                                               trap_filename=trap_filename,trap_lcr=trap_lcr)
            """
            dpa.plot_bond_neighbor_change_oop()
            dpa.draw_bonds.draw_bonds_conditional_bond()
            dpa.draw_bonds.plot_neighbor_change(txyz_stable,nb_change)
            dpa.draw_bonds.plot_traps(trap_filename,LinearCompressionRatio)
            """

class show_tuned_image:
    def __init__(self):
        pass

    def read_image(self):
        import points_analysis_2D as pa
        import particle_tracking as pt
        prefix = '/home/tplab/Downloads/20230321/' 
        image_filename = prefix+'DefaultImage_12.jpg' 
        spe = pt.particle_track()
        spe.single_frame_particle_tracking(image_filename,D=11,minmass=800)#,calibration=True
        print(spe.xy)
        spe.xy[:,1] = -spe.xy[:,1]
        pixel2um = 3/32
        points_um = spe.xy*pixel2um
        spa = pa.static_points_analysis_2d(points=points_um)
        hist_filename = image_filename+'hist.jpg'
        spa.get_first_minima_bond_length_distribution(lattice_constant=1.0,hist_cutoff=4.6,png_filename=hist_filename,x_unit='um')
        check = [2,spa.bond_first_minima_left]
        bond_filename = image_filename+'bond.jpg'
        line = spa.draw_bonds_conditional_bond_oop(check=check,png_filename=bond_filename)
    
    def draw_tuned_image(self):
        import points_analysis_2D as pa
        import particle_tracking as pt
        prefix = '/home/tplab/Downloads/20230321/' 
        image_filename = prefix+'DefaultImage_12.jpg' 
        spe = pt.particle_track()
        spe.single_frame_particle_tracking(image_filename,D=11,minmass=800)#,calibration=True
        pixel2um = 3/32
        um2pxiel = 1/pixel2um
        spa = pa.static_points_analysis_2d(points=spe.xy)
        hist_filename = image_filename+'hist_pix.jpg'
        spa.get_first_minima_bond_length_distribution(lattice_constant=1.0,hist_cutoff=4.6*um2pxiel,png_filename=hist_filename,x_unit='pixel')
        check = [2*um2pxiel,44.16]#44.16#spa.bond_first_minima_left
        bond_filename = image_filename+'bond_pix.jpg'
        f0 = plt.imread(image_filename)
        line = spa.draw_bonds_conditional_bond_for_image_oop(f0,check=check,png_filename=bond_filename,
                            x_unit='(pix)',axis_limit=[400,800,700,250])#,axis_limit=[400,800,250,700]

    def draw_points_with_conditional_bond(self,xy,bond_length=None,bond_length_limmit=[0.9,2.0]):
        R"""
        Introduction:
            copy from points_analysis_2D.bond_plot_module.draw_points_with_conditional_bond()
        Parameters:
            xy: particle positions of a frame, with no one removed.
            bond_length: [particle_id1,particle_id2, bond_length] for txyz.
            bond_length_limmit: limit the shortest and longest bond( in bond_length) to draw.
        weight of shapes:
            bond(blue line) < particles(black circle) < neighbor_change(orange circle) < traps(red cross)
            0   1   2   3   
        Examples:
        """
        if not (bond_length is None):
            bond_check= tuple([bond_length_limmit[0],bond_length_limmit[1]])
            #add lines for edges
            for i in range(np.shape(bond_length)[0]):
                if (bond_length[i,2] > bond_check[0])&(bond_length[i,2] < bond_check[1]) :
                    edge = tuple(bond_length[i,0:2].astype(int))
                    pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
                    line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b',zorder=0)#,lineStyle='dashed'
                    #self.ax.add_line(line)
        return line

class show_waiting_time_brownian:
    R"""
    intro:
        when a neighbor-changing event happened on the particle i, 
        search the first neighbor-changing event happening on the particle i,
        and record the waiting time step. count the number of events for each time step. 
    parameters:
        csv_file: [frame, particle_id, sum_id_neighbors, if_nb_change]
    example:
        import workflow_analysis as wa
        wtb = wa.show_waiting_time_brownian()
        wtb.plot_hist()
    """
    def __init__(self) :
        pass

    def compute(self):
        import pandas as pd
        prefix = '/home/remote/Downloads/4302_9/'
        filename_csv = prefix+'list_sum_id_nb_stable.csv'
        list_sum_id_nb_stable = pd.read_csv(filename_csv)
        data_nb_changed = list_sum_id_nb_stable[list_sum_id_nb_stable['if_nb_change'] == True]
        ids = data_nb_changed['particle_id'].values
        ids = np.unique(ids)
        ids_waiting_time = []
        ids_waiting_time = np.array(ids_waiting_time)
        for id in ids:
            id_data_nb_changed = data_nb_changed[data_nb_changed['particle_id']==id]
            id_frames = id_data_nb_changed['frame'].values
            id_frames_sorted = np.sort(id_frames)
            id_waiting_time = id_frames_sorted[1:] - id_frames_sorted[:-1]
            ids_waiting_time = np.concatenate((ids_waiting_time,id_waiting_time))
        filename = prefix +'list_ids_waiting_time.txt'
        np.savetxt(filename,ids_waiting_time)
        #pd.DataFrame.to_csv()
    
    def plot_hist(self):
        prefix = '/home/remote/Downloads/4302_9/brownian/'
        filename = prefix +'list_ids_waiting_time.txt'
        ids_waiting_time = np.loadtxt(filename)

        fig,ax = plt.subplots()
        count_bins = ax.hist (ids_waiting_time,bins=500)#,log=True,bins=20,range=[0,100]
        #plt.show()
        _count=count_bins[0]
        _bins=count_bins[1]
        fig2,ax2 = plt.subplots()
        ax2.semilogy(_bins[1:],_count)#semilogy,loglog
        ax2.set_xlabel('waiting time brownian(k steps)')
        ax2.set_ylabel('count (1)')
        png_filename = prefix + 'list_waiting_time_1_bin_brownian.png'
        plt.savefig(png_filename)
        txt_filename = prefix + 'list_waiting_time_1_bin_brownian.txt'
        cb = np.zeros((_count.size,2))
        cb[:,0] = _bins[1:]
        cb[:,1] = _count
        np.savetxt(txt_filename,cb)

class show_waiting_time_dynamical_facilitation:
    R"""
    intro:
        when a neighbor-changing event happened on the particle i, 
        search the first neighbor-changing event happening nearby(within rcut),
        and record the waiting time step. count the number of events for each time step. 
    parameters:
        csv_file: [frame, particle_id, sum_id_neighbors, if_nb_change]
                (particle_ids are of txyz_stable, and edge_cut to nb_stable )
        rcut = a*lcr*sqrt(3)*(1+10%). where a is lattice constant, lcr is linear compression ratio, 
            10% for thermal fluctuation. 
    example:
        import workflow_analysis as wa
        wtdf = wa.show_waiting_time_dynamical_facilitation()
        rcut=3*0.81*1.73*1.1#a*lcr*sqrt(3)*(1+10%)
        #fn = wtdf.compute(rcut)
        prefix = '/home/remote/Downloads/4302_9/dynamical_facilitation_nb/'
        fn=prefix +'list_waiting_time_dyfa_rcut4.62.txt'
        wtdf.plot_hist(fn)
    """
    def __init__(self) :
        pass
    def compute(self,rcut):
        R"""
        idea:
            base on data_nb_changed, add 'x','y' cloumns with [0,0].
            search txyz_stable[frame,id] and fill x,y cloumns.
            calculate distance framewise 
            frame+1, calculate dr, select row which is minima, 
                if minima < rcut, record and finish;
                if not, continue.
            list [frame,particle_id,x,y]
        """
        import pandas as pd
        prefix = '/home/remote/Downloads/4302_9/'
        filename_txyz_stable = prefix+'txyz_stable.npy'
        txyz_stable = np.load(filename_txyz_stable)
        txyz_stable = txyz_stable[:,:,:2]
        filename_csv = prefix+'list_sum_id_nb_stable.csv'
        list_sum_id_nb_stable = pd.read_csv(filename_csv)
        data_nb_changed = list_sum_id_nb_stable[list_sum_id_nb_stable['if_nb_change'] == True]
        frames = data_nb_changed['frame'].values
        frames = np.unique(frames)
        frames = np.sort(frames)
        #initialize data-recorder
        list_waiting_time = []
        list_waiting_time = np.array(list_waiting_time)

        for frame in frames:#small to large?
            frame_data_nb_changed = data_nb_changed[data_nb_changed['frame']==frame]
            ids = frame_data_nb_changed['particle_id'].values
            ids = np.sort(ids)
            for id in ids:
                #search following events happening nearby framewisely
                for frame_follow in frames:
                    if frame_follow > frame:
                        frame_follow_data_nb_changed = data_nb_changed[data_nb_changed['frame']==frame_follow]
                        ids_follow = frame_follow_data_nb_changed['particle_id'].values
                        ids_follow = np.sort(ids_follow)
                        list_dxy = txyz_stable[frame_follow,ids_follow] - txyz_stable[frame,id]#before or now?
                        list_dxy2 = np.square(list_dxy)
                        list_dr2 = list_dxy2[:,0]+list_dxy2[:,1]
                        list_dr2 = list_dr2 - rcut*rcut
                        list_dr2 = np.sort(list_dr2)
                        if list_dr2[0]<=0:
                            dframe = frame_follow-frame
                            dframe = np.array([dframe])
                            list_waiting_time = np.concatenate((list_waiting_time,dframe))
                            #[frame1,id1,frame2,id2] to link event chain?
                            #what if one event causes two or more events?
                            break
        filename = prefix +'dynamical_facilitation_nb/'+'list_waiting_time_dyfa_rcut'+str(np.around(rcut,2))+'.txt'
        np.savetxt(filename,list_waiting_time)
        return filename
        
    def plot_hist(self,filename):
        #dynamical_facilitation_based_on_neighber_change
        prefix = '/home/remote/Downloads/4302_9/dynamical_facilitation_nb/'
        ids_waiting_time = np.loadtxt(filename)

        fig,ax = plt.subplots()
        count_bins = ax.hist (ids_waiting_time,bins=500,log=True)#,log=True,bins=20,range=[0,100]
        #plt.show()
        _count=count_bins[0]
        _bins=count_bins[1]
        fig2,ax2 = plt.subplots()
        ax2.loglog(_bins[1:],_count)#semilogy,loglog
        ax2.set_xlabel('waiting time(k steps)')
        ax2.set_ylabel('count (1)')
        png_filename = prefix + 'list_waiting_time_1_bin_nb.png'
        plt.savefig(png_filename)
        txt_filename = prefix + 'list_waiting_time_1_bin_nb.txt'
        cb = np.zeros((_count.size,2))
        cb[:,0] = _bins[1:]
        cb[:,1] = _count
        np.savetxt(txt_filename,cb)

class show_waiting_time_interstitial_motion:
    R"""
    intro:
        when an interstitial motion(or large jump) happened on the particle i, 
        search the first interstitial motion happening nearby(within rcut),
        and record the waiting time step. count the number of events for each time step. 

        to get coarse grained trajectory
        title: Excitations Are Localized and Relaxation Is Hierarchical in Glass-Forming Liquids
        DOI: 10.1103/PhysRevX.1.021013
    parameters:
        csv_file: [frame, particle_id, sum_id_neighbors, if_nb_change]
                (particle_ids are of txyz_stable, and edge_cut to nb_stable )
    """
    def __init__(self):
        pass

    def trajectory_single_particle(self,txyz_stable_id=16):
        import points_analysis_2D as pa
        import numpy as np
        prefix = '/home/remote/Downloads/4302_9/'
        filename_txyz_stable = prefix+'txyz_stable.npy'
        txyz_stable = np.load(filename_txyz_stable)
        dpa = pa.dynamic_points_analysis_2d(txyz_stable)
        dpa.plot_trajectory_single_particle(txyz_stable_id,prefix)

    def trajectory_coarse_grain_single_particle(self,txyz_stable_id=16):
        R"""
        input:
            txyz_stable_id: [frame,x,y] for a single particle
        output:
            method1: compare the i-1 frame an the i+1 frame, merge trajectory in the same region.
            for the same particle id, dr(df) = r(frame+df)-r(frame),
            method2: average 3 frames of positions and draw the trajectory
        """
        import points_analysis_2D as pa
        prefix = '/home/remote/Downloads/4302_9/'
        filename_txyz_stable = prefix+'txyz_stable.npy'
        txyz_stable = np.load(filename_txyz_stable)
        sz = numpy.shape(txyz_stable)
        df = pa.dynamical_facilitation_module()
        id = 16
        dr = df.scan_displacement_t(txyz_stable[:,id,:2],)

        pass
    
    def trajectory_coarse_grain(self,txyz_stable_id):
        R"""
        input:
            txyz_stable_id: [frame,x,y] for a single particle
        output:
            method1: compare the i-1 frame an the i+1 frame, merge trajectory in the same region.
            for the same particle id, dr(df) = r(frame+df)-r(frame),
            method2: average 3 frames of positions and draw the trajectory
        """
        import points_analysis_2D as pa
        prefix = '/home/remote/Downloads/4302_9/'
        filename_txyz_stable = prefix+'txyz_stable.npy'
        txyz_stable = np.load(filename_txyz_stable)
        sz = numpy.shape(txyz_stable)
        df = pa.dynamical_facilitation_module()
        for id in range(sz[1]):
            dr = df.scan_displacement_t(txyz_stable[:,id,:2],)
        pass

class show_polygon_dye:
    def __init__(self):
        pass

    def get_points(self):
        import points_analysis_2D as pa
        import numpy as np
        
        prefix = '/home/remote/Downloads/4302_9/'
        prefix_write='/home/remote/Downloads/4302_9/polygon_6/'
        self.read_data(prefix_write)

    def get_points_plot(self):
        #self.test_patch()
        import points_analysis_2D as pa
        import numpy as np
        
        prefix = '/home/remote/Downloads/4302_9/'
        prefix_write='/home/remote/Downloads/4302_9/polygon_6/'#delauny_voronoi_polygon/'
        filename_txyz = prefix+'txyz.npy'
        txyz = np.load(filename_txyz)
        nframes =2000
        #self.read_data(prefix_write)
        record_polygon_s_rate = np.zeros((20,nframes))
        for frame in range(nframes):
            if frame>1421:
                points = txyz[frame]
                self.p2d = pa.static_points_analysis_2d(points)
                self.get_first_minima_ridge_length_distribution(prefix_write,frame=frame)#,io_only=True
                #draw bonds selected
                #self.p2d.draw_polygon_patch_oop()
                count_polygon_relative = self.draw_bonds_conditional_ridge_oop(prefix_write,frame)#,io_only=True
                record_polygon_s_rate[:,frame] = count_polygon_relative[:,1]
        fname=prefix_write+'polygon_n_vs_frame_rate_as_value2000.txt'
        np.savetxt(fname,record_polygon_s_rate)
        #self.read_data(prefix_write)

    def get_points_plot_xylim(self):
        #self.test_patch()
        import points_analysis_2D as pa
        import numpy as np
        
        prefix = '/home/remote/Downloads/4302_9/'
        prefix_write='/home/remote/Downloads/4302_9/polygon_6_lim/'#delauny_voronoi_polygon/'
        filename_txyz = prefix+'txyz.npy'
        txyz = np.load(filename_txyz)
        nframes =32#2000
        #self.read_data(prefix_write)
        for frame in range(nframes):
            if frame<31:
                points = txyz[frame]
                self.p2d = pa.static_points_analysis_2d(points)
                self.get_first_minima_ridge_length_distribution(prefix_write,frame=frame,io_only=True)
                #draw bonds selected
                #self.p2d.draw_polygon_patch_oop()
                self.draw_bonds_conditional_ridge_oop(prefix_write,frame,limit=True)
                
    def get_point_simulations(self,simu_index=4302,seed=9,account='remote'):
        import points_analysis_2D as pa
        import numpy as np
        str_simu_index = str(int(simu_index))+'_'+str(seed)
        prefix = '/home/'+account+'/Downloads/'+str_simu_index#+'/'
        #check if the folder exists
        isExists=os.path.exists(prefix)
        if not isExists:
            os.makedirs(prefix)
        prefix = prefix+'/'
        prefix_write=prefix+'polygon_6'#delauny_voronoi_polygon/'
        isExists=os.path.exists(prefix_write)
        if not isExists:
            os.makedirs(prefix_write)
        prefix_write = prefix_write+'/'
        filename_txyz = prefix+'txyz.npy'
        txyz = np.load(filename_txyz)

        nframes =50
        #self.read_data(prefix_write)
        record_polygon_s_rate = np.zeros((20,nframes))
        for frame in range(nframes):
            if frame>=0:
                points = txyz[frame]
                self.p2d = pa.static_points_analysis_2d(points)
                self.get_first_minima_ridge_length_distribution(prefix_write,frame=frame)#,io_only=True
                #draw bonds selected
                #self.p2d.draw_polygon_patch_oop()
                count_polygon_relative = self.draw_bonds_conditional_ridge_oop(prefix_write,frame)#,io_only=True
                record_polygon_s_rate[:,frame] = count_polygon_relative[:,1]
        fname=prefix_write+'polygon_n_vs_frame_rate_as_value2000.txt'
        np.savetxt(fname,record_polygon_s_rate)
        #self.read_data(prefix_write)

    def plot_points(self,points,prefix_write):
        #self.plot_points(points,prefix_write)
        fig,ax = plt.subplots()
        #ax.scatter(pos[:,0],pos[:,1],c='k')
        ax.scatter(points[:,0],points[:,1],c='k')
        ax.set_aspect('equal','box')
        png_filename = 'points.png'
        plt.savefig(prefix_write+png_filename)

    def plot_bond_points(self,prefix_write):
       
        self.bond_length = self.p2d.bond_length
        self.points = self.p2d.points
        fig,ax = plt.subplots()
        plt.scatter(self.points[:,0],self.points[:,1],color='k')
        #add lines for edges
        for edge in self.bond_length[:,0:2].astype(int):
            #print(edge)
            pt1,pt2 = [self.points[edge[0]],self.points[edge[1]]]
            #plt.gca().add_line(plt.Line2D(pt1,pt2))
            line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='b')
            plt.gca().add_line(line)
        
        ax.set_aspect('equal','box')
        png_filename = 'bonds.png'
        plt.savefig(prefix_write+png_filename)
    
    def plot_bond_hist(self,prefix_write):
        self.bond_length = self.p2d.bond_length
        fig,ax = plt.subplots()
        hist_cutoff=5
        ax.hist(self.bond_length[:,2],bins=20,range=[0,hist_cutoff])
        png_filename = 'bonds_length_hist.png'
        plt.savefig(prefix_write+png_filename)

    def plot_ridge_hist(self,prefix_write):
        self.ridge_length = self.p2d.voronoi.ridge_length
        fig,ax = plt.subplots()
        hist_cutoff=5
        ax.hist(self.ridge_length[:],bins=20,range=[0,hist_cutoff])#
        png_filename = 'ridge_length_hist.png'
        plt.savefig(prefix_write+png_filename)

    def plot_ridge_hist_log(self,prefix_write):
        self.ridge_length = self.p2d.voronoi.ridge_length
        fig,ax = plt.subplots()
        hist_cutoff=5
        count_bins = ax.hist(self.ridge_length[:],bins=20,range=[0,hist_cutoff])#,bins=20,range=[0,hist_cutoff]
        _count=count_bins[0]
        _bins=count_bins[1]
        fig2,ax2 = plt.subplots()
        ax2.semilogx(_bins[1:],_count)#semilogy,loglog
        #ax2.set_xlabel('waiting time (k steps)')
        #ax2.set_ylabel('count (1)')
        png_filename = 'ridge_length_hist_log.png'
        plt.savefig(prefix_write+png_filename)

    def plot_voronoi(self,prefix_write):
        from scipy.spatial import voronoi_plot_2d
        fig = voronoi_plot_2d(self.p2d.voronoi)
        png_filename = 'voronoi.png'
        plt.savefig(prefix_write+png_filename)
    
    def get_first_minima_ridge_length_distribution(self,prefix_write,hist_cutoff=5,frame=None,io_only=True):
        print(frame)
        if io_only:
            self.p2d.get_first_minima_ridge_length_distribution(hist_cutoff=hist_cutoff)
        else :
            png_filename = prefix_write+"ridge_length_hist"+str(int(frame))+".png"
            self.p2d.get_first_minima_ridge_length_distribution(hist_cutoff=hist_cutoff,png_filename=png_filename)
    def draw_bonds_conditional_ridge_oop(self,prefix_write,frame,io_only=False,limit=False):
        count_polygon_relative = self.p2d.get_conditional_bonds_and_simplices()
        if not io_only:
            png_filename = prefix_write+"bond_vertices_patch"+str(int(frame))+".png"
            fig,ax = plt.subplots()
            fig,ax = self.p2d.draw_polygon_patch_oop(ax)
            if limit:
                ax.set_xlim(-4,12)
                ax.set_ylim(-5,16)
            self.p2d.draw_bonds_simplex_conditional_oop(png_filename=png_filename,x_unit='(1)',fig=fig,ax=ax)#
        #plt.show()
        return count_polygon_relative
    
    def vertices_vs_simplices(self,i):
        R"""
        results:
            the ith vertex whose related simplex is also the ith delauney triangle. that is,
            voronoi.vertices[i]~list_ridge_points == delaunay.simplices[i].
        parameters:
            voronoi.vertices: N_vertices rows of (x,y)
            voronoi.ridge_points: n rows of (start_point_row, end_point_row).
                Indices of the points between which each Voronoi ridge lies
            bond_length: n rows of (start_point_row, end_point_row, bond_length)
            voronoi.ridge_vertices: n rows of (start_vertex_row,end_vertex_row)
            voronoi.ridge_length: n rows of (length)
            delaunay: scipy.spatial.Delaunay
            
            cutted_bonds: n rows of (start_point_row, end_point_row)
            delaunay.simplices: n rows of (point_row_1, point_row_2, point_row_3) vertices of triangles
            --------------------------------
            linked_triangles: n rows of (simplices_row_1,simplices_row_2,simplices_row_3,simplices_row_4,etc)
        example:
            spd = wa.show_polygon_dye()
            spd.get_points()
            list1=[123,201,255]
            for i in list1:
                spd.data_introduction(i)#200
        """
        #relate simplices and vertices
        print('index of vertex')
        print(i)
        #print('position of vertices')
        #print(self.p2d.voronoi.vertices[i])
        list_ridge_index1 = np.where(self.p2d.voronoi.ridge_vertices[:,:]==i)
        #print('index of ridges')
        #print(self.p2d.voronoi.ridge_vertices[list_ridge_index1[0]])
        print('index of points')
        list_points = self.p2d.voronoi.ridge_points[list_ridge_index1[0]]
        list_points = np.unique(list_points)
        print(list_points)
        #print('positions of points of the simplices')
        #print(self.p2d.voronoi.points[list_points])
        list_points_from_simp = self.p2d.delaunay.simplices[i]
        list_points_from_simp = np.sort(list_points_from_simp)
        print(list_points_from_simp)

    def test_patch(self):
        import numpy as np
        from matplotlib.patches import Circle, Wedge, Polygon
        from matplotlib.collections import PatchCollection
        import matplotlib.pyplot as plt

        # Fixing random state for reproducibility
        np.random.seed(19680801)


        fig, ax = plt.subplots()

        resolution = 50  # the number of vertices
        N = 3
        x = np.random.rand(N)
        y = np.random.rand(N)
        radii = 0.1*np.random.rand(N)
        patches = []
        for x1, y1, r in zip(x, y, radii):
            circle = Circle((x1, y1), r)
            patches.append(circle)

        x = np.random.rand(N)
        y = np.random.rand(N)
        radii = 0.1*np.random.rand(N)
        theta1 = 360.0*np.random.rand(N)
        theta2 = 360.0*np.random.rand(N)
        for x1, y1, r, t1, t2 in zip(x, y, radii, theta1, theta2):
            wedge = Wedge((x1, y1), r, t1, t2)
            patches.append(wedge)

        # Some limiting conditions on Wedge
        patches += [
            Wedge((.3, .7), .1, 0, 360),             # Full circle
            Wedge((.7, .8), .2, 0, 360, width=0.05),  # Full ring
            Wedge((.8, .3), .2, 0, 45),              # Full sector
            Wedge((.8, .3), .2, 45, 90, width=0.10),  # Ring sector
        ]

        for i in range(N):
            polygon = Polygon(np.random.rand(N, 2), closed=True)
            patches.append(polygon)

        colors = 100 * np.random.rand(len(patches))
        p = PatchCollection(patches, alpha=0.4)
        p.set_array(colors)
        ax.add_collection(p)
        fig.colorbar(p, ax=ax)

        plt.show()
    
    def read_data(self,prefix_write):
        fname=prefix_write+'polygon_n_vs_frame_rate_as_value2000.txt'
        record_polygon_s_rate = np.loadtxt(fname)
        fig,ax = plt.subplots()
        ax.plot(record_polygon_s_rate[4])
        ax.set_title('polygon_n_vs_frame_rate_as_value')
        ax.set_xlabel('time(k step)')
        ax.set_ylabel('honeycomb(%)')
        
        png_filename = prefix_write+"polygon_n_vs_frame_rate.png"
        plt.savefig(png_filename)

class show_disp_field:
    def __init__(self):
        pass
    
    def get_points_plot(self):
        #self.test_patch()
        import points_analysis_2D as pa
        import numpy as np
        
        prefix = '/home/remote/Downloads/4302_9/'
        prefix_write='/home/remote/Downloads/4302_9/'
        filename_txyz = prefix+'txyz_stable.npy'
        txyz = np.load(filename_txyz)
        png_filename=prefix_write+'displacement_field_xy_0_8_part.png'
        
        self.p2d = pa.dynamic_points_analysis_2d(txyz)
        self.p2d.displacement_field_module()
        self.p2d.displacemnt_field.get_displacement_field_xy(0,8,plot=True,png_filename=png_filename,limit=True)
        


