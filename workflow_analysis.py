import numpy as np
import matplotlib.pyplot as plt
import numpy
import os

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
        fig.savefig(png_filename)
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

class show_bond_image:
    R"""
    import workflow_analysis as wa
    prefix_image='image_to_proceed/'
    sti = wa.show_tuned_image(prefix_image,'DefaultImage_2.jpg')
    sti.draw_tuned_image()
    prefix = '/home/tplab/xiaotian_file/'
    """
    def __init__(self,prefix,filename):
        self.prefix = prefix
        self.filename = filename

    def read_image(self):
        import points_analysis_2D as pa
        import particle_tracking as pt
        prefix = self.prefix#'/home/tplab/Downloads/20230321/' 
        image_filename = prefix+self.filename#'DefaultImage_12.jpg' #'-' can not be recognized by plt.imread()!
        spe = pt.particle_track()
        spe.single_frame_particle_tracking(image_filename,D=11,minmass=800,calibration=True)#
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
        R"""
        image_name,trap,D,minmass,hist_ctoff,bond_length,axis_limit
        20230321-IMAGE12,honey_part,11,800,4.6*32/3,44.16,[400,800,700,250]
        20230113-video8-2246,kagome_part,11,400,90,60,[100,900,900,100]
        
        """
        import points_analysis_2D as pa
        import particle_tracking as pt
        prefix = self.prefix#'/home/tplab/Downloads/20230321/' 
        image_filename = prefix+self.filename#'DefaultImage_12.jpg' 
        spe = pt.particle_track()
        spe.single_frame_particle_tracking(image_filename,D=11,minmass=400,axis_limit=[100,901,100,901])#,calibration=True
        pixel2um = 3/32
        um2pxiel = 1/pixel2um
        spa = pa.static_points_analysis_2d(points=spe.xy)
        hist_filename = image_filename+'hist_pix.jpg'
        spa.get_first_minima_bond_length_distribution(lattice_constant=1.0,hist_cutoff=120,png_filename=hist_filename,x_unit='pixel')
        check = [2*um2pxiel,spa.bond_first_minima_left]#44.16#spa.bond_first_minima_left
        bond_filename = image_filename+'bond_pix.jpg'
        #f0 = plt.imread(image_filename)
        line = spa.draw_bonds_conditional_bond_for_image_oop(spe.image,check=check,png_filename=bond_filename,
                            x_unit='(pix)')#,axis_limit=[400,800,250,700]

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

    def trajectory_coarse_grain_single_particle(self,txyz_stable_id=16,n_frame_to_coarse=2):
        R"""
        input:
            txyz_stable_id: [frame,x,y] for a single particle
        output:
            method1: compare the i-1 frame an the i+1 frame, merge trajectory in the same region.
            for the same particle id, dr(df) = r(frame+df)-r(frame),
            method2: average 3 frames of positions and draw the trajectory
        results:
            n_frame_to_coarse=1,plot: lines + clusters
            n_frame_to_coarse=2,plot: lines + clusters
        """
        import points_analysis_2D as pa
        prefix = '/home/remote/Downloads/4302_9/'
        filename_txyz_stable = prefix+'txyz_stable.npy'
        txyz_stable = np.load(filename_txyz_stable)
        sz = numpy.shape(txyz_stable)
        df = pa.dynamical_facilitation_module()
        id = txyz_stable_id
        dr = df.scan_displacement_t(txyz_stable[:,id,:2])
        self.plot_trajectory(txyz_stable[:,id,:2])
        self.plot_disp_hist(dr)
        #self.plot_trajectory(txyz_stable[:,id,:2])
        plt.show()

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
    
    def plot_trajectory(self,txy):
        fig,ax = plt.subplots()
        #ax.scatter(pos[:,0],pos[:,1],c='k')
        ax.plot(txy[:,0],txy[:,1],c='k')
        ax.set_aspect('equal','box')
        
        #png_filename = 'points.png'
        #plt.savefig(prefix_write+png_filename)
    
    def plot_disp_hist(self,disp):
        fig,ax = plt.subplots()
        #ax.scatter(pos[:,0],pos[:,1],c='k')
        ax.hist(disp,bins=50,log=True)
        ax.set_ylabel('count(1)')
        ax.set_xlabel('$dr(\sigma)$')
        #ax.set_aspect('equal','box')
        

class show_polygon_dye:
    def __init__(self):
        pass
    
    def plot_polygon_bond_xylim(self):
        """
        example:
            import workflow_analysis as wa
            spd = wa.show_polygon_dye() 
            spd.plot_polygon_bond_xylim()
        """
        import points_analysis_2D as pa
        account = 'tplab'#'remote'
        prefix = '/home/'+account+'/Downloads/4302_9/'
        prefix_write='/home/'+account+'/Downloads/4302_9/polygon_6/'#_lim,delauny_voronoi_polygon/'
        filename_txyz = prefix+'txyz.npy'
        txyz = np.load(filename_txyz)
        nframes =32#2000
        #self.read_data(prefix_write)
        for frame in [8,104,1427]:#range(nframes):
            #if frame<6:
            points = txyz[frame]
            self.p2d = pa.static_points_analysis_2d(points)#,hide_figure=False
            self.get_first_minima_ridge_length_distribution(prefix_write,frame=frame,io_only=True)
            #draw bonds selected
            #self.p2d.draw_polygon_patch_oop()
            self.draw_bonds_conditional_ridge_oop(prefix_write,frame,limit=True)

    def get_points_plot(self):
        import points_analysis_2D as pa
        import numpy as np
        import proceed_file as pf
        pfile = pf.proceed_file_shell()
        prefix_simu, str_simu_index = pfile.create_prefix_results(simu_index=4302,seed=9)
        prefix_write = pfile.create_folder(prefix_simu,'polygon_6')
        print(prefix_write)
        #delauny_voronoi_polygon/'
        filename_txyz = prefix_simu+'txyz.npy'
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
                
    def get_points_plot_loop(self,seed):
        import points_analysis_2D as pa
        import numpy as np
        import proceed_file as pf
        pfile = pf.proceed_file_shell()
        prefix_simu, str_simu_index = pfile.create_prefix_results(account='tplab',simu_index=4302,seed=seed)
        prefix_write = pfile.create_folder(prefix_simu,'polygon_6')
        print(prefix_write)
        #delauny_voronoi_polygon/'
        filename_txyz = prefix_simu+'txyz.npy'
        txyz = np.load(filename_txyz)
        nframes =40
        #self.read_data(prefix_write)
        for frame in [8,104,1427]:#range(nframes):
            #if frame>=20:
            self.p2d = pa.static_points_analysis_2d(txyz[frame])
            self.get_first_minima_ridge_length_distribution(prefix_write,frame=frame,io_only=True)#
            #draw bonds selected
            #self.p2d.draw_polygon_patch_oop()
            #self.draw_bonds_conditional_ridge_oop(prefix_write,frame,limit=True)#,io_only=True
            png_filename = prefix_write+"bond_vertices_patch"+str(int(frame))+".png"
            fig,ax = plt.subplots()
            self.p2d.draw_polygon_patch_oop(fig,ax)#fig,ax = 
            a=10
            lim = [[-a,a],[-a,a]]#-4,12,-5,16
            self.p2d.draw_bonds_simplex_conditional_oop(png_filename=png_filename,x_unit='($\sigma$)',axis_limit=lim,fig=fig,ax=ax)

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
    
    def get_first_minima_ridge_length_distribution(self,prefix_write,hist_cutoff=5,frame=None,io_only=False):
        print(frame)
        if io_only:
            self.p2d.get_first_minima_ridge_length_distribution(hist_cutoff=hist_cutoff)
        else :
            png_filename = prefix_write+"ridge_length_hist"+str(int(frame))+".png"
            self.p2d.get_first_minima_ridge_length_distribution(hist_cutoff=hist_cutoff,png_filename=png_filename)
    
    def draw_bonds_conditional_ridge_oop(self,prefix_write,frame,io_only=False,limit=False):
        count_polygon_relative = self.p2d.get_conditional_bonds_and_simplices()
        if not io_only:
            png_filename = prefix_write+"bond_vertices_patch"+str(int(frame))+".pdf"
            fig,ax = plt.subplots()
            fig,ax = self.p2d.draw_polygon_patch_oop(ax)
            if limit:
                ax.set_xlim([-18,18])#-4,12
                ax.set_ylim([-18,18])#-5,16
            self.p2d.draw_bonds_simplex_conditional_oop(png_filename=png_filename,x_unit='($\sigma$)',fig=fig,ax=ax)#   
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
        fig.savefig(png_filename)

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
    
    def get_disp(self,seed,frame1,frame2):
        import points_analysis_2D as pa
        import numpy as np
        import proceed_file as pf
        pfile = pf.proceed_file_shell()
        prefix_simu, str_simu_index = pfile.create_prefix_results(simu_index=4302,seed=seed)
        prefix_write = pfile.create_folder(prefix_simu,'polygon_6')
        filename_txyz = prefix_simu+'txyz_stable.npy'
        txyz = np.load(filename_txyz)
        png_filename=prefix_write+'displacement_field_xy_'+str(int(frame1))+'_'+str(int(frame2))+'.png'#'_part.png'

        self.p2d = pa.dynamic_points_analysis_2d(txyz)
        self.p2d.displacement_field_module()
        self.p2d.displacemnt_field.get_displacement_field_xy(frame1,frame2,plot=True,png_filename=png_filename)#,limit=True
    
    def get_disp_lim(self,seed,frame):
        import points_analysis_2D as pa
        import numpy as np
        import proceed_file as pf
        pfile = pf.proceed_file_shell()
        prefix_simu, str_simu_index = pfile.create_prefix_results(simu_index=4302,seed=seed)
        prefix_write = pfile.create_folder(prefix_simu,'polygon_6')
        filename_txyz = prefix_simu+'txyz_stable.npy'
        txyz = np.load(filename_txyz)
        png_filename=prefix_write+'displacement_field_xy_0_'+str(int(frame))+'_part.png'

        self.p2d = pa.dynamic_points_analysis_2d(txyz)
        self.p2d.displacement_field_module()
        self.p2d.displacemnt_field.get_displacement_field_xy(0,frame,plot=True,png_filename=png_filename,limit=[[-4,8],[-3,16]])
                
    def get_bicolor_disp(self,seed,frame1,frame2,lim):
        R"""
        import workflow_analysis as wa
        sdf = wa.show_disp_field()
        seed=[0,1,2,8,9]
        frame=[12,10,29,6,8]
        lim=[[[-5,9],[-5,10]],[[11,21],[-13,1]],[[9,21],[-7,7]],[[-9,4],[-12,3]],[[-4,8],[-3,16]]]
        for i in range(5):
            print(seed[i])
            sdf.get_bicolor_disp(seed[i],frame[i],lim[i])
        """
        #seed=9
        lcr=0.81
        trap_filename='/home/remote/hoomd-examples_0/testhoneycomb3-8-12-part1'
        traps = numpy.loadtxt(trap_filename)*lcr
        import proceed_file as pf
        pfile = pf.proceed_file_shell()
        prefix_simu, str_simu_index = pfile.create_prefix_results(simu_index=4302,seed=seed)
        prefix_write = pfile.create_folder(prefix_simu,'pin_check')
        file_txyz_stable = prefix_simu + 'txyz_stable.npy'
        txyz_stable = numpy.load(file_txyz_stable)
        import points_analysis_2D as pa
        df = pa.dynamical_facilitation_module()
        #df.get_pin_bool(traps,txyz_stable,prefix_write,1.0)
        #plot
        #frame=8
        file_t_pin_bool = '/home/remote/Downloads/'+str_simu_index+'/pin_check/t_pin_bool.npy'
        t_pin_bool = numpy.load(file_t_pin_bool)
        
        self.p2d = pa.dynamic_points_analysis_2d(txyz_stable)
        self.p2d.displacement_field_module()
        #lim=[[-25,25],[-20,20]]
        png_filename=prefix_write+'displacement_field_xy_'+str(int(frame1))+'_'+str(int(frame2))+'.png'#_part
        self.p2d.displacemnt_field.get_bicolor_disp(t_pin_bool[frame2],frame1,frame2,plot=True,png_filename=png_filename,limit=lim,traps=traps)
        #png_filename=prefix_write+'displacement_field_xy_0_2000.png'#'_part.png'
        #self.p2d.displacemnt_field.get_bicolor_disp(t_pin_bool[2000],0,2000,plot=True,png_filename=png_filename)

class show_dual_lattice:
    R"""
    import getDataAndScatter
    getDataAndScatter.get_dual_lattice()
    """
    def __init__(self) -> None:
        pass

    def go(self):
        trap_filename='/home/remote/hoomd-examples_0/testhoneycomb3-8-12'
        traps=np.loadtxt(trap_filename)
        LinearCompressionRatio=1
        traps=np.multiply(traps,LinearCompressionRatio)
        fig,ax = plt.subplots()
        #ax.scatter(pos[:,0],pos[:,1],c='k')
        ax.scatter(traps[:,0], 
                traps[:,1],
                c='r')#,marker = 'x'
        ax.set_aspect('equal','box')
        """
        limit=[[],[]]
        ax.set_xlim(limit[0])
        ax.set_ylim(limit[1])
        """
            
class show_cairo_order_parameter:
    def __init__(self) -> None:
        pass
    def compute_theoretical_diagram(self):
        R"""
        output:
            record: [cn3,cn4,order_parameter_for_cairo]
        example:
            import workflow_analysis as wa
            ca = wa.show_cairo_order_parameter()
            xyc = ca.compute()
            ca.plot_diagram(xyc)
        """
        num_x = 99
        list_num = np.linspace(0.01,0.99,num_x)
        num_scan = 0.5*num_x*(num_x+1)
        record = np.zeros((int(num_scan),3))
        count = 0
        for cn3 in list_num:
            for cn4 in list_num:
                p1 = cn3+cn4
                if p1<=1:
                    cn4_normal = cn4*2
                    r34 = cn3/cn4_normal
                    r43 = cn4_normal/cn3
                    p2 = np.minimum(r34,r43)
                    p = p1*p2
                    record[count] = [cn3,cn4,p]
                    count = count + 1
        return record
    
    def compute_cairo_order_parameter(self,cn3,cn4):
        R"""
        output:
            p: order_parameter_for_cairo
        """
        p1 = cn3+cn4
        if p1<=1:
            cn4_normal = cn4*2
            r0 = cn3*cn4
            if r0>0:
                r34 = cn3/cn4_normal
                r43 = cn4_normal/cn3
                p2 = np.minimum(r34,r43)
                p = p1*p2
            else:
                p = 0
        else:
            print('error: cn3+cn4>1')
            p = -1
                    
        return p

    def plot_diagram(self,xyc,account='remote'):
        fig,ax = plt.subplots()
        #ax.scatter(pos[:,0],pos[:,1],c='k')
        scat = ax.scatter(xyc[:,0],xyc[:,1],c=xyc[:,2])
        ax.set_aspect('equal','box')
        ax.set_xlabel('cn3(1)')
        ax.set_ylabel('cn4(1)')
        ax.set_title('order_parameter_for_cairo')
        fig.colorbar(scat,ax=ax)
        prefix = "/home/"+account+"/Downloads/"
        png_filename=prefix+'test_order_parameter_for_cairo.png'
        plt.savefig(png_filename)
    
    def plot_depin_order(self):
        prefix = '/home/tplab/Downloads/record_results/cairo/depin_from_Cairo/Pcairo_egct2lcra/'
        filename = '2573-2582klcn346'#'2563-2572klcn346'#
        import opertateOnMysql as osql
        osql.loadDataToMysql(prefix+filename,'depin_from_cairo_egct2lcra')
        #data = np.loadtxt(prefix+filename)

class show_cn_k:
    def __init__(self) -> None:
        pass
    
    def get_cn_k(self,a_frame,lattice_constant,gsd_data,time_steps,i):
        R"""
        CN0 % should be 0 for all the particles must be linked by bond.
        CN1 % is likely to be edge?
        CN2 % in body(edge-cutted) shows the mechanical unstability
        CN3 % shows the proportion of honeycomb.
        CN4 % shows the proportion of kagome.
        CN6 % shows the proportion of hexagonal.
        CN5/7 % shows the proportion of disclination.
        """
        #print('index '+str(i))
        #print(snap.particles.position[137])
        a_frame.get_coordination_number_conditional(lattice_constant=lattice_constant)#cut edge to remove CN012
        ccn = a_frame.count_coordination_ratio#[time_steps,psi3,psi6]
        ccn = numpy.transpose(ccn)
        if not "record_cn" in locals():#check if the variable exists
            #load CN_k s
            record_cn = numpy.zeros((gsd_data.num_of_frames,numpy.shape(ccn)[1]+1))
            record_cn[:,0] = time_steps#range(10)##gsd frame is different from log frame for period set 100 vs 2e3
        #print(numpy.shape(ccn)[1])
        record_cn[i,1:numpy.shape(ccn)[1]+1] = ccn
    
    def show_cn_k(self,frame_cut,record_cn,prefix,str_index):
        plt.figure()
        if frame_cut == 0:#frame_cut is set to abstract a part of the process to watch in detail
            #plt.plot(record_cn[:,0],record_cn[:,1],label='CN_0')
            #plt.plot(record_cn[:,0],record_cn[:,2],label='CN_1')
            #plt.plot(record_cn[:,0],record_cn[:,3],label='CN_2')
            plt.plot(record_cn[:,0],record_cn[:,4],label='CN_3')
            plt.plot(record_cn[:,0],record_cn[:,5],label='CN_4')
            plt.plot(record_cn[:,0],record_cn[:,6],label='CN_5')
            plt.plot(record_cn[:,0],record_cn[:,7],label='CN_6')
            plt.plot(record_cn[:,0],record_cn[:,8],label='CN_7')
            #plt.plot(record_cn[:,0],record_cn[:,9],label='CN_8')
            #plt.plot(record_cn[:,0],record_cn[:,-1],label='CN_9')
            png_filename = prefix +'T_VS_CN_k'+'index'+str_index+'egcut'+'.png'
        else:
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,1],label='CN_0')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,2],label='CN_1')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,3],label='CN_2')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,4],label='CN_3')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,5],label='CN_4')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,6],label='CN_5')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,7],label='CN_6')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,8],label='CN_7')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,9],label='CN_8')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,-1],label='CN_9')
            png_filename = prefix +'T_VS_CN_k_tcut'+'index'+str_index+'egcut'+'.png'
        plt.legend()
        plt.title('CN_k '+'index:'+str_index)
        plt.xlabel('time(steps)')
        plt.ylabel('CN_k(1)')
        #plt.show()
        plt.savefig(png_filename)
        record_filename = prefix +'T_VS_CN_k_cut'+'index'+str_index+'.txt'
        numpy.savetxt(record_filename,record_cn)
        plt.close()

    def reorganize_cn_k_txt_and_show(self):
        filename1 = '/home/remote/Downloads/record_results/pin_hex_to_honeycomb_2m/'+'T_VS_CN_k_cutindex5238_9.txt'
        record_cn1 = np.loadtxt(filename1)
        cn3h = record_cn1[:,4]#4634'CN_3'
        t = record_cn1[:,0]#steps

        filename2 = '/home/remote/Downloads/4302_9/cn_k/'+'T_VS_CN_k_cutindex4302_9.txt'
        record_cn2 = np.loadtxt(filename2)
        cn3hp = record_cn2[:,4]#4302'CN_3'

        fig,ax = plt.subplots()
        ax.plot(t,cn3h,label='honeycomb traps')
        ax.plot(t,cn3hp,label='honeycomb_part traps')
        ax.legend()
        ax.set_title('honeycomb vs honeycomb_part')
        ax.set_xlabel('time(steps)')
        ax.set_ylabel('$CN_3(1)$')
        #plt.show()
        prefix_write = '/home/remote/Downloads/4302_9/cn_k/'
        png_filename = prefix_write+'T_VS_CN_kindex_5238_4302_egcut.png'
        fig.savefig(png_filename)
        #fig.savefig(png_filename)

    def show_list_of_cn_k(self,prefix,simu_index,seed):
        prefix = '/home/remote/Downloads/record_results/pin_hex_to_honeycomb_2m/'
        simu_index = 4302

        fig,ax = plt.subplots()

        for seed in range(10):
            cn3h,t = self.extract_cn_k_txt(prefix,simu_index,seed)
            
            ax.plot(t,cn3h,label='honeycomb traps')


        ax.legend()
        ax.set_title('honeycomb vs honeycomb_part')
        ax.set_xlabel('time(steps)')
        ax.set_ylabel('$CN_3(1)$')
        #plt.show()
        prefix_write = '/home/remote/Downloads/4302_9/cn_k/'
        png_filename = prefix_write+'T_VS_CN_kindex_5238_4302_egcut.png'
        fig.savefig(png_filename)

        
    
    def extract_cn_k_txt(self,prefix,simu_index,seed):
        str_index = str(int(simu_index))+'_'+str(int(seed))
        filename = prefix+'T_VS_CN_k'+'index'+str_index+'egcut'+'.png'
        #filename = '/home/remote/Downloads/record_results/pin_hex_to_honeycomb_2m/'+'T_VS_CN_k_cutindex5238_9.txt'
        record_cn = np.loadtxt(filename)
        cn3h = record_cn[:,4]#4634'CN_3'
        t = record_cn[:,0]#steps
        return cn3h,t

    def show_cn_k_compare(self,cni,cnj,time_steps,record_cn,prefix,str_index,frame_cut=0):
        plt.figure()
        if frame_cut == 0:#frame_cut is set to abstract a part of the process to watch in detail
            #plt.plot(record_cn[:,0],record_cn[:,1],label='CN_0')
            #plt.plot(record_cn[:,0],record_cn[:,2],label='CN_1')
            #plt.plot(record_cn[:,0],record_cn[:,3],label='CN_2')
            plt.plot(record_cn[:,0],record_cn[:,4],label='CN_3')
            plt.plot(record_cn[:,0],record_cn[:,5],label='CN_4')
            plt.plot(record_cn[:,0],record_cn[:,6],label='CN_5')
            plt.plot(record_cn[:,0],record_cn[:,7],label='CN_6')
            plt.plot(record_cn[:,0],record_cn[:,8],label='CN_7')
            #plt.plot(record_cn[:,0],record_cn[:,9],label='CN_8')
            #plt.plot(record_cn[:,0],record_cn[:,-1],label='CN_9')
            png_filename = prefix +'T_VS_CN_k'+'index'+str_index+'egcut'+'.png'
        else:
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,1],label='CN_0')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,2],label='CN_1')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,3],label='CN_2')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,4],label='CN_3')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,5],label='CN_4')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,6],label='CN_5')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,7],label='CN_6')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,8],label='CN_7')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,9],label='CN_8')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,-1],label='CN_9')
            png_filename = prefix +'T_VS_CN_k_tcut'+'index'+str_index+'egcut'+'.png'
        plt.legend()
        plt.title('CN_k '+'index:'+str_index)
        plt.xlabel('time(steps)')
        plt.ylabel('CN_k(1)')
        #plt.show()
        plt.savefig(png_filename)
        record_filename = prefix +'T_VS_CN_k_cut'+'index'+str_index+'.txt'
        numpy.savetxt(record_filename,record_cn)
        plt.close()

class compare_diagram:
    def __init__(self):
        pass

    def workflow_simu_to_bond_plot(self):
        table_return = self.select_a_list_of_simulations_from_mysql()
        n_simu = len(table_return)

        for i in range(n_simu):
            txyz,prefix_simu, str_simu_index = self.from_gsd_to_txyz(table_return[i][0],table_return[i][6])
            lcr=table_return[i][2]
            prefix_simu = '/home/remote/Downloads/'
            self.show_bond_plot(xy=txyz[-1,:,0:2],folder_name=prefix_simu,str_index=str_simu_index,lcr=lcr)


    def select_a_list_of_simulations_from_mysql(self):
        R"""
        introduction:
            scan a list of lcrs in dynamical diagram to see 
            what the system meet when lcr is not so suitable.
        example:
            select * from pin_hex_to_honeycomb_part_klt_2m 
            where HarmonicK=800 and RandomSeed=8;
            +-----------+-----------+------------------------+------+
            | SimuIndex | HarmonicK | LinearCompressionRatio | kT   
            | Psi3     | Psi6     | RandomSeed |


        Caution:
            simulations where kT=0.1 are included when RandomSeed=9. 
            RandomSeed=0-8 is ok.
        """
        import opertateOnMysql as osql
        table_name = 'pin_hex_to_honeycomb_part_klt_2m'
        cond = 'where HarmonicK = 800 and RandomSeed=8'
        table_return = osql.getDataFromMysql(table_name=table_name,search_condition=cond)
        data = np.array(table_return)
        list_simu_index = data[:,0].astype(int)
        list_seed = data[:,-1].astype(int)
        #print(list_simu_index,list_seed)
        return table_return

    def from_gsd_to_txyz(self,simu_index,seed,io_only=True):
        import proceed_file as pf
        pfile = pf.proceed_file_shell()
        prefix_simu, str_simu_index = pfile.create_prefix_results(account='remote',simu_index=simu_index,seed=seed)
        if not io_only:
            self.prefix=prefix_simu
            #save data
            gsd_data = pf.proceed_gsd_file(account='remote',simu_index=simu_index,seed=seed)
            gsd_data.get_trajectory_data(self.prefix)
            gsd_data.get_trajectory_stable_data(self.prefix)
            
        #load data
        filename_txyz = prefix_simu+'txyz.npy'
        txyz = np.load(filename_txyz)
        #nframes =2000
        #file_txyz_npy = self.prefix+'txyz_stable.npy'#_stable
        #txyz_stable = np.load(file_txyz_npy)
        return txyz,prefix_simu, str_simu_index
        
    def show_bond_plot(self,xy,folder_name,str_index,lcr):
        import points_analysis_2D as pa
        a_frame = pa.static_points_analysis_2d(points=xy)
        i=2000
        bond_cut_off=6
        trap_lcr=lcr
        trap_filename='/home/remote/hoomd-examples_0/testhoneycomb3-8-12-part1'
        self.x_unit = '($\sigma$)'
        png_filename1 = folder_name +'bond_hist_index'+str_index+'_'+str(int(i))+'.png'#+"/"
        png_filename2 = folder_name +'bond_plot_1st_minima_index'+str_index+'_'+str(int(i))+'.png'#+"/"
        
        a_frame.get_first_minima_bond_length_distribution(lattice_constant=1,hist_cutoff=bond_cut_off,png_filename=png_filename1)#,png_filename=png_filename1
        #a_frame.draw_bonds_conditional_bond(check=[0.4, a_frame.bond_first_minima_left], png_filename=png_filename2,
        #                               show_traps=show_traps,LinearCompressionRatio=trap_lcr,trap_filename=trap_filename,
        #                                nb_change=ids,x_unit=self.x_unit)
        a_frame.draw_bonds_conditional_bond_oop(check=[0.4, a_frame.bond_first_minima_left], png_filename=png_filename2,
                                                xy_stable=xy,x_unit=self.x_unit,
                                        LinearCompressionRatio=trap_lcr, trap_filename=trap_filename)
    
    def select_a_list_of_simulations_from_mysql_line(self,table_name = 'pin_hex_to_honeycomb_part_klt_2m'):
        R"""
        introduction:
            scan a list of lcrs in dynamical diagram to see 
            what the system meet when lcr is not so suitable.
        example:
            select * from pin_hex_to_honeycomb_part_klt_2m 
            where HarmonicK=800 and kT=1;
            +-----------+-----------+------------------------+------+
            | SimuIndex | HarmonicK | LinearCompressionRatio | kT   
            | Psi3     | Psi6     | RandomSeed |
        """
        import opertateOnMysql as osql
        #table_name = 'pin_hex_to_honeycomb_part_klt_2m'
        cond = 'where HarmonicK=800 and kT=1'
        table_return = osql.getDataFromMysql(table_name=table_name,search_condition=cond)
        return table_return
        
    def workflow_simu_to_line_honey_part(self):
        import pandas as pd
        table_return = self.select_a_list_of_simulations_from_mysql_line()
        n_simu = len(table_return)
        #SimuIndex | HarmonicK | LinearCompressionRatio | kT   | Psi3     | Psi6     | RandomSeed
        columns_name = ["SimuIndex","HarmonicK","LinearCompressionRatio","kT","Psi3", "Psi6","RandomSeed"]
        df_return = pd.DataFrame(table_return)
        df_return.columns = columns_name
        simu_index = df_return["SimuIndex"].unique()
        record_lcr_psi3_std = np.zeros((9,3))
        for i in range(9):
            df_simu1 = df_return[df_return["SimuIndex"]==simu_index[i]]
            list_psi3 = df_simu1["Psi3"]
            av_psi3 = np.average(list_psi3)
            std_psi3 = np.std(list_psi3)
            list_lcr = df_simu1["LinearCompressionRatio"].values
            record_lcr_psi3_std[i]=[list_lcr[0],av_psi3,std_psi3]
        
        
        columns_name = ["LinearCompressionRatio","Psi3", "std"]
        df_lps = pd.DataFrame(record_lcr_psi3_std)
        df_lps.columns = columns_name
        #https://blog.csdn.net/MsSpark/article/details/83154128
        df_lps.sort_values(by='LinearCompressionRatio',ascending=True,inplace=True)
        record_lcr_psi3_std = df_lps.values

        fig,ax = plt.subplots()
        #ax.plot(record_lcr_psi3_std[:,0],record_lcr_psi3_std[:,1])
        ax.errorbar(record_lcr_psi3_std[:,0],record_lcr_psi3_std[:,1],record_lcr_psi3_std[:,2])
        ax.set_xlabel("lcr(1)")
        ax.set_ylabel("$\psi_3(1)$")
        ax.set_title("scan $U_{trap}=400 k_BT_m$ in diagram under honey_part")
        prefix = '/home/remote/Downloads/'
        png_filename = prefix+'honey_part_scan_k800.png'
        fig.savefig(png_filename)
        plt.close()
    
    def workflow_simu_to_line_honey(self):
        import pandas as pd
        table_return = self.select_a_list_of_simulations_from_mysql_line('pin_hex_to_honeycomb_klt_2m')
        n_simu = len(table_return)
        #SimuIndex | HarmonicK | LinearCompressionRatio | kT   | Psi3     | Psi6     | RandomSeed
        columns_name = ["SimuIndex","HarmonicK","LinearCompressionRatio","kT","Psi3", "Psi6","RandomSeed"]
        df_return = pd.DataFrame(table_return)
        df_return.columns = columns_name
        simu_index = df_return["SimuIndex"].unique()
        nr = len(simu_index)
        record_lcr_psi3_std = np.zeros((nr,3))
        for i in range(nr):
            df_simu1 = df_return[df_return["SimuIndex"]==simu_index[i]]
            list_psi3 = df_simu1["Psi3"]
            av_psi3 = np.average(list_psi3)
            std_psi3 = np.std(list_psi3)
            list_lcr = df_simu1["LinearCompressionRatio"].values
            record_lcr_psi3_std[i]=[list_lcr[0],av_psi3,std_psi3]
        
        
        columns_name = ["LinearCompressionRatio","Psi3", "std"]
        df_lps = pd.DataFrame(record_lcr_psi3_std)
        df_lps.columns = columns_name
        #https://blog.csdn.net/MsSpark/article/details/83154128
        df_lps.sort_values(by='LinearCompressionRatio',ascending=True,inplace=True)
        record_lcr_psi3_std = df_lps.values

        fig,ax = plt.subplots()
        #ax.plot(record_lcr_psi3_std[:,0],record_lcr_psi3_std[:,1])
        ax.errorbar(record_lcr_psi3_std[:,0],record_lcr_psi3_std[:,1],record_lcr_psi3_std[:,2])
        ax.set_xlabel("lcr(1)")
        ax.set_ylabel("$\psi_3(1)$")
        ax.set_title("scan $U_{trap}=400 k_BT_m$ in diagram under honey")
        prefix = '/home/remote/Downloads/'
        png_filename = prefix+'honey_scan_k800.png'
        fig.savefig(png_filename)
        plt.close()
    
    def simu_to_line(self,table_name):
        import pandas as pd
        table_return = self.select_a_list_of_simulations_from_mysql_line(table_name)
        n_simu = len(table_return)
        #SimuIndex | HarmonicK | LinearCompressionRatio | kT   | Psi3     | Psi6     | RandomSeed
        columns_name = ["SimuIndex","HarmonicK","LinearCompressionRatio","kT","Psi3", "Psi6","RandomSeed"]
        df_return = pd.DataFrame(table_return)
        df_return.columns = columns_name
        simu_index = df_return["SimuIndex"].unique()
        nr = len(simu_index)
        record_lcr_psi3_std = np.zeros((nr,3))
        for i in range(nr):
            df_simu1 = df_return[df_return["SimuIndex"]==simu_index[i]]
            list_psi3 = df_simu1["Psi3"]
            av_psi3 = np.average(list_psi3)
            std_psi3 = np.std(list_psi3)
            list_lcr = df_simu1["LinearCompressionRatio"].values
            record_lcr_psi3_std[i]=[list_lcr[0],av_psi3,std_psi3]
        
        columns_name = ["LinearCompressionRatio","Psi3", "std"]
        df_lps = pd.DataFrame(record_lcr_psi3_std)
        df_lps.columns = columns_name
        #https://blog.csdn.net/MsSpark/article/details/83154128
        df_lps.sort_values(by='LinearCompressionRatio',ascending=True,inplace=True)
        record_lcr_psi3_std = df_lps.values
        return record_lcr_psi3_std
        
    def compare_honey_and_part(self):
        R"""
        introduction:
            scan a list of lcrs in dynamical diagram to see 
            what the system meet when lcr is not so suitable.
        example:
            select * from pin_hex_to_honeycomb_part_klt_2m 
            where HarmonicK=800 and kT=1;
            +-----------+-----------+------------------------+------+
            | SimuIndex | HarmonicK | LinearCompressionRatio | kT   
            | Psi3     | Psi6     | RandomSeed |
        """
        record_lcr_psi3_std_h = self.simu_to_line('pin_hex_to_honeycomb_klt_2m')
        record_lcr_psi3_std_hp = self.simu_to_line('pin_hex_to_honeycomb_part_klt_2m')

        fig,ax = plt.subplots()
        #ax.plot(record_lcr_psi3_std[:,0],record_lcr_psi3_std[:,1])
        err = True#False #True
        start = 0#3
        if err:
            ax.errorbar(record_lcr_psi3_std_hp[:,0],record_lcr_psi3_std_hp[:,1],record_lcr_psi3_std_hp[:,2],label='honey_part')
            ax.errorbar(record_lcr_psi3_std_h[start:,0],record_lcr_psi3_std_h[start:,1],record_lcr_psi3_std_h[start:,2],label='honey')

        else:
            ax.plot(record_lcr_psi3_std_hp[:,0],record_lcr_psi3_std_hp[:,1],label='honey_part')
            ax.plot(record_lcr_psi3_std_h[start:,0],record_lcr_psi3_std_h[start:,1],label='honey')
        ax.legend()
        ax.set_xlabel("lcr(1)")
        ax.set_ylabel("$\psi_3(1)$")
        ax.set_title("scan $U_{trap}=400 k_BT_m$ in diagram under honey&part")
        prefix = '/home/remote/Downloads/'
        png_filename = prefix+'honey+part_scan_k800_err.png'
        fig.savefig(png_filename)
        plt.close()

    def plot_hp_diagram_cn3newmethod(self):
        #trial 3
        #get simu index from mysql
        #lcr0.8 seed9 scan k
        import opertateOnMysql as osql
        table_name = 'pin_hex_to_honeycomb_part_klt_2m'
        cond = 'where HarmonicK=800 and kT=1'#and RandomSeed=9
        table_return = osql.getDataFromMysql(table_name=table_name,search_condition=cond)
        
        #sort index
        import pandas as pd
        columns_name = ["SimuIndex","HarmonicK","LinearCompressionRatio","kT","Psi3", "Psi6","RandomSeed"]
        df_return = pd.DataFrame(table_return)
        df_return.columns = columns_name
        df_return.sort_values(by='LinearCompressionRatio',ascending=True,inplace=True)

        #get (steps,cn3)
        list_simu_index = df_return["SimuIndex"].values
        list_seed = df_return["RandomSeed"].values
        account = 'remote'
        import points_analysis_2D as pa
        prefix='/home/'+account+'/Downloads/'
        n_simu = len(list_simu_index)
        record_cn3bond = []
        record_cn3ridge = []
        for i in range(n_simu):
            data_filename=prefix+'index'+str(list_simu_index[i])+"_"+str(list_seed[i])
            obj_of_simu_index = pa.static_points_analysis_2d(filename=data_filename)
            obj_of_simu_index.get_coordination_number_conditional()
            ccn=obj_of_simu_index.count_coordination_ratio
            record_cn3bond.append(ccn[3])
            obj_of_simu_index.get_coordination_number_conditional(method="ridge_length_method")
            ccn=obj_of_simu_index.count_coordination_ratio
            record_cn3ridge.append(ccn[3])
        
        #plot
        plot_cn3=False
        if plot_cn3:
            list_lcr = df_return["LinearCompressionRatio"].values
            fig,ax = plt.subplots()
            ax.plot(list_lcr,record_cn3bond)
        
            ax.set_xlabel("LCR(1)")
            ax.set_ylabel("$CN_3(1)$")
            ax.set_title("scan LCR in diagram under honey part")
            png_filename = prefix+'cn3ri_honey_part_k800seed9_scan_lcr.png'
            fig.savefig(png_filename)
            plt.close()
            
        #save data
        df_fn = pd.DataFrame(record_cn3bond)
        df_return['cn3bond'] = df_fn
        df_fn2 = pd.DataFrame(record_cn3ridge)
        df_return['cn3ridge'] = df_fn2
        csv_filename = prefix + 'cn3bond_ridge_k800_scan_lcr.csv'
        pd.DataFrame.to_csv(df_return,csv_filename) 
    
    def compare_hp_diagram_cn3newmethod_bond_ridge(self):
        #load data
        import pandas as pd
        account = 'remote'
        prefix='/home/'+account+'/Downloads/'
        csv_filename = prefix + 'cn3bond_ridge_k800_scan_lcr.csv'
        df_return = pd.read_csv(csv_filename)
        #columns_name = ["SimuIndex","HarmonicK","LinearCompressionRatio","kT","Psi3", "Psi6",
        #                "RandomSeed","cn3bond","cn3ridge"]
        #average of 10 seeds
        df_return_single_seed = df_return[df_return["RandomSeed"]==9]
        
        #list_simu_index = df_return["SimuIndex"].values
        #list_simu_index = np.unique(list_simu_index) 
        #columns_name = ["SimuIndex","LinearCompressionRatio","cn3bond","cn3ridge"]
        df_total = df_return_single_seed[["SimuIndex","LinearCompressionRatio"]]#table_total = pd.DataFrame(list_simu_index) 
        list_simu_index = df_total["SimuIndex"].values
        #columns_name = ["SimuIndex"]
        #table_total.columns = columns_name

        cn3bond_avg = []
        cn3bond_std = []
        cn3ridge_avg = []
        cn3ridge_std = []
        for index1 in list_simu_index:
            group = df_return[df_return['SimuIndex']==index1] 
            list_cn3bond = group["cn3bond"].values
            list_cn3ridge = group["cn3ridge"].values
            cn3bond_avg.append( np.average(list_cn3bond))
            cn3bond_std.append( np.std(list_cn3bond))
            cn3ridge_avg.append( np.average(list_cn3ridge))
            cn3ridge_std.append( np.std(list_cn3ridge))
            
        df_total["cn3bond_avg"] = cn3bond_avg
        df_total["cn3bond_std"] = cn3bond_std
        df_total["cn3ridge_avg"] = cn3ridge_avg
        df_total["cn3ridge_std"] = cn3ridge_std
        #save data
        csv_filename = prefix + 'cn3bond_ridge_avg_k800_scan_lcr.csv'
        pd.DataFrame.to_csv(df_total,csv_filename) 

        #plot
        plot_cn3=True
        if plot_cn3:
            list_lcr = df_total["LinearCompressionRatio"].values
            list_cn3bond_avg = df_total["cn3bond_avg"].values
            list_cn3bond_std = df_total["cn3bond_std"].values
            list_cn3ridge_avg = df_total["cn3ridge_avg"].values
            list_cn3ridge_std = df_total["cn3ridge_std"].values
            fig,ax = plt.subplots()
            ax.errorbar(list_lcr,list_cn3bond_avg,list_cn3bond_std,label='cn3bond')
            ax.errorbar(list_lcr,list_cn3ridge_avg,list_cn3ridge_std,label='cn3ridge')
            ax.legend()
            ax.set_xlabel("LCR(1)")
            ax.set_ylabel("$CN_3(1)$")
            ax.set_title("scan LCR in diagram under honey part")
            png_filename = prefix+'cn3compare_honey_part_k800seed9_scan_lcr.png'
            fig.savefig(png_filename)
            plt.close()

class compare_diagram_dynamics:
    def __init__(self):
        pass

    def plot_cn3(self,simu_index,seed,account='remote'):
        import data_analysis_cycle as dac
        record_filename = dac.save_from_gsd_to_cn3(simu_index=simu_index, seed=seed,
                        coordination_number=True,account=account)

    def plot_hp_cn3_scan_k(self):
        #trial 3
        #get simu index from mysql
        #lcr0.8 seed9 scan k
        import opertateOnMysql as osql
        table_name = 'pin_hex_to_honeycomb_part_klt_2m'
        cond_lcr = 'LinearCompressionRatio>0.795 and LinearCompressionRatio<0.805'
        cond = 'where kT=1 and RandomSeed=9 and SimuIndex <4300 and '+cond_lcr
        table_return = osql.getDataFromMysql(table_name=table_name,search_condition=cond)
        
        #sort index
        import pandas as pd
        columns_name = ["SimuIndex","HarmonicK","LinearCompressionRatio","kT","Psi3", "Psi6","RandomSeed"]
        df_return = pd.DataFrame(table_return)
        df_return.columns = columns_name
        df_return.sort_values(by='HarmonicK',ascending=True,inplace=True)

        #get (steps,cn3)
        list_simu_index = df_return["SimuIndex"].values
        list_seed = df_return["RandomSeed"].values
        account = 'remote'
        import data_analysis_cycle as dac
        n_simu = len(list_simu_index)
        list_record_filename = []
        for i in range(n_simu):
            record_filename = dac.save_from_gsd_to_cn3(simu_index=list_simu_index[i], seed=list_seed[i],
                        coordination_number=True,account=account)
            list_record_filename.append(record_filename)
        
        #plot
        n_fname = len(list_record_filename)
        list_k = df_return["HarmonicK"].values
        frame_cut = 2000
        fig,ax = plt.subplots()
        for i in range(n_fname):
            file_name_npy = list_record_filename[i]
            record_cn = np.load(file_name_npy)
            ax.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,4],label=str(list_k[i]/2))
    
        ax.legend()
        ax.set_xlabel("time(steps)")
        ax.set_ylabel("$CN_3(1)$")
        ax.set_title("scan lcr=0.8 seed=9 in diagram under honey part")
        prefix='/home/'+account+'/Downloads/'
        png_filename = prefix+'cn3_honey_part_lcr080seed9_scan_k.png'
        fig.savefig(png_filename)
        plt.close()
        
        #save data
        df_fn = pd.DataFrame(list_record_filename)
        df_return['record_filename'] = df_fn
        csv_filename = prefix + 'lcr080seed9_scan_k.csv'
        pd.DataFrame.to_csv(df_return,csv_filename) 

    def plot_hp_cn3_scan_k_fast(self,account = 'remote'):
        #trial 3
        #get simu index from mysql
        #lcr0.8 seed9 scan k
        
        #read file
        columns_name = ["SimuIndex","HarmonicK","LinearCompressionRatio",
        "kT","Psi3", "Psi6","RandomSeed","record_filename"]
        #account = 'remote'
        prefix='/home/'+account+'/Downloads/'
        csv_filename = prefix + 'lcr080seed9_scan_k.csv'
        #get data
        import pandas as pd
        df_return = pd.read_csv(csv_filename)
        list_record_filename = df_return["record_filename"].values        
        #plot
        n_fname = len(list_record_filename)
        list_k = df_return["HarmonicK"].values
        frame_cut = 2000
        fig,ax = plt.subplots()
        for i in range(n_fname):
            file_name_npy = list_record_filename[i]
            record_cn = np.load(file_name_npy)
            ax.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,4],label=str(list_k[i]/2))#semilogx
    
        ax.legend()
        ax.set_xlabel("time(step)")
        ax.set_ylabel("$CN_3(1)$")
        ax.set_title("scan lcr=0.8 seed=9 in diagram under honey part")
        prefix='/home/'+account+'/Downloads/'
        png_filename = prefix+'cn3_honey_part_lcr080seed9_scan_k.png'#_logx
        plt.show()
        fig.savefig(png_filename)
        plt.close()

    def plot_honey_cn3_scan_k(self):
        #trial 3
        #get simu index from mysql
        #lcr0.8 seed9 scan k
        import opertateOnMysql as osql
        table_name = 'pin_hex_to_honeycomb_klt_2m'
        cond_lcr = 'LinearCompressionRatio>0.795 and LinearCompressionRatio<0.805'
        cond = 'where kT=1 and RandomSeed=9 and SimuIndex <4646 and '+cond_lcr
        table_return = osql.getDataFromMysql(table_name=table_name,search_condition=cond)
        
        #sort index
        import pandas as pd
        columns_name = ["SimuIndex","HarmonicK","LinearCompressionRatio","kT","Psi3", "Psi6","RandomSeed"]
        df_return = pd.DataFrame(table_return)
        df_return.columns = columns_name
        df_return.sort_values(by='HarmonicK',ascending=True,inplace=True)

        #get (steps,cn3)
        list_simu_index = df_return["SimuIndex"].values
        list_seed = df_return["RandomSeed"].values
        account = 'remote'
        import data_analysis_cycle as dac
        n_simu = len(list_simu_index)
        list_record_filename = []
        for i in range(n_simu):
            record_filename = dac.save_from_gsd_to_cn3(simu_index=list_simu_index[i], seed=list_seed[i],
                        coordination_number=True,account=account)
            list_record_filename.append(record_filename)
        
        #plot
        n_fname = len(list_record_filename)
        list_k = df_return["HarmonicK"].values
        frame_cut = 2000
        fig,ax = plt.subplots()
        for i in range(n_fname):
            file_name_npy = list_record_filename[i]
            record_cn = np.load(file_name_npy)
            ax.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,4],label=str(list_k[i]/2))
    
        ax.legend()
        ax.set_xlabel("time(steps)")
        ax.set_ylabel("$CN_3(1)$")
        ax.set_title("scan lcr=0.8 seed=9 in diagram under honey")
        prefix='/home/'+account+'/Downloads/'
        png_filename = prefix+'cn3_honey_lcr080seed9_scan_k.png'
        fig.savefig(png_filename)
        plt.close()
        
        #save data
        df_fn = pd.DataFrame(list_record_filename)
        df_return['record_filename'] = df_fn
        csv_filename = prefix + 'h_lcr080seed9_scan_k.csv'
        pd.DataFrame.to_csv(df_return,csv_filename) 

    def plot_h_hp_cn3_scan_k_fast(self,account = 'tplab'):
        """
        import workflow_analysis as wa
        spd = wa.compare_diagram_dynamics()#show_polygon_dye() 
        spd.plot_h_hp_cn3_scan_k_fast()#get_data_cn3_scan_k_fast()#plot_polygon_bond_xylim()
        """
        #trial 3
        #get simu index from mysql
        #lcr0.8 seed9 scan k
        fig,ax = plt.subplots()

        n_coarse=10
        self.get_data_cn3_scan_k_fast(ax,'lcr080seed9_scan_k.csv','-<',n_coarse)
        self.get_data_cn3_scan_k_fast(ax,'h_lcr080seed9_scan_k.csv','-H',n_coarse)
        
        ax.legend(loc='lower right',ncol=2)#
        ax.set_xlabel("time(steps)")
        ax.set_ylabel("$CN_3(1)$")
        ax.set_title("scan lcr=0.8 seed=9 in diagram under honey & part")
        #account = 'tplab'#'remote'
        prefix='/home/'+account+'/Downloads/'
        png_filename = prefix+'cn3_honey_and_part_lcr080seed9_scan_k_plot_coarse.pdf'#_logx _select3
        plt.show()
        fig.savefig(png_filename)
        plt.close()
    
    def get_data_cn3_scan_k_fast(self,ax,csv_filename,marker,n_coarse=None):
        R"""
        ax:
            the axe to plot lines.
        csv_filename:
            the .csv formatted filename containing data. 
        marker: 
            '<': 'triangle_left';'H': 'hexagon2'
            see matplotlib.markers
        """
        #read file
        """
        columns_name = ["SimuIndex","HarmonicK","LinearCompressionRatio",
        "kT","Psi3", "Psi6","RandomSeed","record_filename"]
        """
        account = 'remote'
        prefix='/home/'+account+'/Downloads/'#+'cn3t/'
        csv_fname = prefix + csv_filename#'lcr080seed9_scan_k.csv'
        #get data
        import pandas as pd
        df_return = pd.read_csv(csv_fname)
        list_record_filename = df_return["record_filename"].values        
        #plot
        n_fname = len(list_record_filename)
        list_k = df_return["HarmonicK"].values
        frame_cut = 2000
        #fig,ax = plt.subplots()
        for i in [0,2,5]:#[0,2,5]:#range(n_fname)：
            file_name_npy = list_record_filename[i]
            record_cn = np.load(file_name_npy)
            
            if n_coarse is None:
                ax.semilogx(record_cn[0:frame_cut,0],record_cn[0:frame_cut,4],marker,label=str(list_k[i]/2))#semilogx
            else:
                list_index,record_cn_coarse = self.data_decorating(record_cn[:frame_cut,4],n_coarse)
                ax.plot(record_cn[list_index,0],record_cn_coarse,marker,label=str(list_k[i]/2))#semilogx

    def plot_hp_cn3_scan_seed(self):
        #trial 2
        #get simu index from mysql
        ##lcr0.81 k700 scan seed
        import opertateOnMysql as osql
        table_name = 'pin_hex_to_honeycomb_part_klt_2m'
        cond_lcr = 'LinearCompressionRatio>0.805 and LinearCompressionRatio<0.815'
        cond = 'where HarmonicK=700 and kT=1 and '+cond_lcr
        table_return = osql.getDataFromMysql(table_name=table_name,search_condition=cond)
        
        #sort index
        import pandas as pd
        columns_name = ["SimuIndex","HarmonicK","LinearCompressionRatio","kT","Psi3", "Psi6","RandomSeed"]
        df_return = pd.DataFrame(table_return)
        df_return.columns = columns_name
        df_return.sort_values(by='RandomSeed',ascending=True,inplace=True)#edit

        #get (steps,cn3)
        list_simu_index = df_return["SimuIndex"].values
        list_seed = df_return["RandomSeed"].values
        account = 'remote'
        import data_analysis_cycle as dac
        n_simu = len(list_simu_index)
        list_record_filename = []
        for i in range(n_simu):
            record_filename = dac.save_from_gsd_to_cn3(simu_index=list_simu_index[i], seed=list_seed[i],
                        coordination_number=True,account=account)
            list_record_filename.append(record_filename)
        
        #plot
        n_fname = len(list_record_filename)
        frame_cut = 2000
        fig,ax = plt.subplots()
        for i in range(n_fname):
            file_name_npy = list_record_filename[i]
            record_cn = np.load(file_name_npy)
            ax.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,4],label=str(list_seed[i]))
        
        ax.legend()
        ax.set_xlabel("time(steps)")
        ax.set_ylabel("$CN_3(1)$")
        ax.set_title("scan $U_{trap}=350 k_BT_m$ lcr=0.81 in diagram under honey part")
        prefix='/home/'+account+'/Downloads/'
        png_filename = prefix+'cn3_honey_part_lcr081k700_scan_seed.png'
        fig.savefig(png_filename)
        plt.close()
        
        #save data
        df_fn = pd.DataFrame(list_record_filename)
        df_return['record_filename'] = df_fn
        csv_filename = prefix + 'lcr081k700_scan_seed.csv'
        pd.DataFrame.to_csv(df_return,csv_filename) 

    def plot_hp_cn3_scan_lcr(self):
        #trial 1
        #get simu index from mysql
        import opertateOnMysql as osql
        table_name = 'pin_hex_to_honeycomb_part_klt_2m'
        cond = 'where HarmonicK=800 and kT=1 and RandomSeed=9'
        table_return = osql.getDataFromMysql(table_name=table_name,search_condition=cond)
        
        #sort index
        import pandas as pd
        columns_name = ["SimuIndex","HarmonicK","LinearCompressionRatio","kT","Psi3", "Psi6","RandomSeed"]
        df_return = pd.DataFrame(table_return)
        df_return.columns = columns_name
        df_return.sort_values(by='LinearCompressionRatio',ascending=True,inplace=True)

        #get (steps,cn3)
        list_simu_index = df_return["SimuIndex"].values
        seed = 9
        account = 'remote'
        import data_analysis_cycle as dac
        n_simu = len(list_simu_index)
        list_record_filename = []
        for simu_index in list_simu_index:
            record_filename = dac.save_from_gsd_to_cn3(simu_index=simu_index, seed=seed,
                        coordination_number=True,account=account)
            list_record_filename.append(record_filename)
        
        n_fname = len(list_record_filename)
        list_lcr = df_return["LinearCompressionRatio"].values
        frame_cut = 2000
        fig,ax = plt.subplots()
        for i in range(n_fname):
            file_name_npy = list_record_filename[i]
            record_cn = np.load(file_name_npy)
            ax.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,4],label=str(list_lcr[i]))
        
        ax.legend()
        ax.set_xlabel("time(steps)")
        ax.set_ylabel("$CN_3(1)$")
        ax.set_title("scan $U_{trap}=400 k_BT_m$ in diagram under honey part")
        prefix = '/home/remote/Downloads/'
        png_filename = prefix+'cn3_honey_part_k800seed9_scan_lcr.png'
        fig.savefig(png_filename)
        plt.close()

        #save data
        df_fn = pd.DataFrame(list_record_filename)
        df_return['record_filename'] = df_fn
        csv_filename = prefix + 'k800seed9_scan_lcr.csv'
        pd.DataFrame.to_csv(df_return,csv_filename) 

    def data_decorating(self,data,coarse_grain_to_n_points=10,navg=5):
        R"""
        parameters:
        data: 1d array with Ndata_points (float)
        Ndata_points: num of data points
        coarse_grain_to_n_points(n2c): (int)coarse_grain_to_n_points
        Navg(navg): num of data points to average, 
            positive odd integral only(1,3,5...)
        list_index:
        data_decorated: n2c rows of data points averaged with (Navg-1) neighbors. 
        
        introduction:
        the 1st data point is not in need of average, but others are. 
        hence the index of last data point must be set smaller than 
        Ndata_points - (Navg-1)/2

        caution:
        the shape of data points is not log function, 
        so fitting is of nonsense.
        """
        n2c = coarse_grain_to_n_points
        sz = np.shape(data)#rows for Ndata
        log_max = np.log(sz[0]) 
        list_log_index = np.linspace(0,log_max,n2c)
        list_index = np.zeros((n2c,),dtype=int)
        list_index[1:] = np.exp(list_log_index[1:]).astype(int)
        list_index[-1] = list_index[-1]-(navg-1)/2
        data_decorated = np.zeros((n2c,))
        print(data[sz[0]-2:])
        for i in range(n2c):
            if i==0:
                data_decorated[i] = data[0] 
                """
                elif i==n2c-1:
                in_st = list_index[i]-(navg-1)
                #in_ed = list_index[i]
                print(i,data[in_st:])
                data_decorated[i] =np.average(data[in_st:])
                """
            else:
                in_st = list_index[i]-(navg-1)/2
                in_ed = list_index[i]+(navg-1)/2
                in_st = in_st.astype(int)
                in_ed = in_ed.astype(int)+1#for numpy +1 is of necessity
                print(i,data[in_st:in_ed])
                data_decorated[i] =np.average(data[in_st:in_ed]) 
        return list_index,data_decorated

    
