import symmetry_transformation_v4_3.analysis_controller as ac
import numpy as np
class analyze_a_gsd_file:
    def __init__(self):
        pass

    def analyze_gsd_files_and_print_6013(self):
        n_file = 30
        list_index = np.linspace(6013,6042,n_file,dtype=int)
        list_gauss_epsilon = np.linspace(-3,-90,n_file,dtype=int)
        list_cn3 = np.zeros((n_file,))
        seed = 19
        print(seed)
        for i in range(n_file):
            simu_index = list_index[i]
            s2g = ac.get_a_gsd_from_setup()
            s2g.set_file_parameters(simu_index,seed)
            
            g2d = ac.get_data_from_a_gsd_frame(s2g.gsd_data[-1])
            list_cn3[i] = g2d.get_cn3_from_a_gsd_frame()
            print(list_gauss_epsilon[i],list_cn3[i])
    
    def analyze_gsd_files_and_print_6043(self):
        n_file = 10
        list_index = np.linspace(6043,6052,n_file,dtype=int)
        list_gauss_epsilon = np.linspace(-100,-1000,n_file,dtype=int)
        list_cn3 = np.zeros((n_file,))
        seed = 19
        print(seed)
        for i in range(n_file):
            simu_index = list_index[i]
            s2g = ac.get_a_gsd_from_setup()
            s2g.set_file_parameters(simu_index,seed)
            
            g2d = ac.get_data_from_a_gsd_frame(s2g.gsd_data[-1])
            list_cn3[i] = g2d.get_cn3_from_a_gsd_frame()
            print(list_gauss_epsilon[i],list_cn3[i])

    def analyze_gsd_files_and_record_as_csv(self):
        import pandas as pd
        prefix_write = "/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/"
        
        output_file_csv = prefix_write+'honeycomb_pin_scan_k100_1000'+'.csv'
        n_file = 10
        list_index = np.linspace(5793,5802,n_file,dtype=int)
        list_lcr = [0.81]*n_file
        list_trap_gauss_epsilon = np.linspace(-100,-1000,n_file,dtype=int)
        list_cn3 = np.zeros((n_file,))
        #df_all = pd.DataFrame()
        for seed in [0,1,2]:#,3,4,7,8,9
            #seed = 7
            list_seed_record = [seed]*n_file
            for i in range(n_file):
                simu_index = list_index[i]
                s2g = ac.get_a_gsd_from_setup()
                s2g.set_file_parameters(simu_index,seed)
                
                g2d = ac.get_data_from_a_gsd_frame()
                list_cn3[i] = g2d.get_cn3_from_a_gsd_frame(s2g.gsd_data[-1])
            
            
            df = pd.DataFrame(list_index,columns=['simu_index'])
            df['seed'] = list_seed_record
            df['lcr'] = list_lcr
            df['trap_gauss_epsilon'] = list_trap_gauss_epsilon
            df['cn3'] = list_cn3
            if seed == 0:
                df_all = df
            elif seed > 0:
                df_all = pd.concat([df_all,df])#.append(df)
        pd.DataFrame.to_csv(df_all,output_file_csv)

    def analyze_gsd_files_and_print(self):
        n_file = 20
        list_index = np.linspace(5773,5792,n_file,dtype=int)
        list_gauss_epsilon = np.linspace(-3,-60,n_file,dtype=int)
        list_cn3 = np.zeros((n_file,))
        seed = 7
        print(seed)
        for i in range(n_file):
            simu_index = list_index[i]
            s2g = ac.get_a_gsd_from_setup()
            s2g.set_file_parameters(simu_index,seed)
            
            g2d = ac.get_data_from_a_gsd_frame()
            list_cn3[i] = g2d.get_cn3_from_a_gsd_frame(s2g.gsd_data[-1])
            print(list_cn3[i])

    def get_bond_plot_from_a_gsd(self):
        ss = ac.get_a_gsd_from_setup()
        for index1 in np.linspace(5971,5974,4):#5950,5959
            ss.set_file_parameters(index1,9)
            ss.get_gsd_data_from_file()
            lf = ac.get_data_from_a_gsd_frame(ss.gsd_data[-1])
            lf.get_bonds_from_a_gsd_frame()
    
    def get_bond_plot_with_time_from_a_gsd(self,index1):
        ss = ac.get_a_gsd_from_setup()
        ss.set_file_parameters(index1,9)
        ss.get_gsd_data_from_file()
        num_list = [100,500,1000,1500]#np.linspace(10,1999,10,dtype=int)
        for i in num_list:#range(2000-1):
            print(ss.gsd_data[i].configuration.box)
            #print(ss.gsd_data[i].configuration.box.dimensions)
            lf = ac.get_data_from_a_gsd_frame(ss.gsd_data[i])
            png_filename = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'+str(ss.simu_index)+'_'+str(i)+'.png'
            lf.get_bonds_png_from_a_gsd_frame(png_filename)
        
class analyze_a_series_of_gsd_file:
    def __init__(self):
        pass

    def get_cn3s_from_mysql_honeycomb(self):
        R"""
            lcr0.81 honeycomb pin
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        ggfs = ac.get_gsds_from_mysql_or_csv()
        record = ggfs.get_record_from_sql_by_lcr()
        #record_np = np.array(record)
        gsds = ggfs.get_gsds_from_mysql_record(record)
        n_simu = len(gsds)
        list_cn3 = np.zeros((n_simu,))
        #list_cn3[:,:2] =  record_np[:,2:4]
        gf = ac.get_data_from_a_gsd_frame()
        for i in range(n_simu):
            frame = ggfs.get_last_frame_from_gsd(gsds[i])
            tune_dis = record[i][2]*3#list_cn3[i,0]*3
            cn3 = gf.get_cn3_from_a_gsd_frame(frame,tune_dis)
            list_cn3[i] = cn3
            #del gf
        
        """
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use(backend="QtAgg")#Backend agg is non-interactive backend. Turning interactive mode off. 'QtAgg' is interactive mode
        fig,ax = plt.subplots()
        ax.plot(list_lcr_k_cn3[:,1],list_lcr_k_cn3[:,2])
        #ax.plot()
        ax.set_aspect('equal','box')
        plt.show()
        plt.close()"""

        #save as csv
        import pandas as pd
        df = pd.DataFrame(record,columns=[ 'simu_index','seed',"lcr",'trap_gauss_epsilon','temperature'])
        df['cn3'] = list_cn3
        #df.sort_values(by=['seed','trap_gauss_epsilon'])
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/'
        csv_filename = prefix_write + 'honeycomb_pin_scan_k3_1000.csv'
        pd.DataFrame.to_csv(df,csv_filename)
        return list_cn3
    
    def get_cn3s_from_mysql_honeycomb_part(self):
        R"""
            lcr0.81 honeycomb pin
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        ggfs = ac.get_gsds_from_mysql_or_csv()
        record = ggfs.get_record_from_sql_by_lcr(lcr1=0.80,table_name='pin_hex_to_honeycomb_part_klt_2m_gauss')
        #record_np = np.array(record)
        gsds = ggfs.get_gsds_from_mysql_record(record)
        n_simu = len(gsds)
        list_cn3 = np.zeros((n_simu,))
        #list_cn3[:,:2] =  record_np[:,2:4]

        for i in range(n_simu):
            frame = ggfs.get_last_frame_from_gsd(gsds[i])
            tune_dis = record[i][2]*3#list_cn3[i,0]*3
            gf = ac.get_data_from_a_gsd_frame(frame)#error:missing last frame
            cn3 = gf.get_cn3_from_a_gsd_frame(frame,tune_dis)
            list_cn3[i] = cn3
            del gf
        
        #save as csv
        import pandas as pd
        df = pd.DataFrame(record,columns=[ 'simu_index','seed',"lcr",'trap_gauss_epsilon','temperature'])
        df['cn3'] = list_cn3


        """import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use(backend="QtAgg")#Backend agg is non-interactive backend. Turning interactive mode off. 'QtAgg' is interactive mode
        fig,ax = plt.subplots()
        ax.semilogx(-df['trap_gauss_epsilon'],list_cn3)
        #ax.plot()
        #ax.set_aspect('equal','box')
        plt.show()
        plt.close()"""

        #df.sort_values(by=['seed','trap_gauss_epsilon'])
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
        csv_filename = prefix_write + 'honeycomb_part_pin_lcr80_scan_k3_1000.csv'
        pd.DataFrame.to_csv(df,csv_filename)
        return list_cn3
    
    def get_cn3s_from_csv_honeycomb_part_gauss_eq(self,output_file_csv):
        R"""
            lcr0.81 honeycomb pin
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        """prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
        index1 = 6013
        output_file_csv = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_'+str(int(index1))+'.csv'"""
        ggfs = ac.get_gsds_from_mysql_or_csv()
        gsds = ggfs.get_record_from_csv(output_file_csv)
        n_simu = len(gsds)
        list_cn3 = np.zeros((n_simu,))
        #list_cn3[:,:2] =  record_np[:,2:4]

        for i in range(n_simu):
            frame = ggfs.get_last_frame_from_gsd(gsds[i])
            lcr1 = ggfs.record['lcr'].values[i]
            tune_dis = lcr1*3#list_cn3[i,0]*3
            gf = ac.get_data_from_a_gsd_frame(frame)#error:missing last frame
            cn3 = gf.get_cn3_from_a_gsd_frame(tune_dis)
            list_cn3[i] = cn3
            del gf
        
        #save as csv
        import pandas as pd
        df = pd.DataFrame(ggfs.record,columns=[ 'simu_index','seed',"lcr",'trap_gauss_epsilon','temperature'])
        df['cn3'] = list_cn3
        df['U_eq'] = df['trap_gauss_epsilon'].values*0.99613
        """import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use(backend="QtAgg")#Backend agg is non-interactive backend. Turning interactive mode off. 'QtAgg' is interactive mode
        fig,ax = plt.subplots()
        ax.semilogx(-df['trap_gauss_epsilon'],list_cn3)
        #ax.plot()
        #ax.set_aspect('equal','box')
        plt.show()
        plt.close()"""
        #df.sort_values(by=['seed','trap_gauss_epsilon'])
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
        csv_filename = prefix_write+'pin_hex_to_honeycomb_part_klt_2m_gauss_6013_19_res.csv'
        #csv_filename = prefix_write + 'honeycomb_part_pin_lcr81_scan_k3_1000_fill_brownian.csv'
        pd.DataFrame.to_csv(df,csv_filename)
        #print(df.head(35))
        return list_cn3

    def get_cn3s_from_csv_files(self):
        import symmetry_transformation_v4_3.analysis_controller as ac
        #import symmetry_transformation_v4_3.simulation_controller as sc
        #sct = sc.simulation_controller_honeycomb_part_traps()
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/'
        output_file_csv = prefix_write + 'pin_hex_to_honeycomb_klt_2m_gauss_5773_5812.csv'

        #print(cm.csv_merged.shape)
        #sct.generate_initial_state_hexagonal_particle_honeycomb_part_trap_scan_csv(output_file_csv)
        gg = ac.get_gsds_from_mysql_or_csv()
        gsds_filename = gg.get_record_from_csv(output_file_csv)
        """lcrs = gg.record['lcr'].values
        ks = gg.record['trap_gauss_epsilon'].values"""
        pdata = ac.get_a_gsd_from_setup()
        cn3s = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            pdata.get_gsd_data_from_filename(gsds_filename[i])
            gdata = ac.get_data_from_a_gsd_frame(pdata.gsd_data[-1])
            cn3 = gdata.get_cn3_from_a_gsd_frame()
            cn3s[i] = cn3
        gg.record['cn3'] = cn3s
        import pandas as pd
        pd.DataFrame.to_csv(gg.record,output_file_csv)

    def get_cn3s_from_two_csv_files(self):
        import symmetry_transformation_v4_3.analysis_controller as ac
        #import symmetry_transformation_v4_3.simulation_controller as sc
        #sct = sc.simulation_controller_honeycomb_part_traps()
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
        index1 = 6013
        output_file_csv = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_'+str(int(index1))+'_not8.csv'
        output_file_csv2 = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_b_'+str(int(index1))+'.csv'
        output_file_csv3 = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_'+str(int(index1))+'_09.csv'
        import proceed_file as pf
        cm = pf.merge_two_csvs(output_file_csv,output_file_csv2,output_file_csv3)
        #print(cm.csv_merged.shape)
        #sct.generate_initial_state_hexagonal_particle_honeycomb_part_trap_scan_csv(output_file_csv)
        gg = ac.get_gsds_from_mysql_or_csv()
        gsds_filename = gg.get_record_from_csv(output_file_csv3)
        """lcrs = gg.record['lcr'].values
        ks = gg.record['trap_gauss_epsilon'].values"""
        pdata = ac.get_a_gsd_from_setup()
        cn3s = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            pdata.get_gsd_data_from_file(gsds_filename[i])
            gdata = ac.get_data_from_a_gsd_frame(pdata.gsd_data[-1])
            cn3 = gdata.get_cn3_from_a_gsd_frame()
            cn3s[i] = cn3
        gg.record['cn3'] = cn3s
        import pandas as pd
        pd.DataFrame.to_csv(gg.record,output_file_csv3)

    def plot_from_csv(self):
        import pandas as pd
        """prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/'
        csv_filename = prefix_write + 'honeycomb_pin_scan_k3_1000.csv'"""
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
        csv_filename = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_6013_09.csv'
        df = pd.read_csv(csv_filename)        
        df.sort_values(by=['seed','trap_gauss_epsilon'],ascending=False,inplace=True)#new sorted df replace old unsorted df 
        #print(df.head(50))#df['seed']

        df_seed = df['seed'].values
        ar_seed = np.array(df_seed,dtype=int)
        list_seeds = np.unique(ar_seed) 
        
        
        n_seed = len(list_seeds)
        n_sim = int(len(df)/n_seed)
        record_cn3s = np.zeros((n_sim ,n_seed))
        record_k_cn3_std = np.zeros((n_sim ,3))
        for seed in list_seeds:#seed=0#
            new_df = df[df['seed']==seed]
            cn3 = np.array(new_df['cn3'].values)#cn3=[0.1,0.5,0.9]#
            if seed == 0:
                k = np.array(-new_df['trap_gauss_epsilon'].values)#k=[1,2,3]#
            record_cn3s[:,seed] = cn3
        record_k_cn3_std[:,0] = k
        record_k_cn3_std[:,1] = np.average(record_cn3s,axis=1)
        record_k_cn3_std[:,2] = np.std(record_cn3s,axis=1)
        
        #plot 10 seeds
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots()
        ax.semilogx(k,record_k_cn3_std[:,1])#,label=str(seed) x and y shoundn't be larger than 10 times, such as 10 vs 1, or plot will be invisible
        plt.errorbar(k,record_k_cn3_std[:,1],yerr=record_k_cn3_std[:,2])
        plt.show()
        print('ok')
        """
        #plot 10 seeds
        import matplotlib.pyplot as plt
        #import matplotlib
        #matplotlib.use(backend="QtAgg")#Backend agg is non-interactive backend. Turning interactive mode off. 'QtAgg' is interactive mode
        fig,ax = plt.subplots()
        for seed in list_seeds:#seed=0#
            new_df = df[df['seed']==seed]
            k = np.array(-new_df['trap_gauss_epsilon'].values)#k=[1,2,3]#
            cn3 = np.array(new_df['cn3'].values)#cn3=[0.1,0.5,0.9]#
            ax.semilogx(k,cn3)#,label=str(seed) x and y shoundn't be larger than 10 times, such as 10 vs 1, or plot will be invisible
        #ax.set_aspect('equal','box')
        plt.legend()
        plt.show()
        print('ok')"""
    
    def get_k_cn3_stds_from_csv(self,csv_filename):
        import pandas as pd
        df = pd.read_csv(csv_filename)        
        df.sort_values(by=['seed','trap_gauss_epsilon'],ascending=False,inplace=True)#new sorted df replace old unsorted df 
        #print(df.head(50))#df['seed']

        df_seed = df['seed'].values
        ar_seed = np.array(df_seed,dtype=int)
        list_seeds = np.unique(ar_seed) 
        
        
        n_seed = len(list_seeds)
        n_sim = int(len(df)/n_seed)
        record_cn3s = np.zeros((n_sim ,n_seed))
        record_k_cn3_std = np.zeros((n_sim ,3))
        for seed in list_seeds:#seed=0#
            new_df = df[df['seed']==seed]
            cn3 = np.array(new_df['cn3'].values)#cn3=[0.1,0.5,0.9]#
            if seed == 0:
                k = np.array(-new_df['trap_gauss_epsilon'].values)#k=[1,2,3]#
            record_cn3s[:,seed] = cn3
        record_k_cn3_std[:,0] = k
        record_k_cn3_std[:,1] = np.average(record_cn3s,axis=1)
        record_k_cn3_std[:,2] = np.std(record_cn3s,axis=1)
        #get averaged and std data
        #sta_df = df[df['seed']==0]
        
        return record_k_cn3_std
    
    def get_k_cn3_stds_csv_from_csv(self,csv_filename):
        import pandas as pd
        df = pd.read_csv(csv_filename)        
        df.sort_values(by=['seed','trap_gauss_epsilon'],ascending=False,inplace=True)#new sorted df replace old unsorted df 
        #print(df.head(50))#df['seed']

        df_seed = df['seed'].values
        ar_seed = np.array(df_seed,dtype=int)
        list_seeds = np.unique(ar_seed) 
        
        
        n_seed = len(list_seeds)
        n_sim = int(len(df)/n_seed)
        record_cn3s = np.zeros((n_sim ,n_seed))
        for seed in list_seeds:#seed=0#
            new_df = df[df['seed']==seed]
            cn3 = np.array(new_df['cn3'].values)#cn3=[0.1,0.5,0.9]#
            if seed == 0:
                k = np.array(-new_df['trap_gauss_epsilon'].values)#k=[1,2,3]#
            record_cn3s[:,seed] = cn3
        df_single_seed = df[df['seed']==0]
        df_single_seed['U_eq'] = k*0.39347#sigma=1,rcut=1, equalvalence U trap
        df_single_seed['cn3avg'] = np.average(record_cn3s,axis=1)
        df_single_seed['cn3std'] = np.std(record_cn3s,axis=1)
        
        
        pd.DataFrame.to_csv(df_single_seed,'pin_hex_to_honeycomb_klt_2m_gauss_5773_5812_res.csv')
        #'pin_hex_to_honeycomb_part_klt_2m_gauss_6013_09_res.csv' 
    
    def plot_k_cn3_stds(self,fig,ax,record_k_cn3_std,label=None):
        import matplotlib.pyplot as plt
        #fig,ax = plt.subplots()
        ax.semilogx(record_k_cn3_std[:,0],record_k_cn3_std[:,1])#,label=str(seed) x and y shoundn't be larger than 10 times, such as 10 vs 1, or plot will be invisible
        ax.errorbar(record_k_cn3_std[:,0],record_k_cn3_std[:,1],yerr=record_k_cn3_std[:,2],capsize=6,label=label)
        return fig,ax

    def get_2_k_cn3_stds_and_plot(self):
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/'
        csv_filename = prefix_write + 'honeycomb_pin_scan_k3_1000.csv'
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
        csv_filename2 = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_6013_09.csv'
        record_k_cn3_std1 = self.get_k_cn3_stds_from_csv(csv_filename)
        record_k_cn3_std2 = self.get_k_cn3_stds_from_csv(csv_filename2)
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots()
        fig,ax = self.plot_k_cn3_stds(fig,ax,record_k_cn3_std1,'hc')
        fig,ax = self.plot_k_cn3_stds(fig,ax,record_k_cn3_std2,'hp')
        plt.legend()
        plt.show()

    def get_bonds_from_simu_indices_type_n(self):
        R"""
            lcr0.81 honeycomb pin
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/'
        index0 = 6302
        for i in range(8):
            index1 = index0+10*i
            #gsd_file = prefix_write + 'trajectory_auto'+str(int(index1))+'_9.gsd'
            ggsd = ac.get_a_gsd_from_setup()
            ggsd.set_file_parameters(index1,9)
            ggsd.get_gsd_data_from_file()
            frame = ggsd.gsd_data[-1]
            
            gf = ac.get_data_from_a_gsd_frame(frame)#error:missing last frame
            gf.get_bonds_png_from_a_gsd_frame('depin_from_type_'+str(4+i)+'_part.png')
            del gf
    
    def get_bonds_from_simu_indices_type_n_from_csv(self):
        R"""
            column = ['simu_index','seed','lcr','trap_gauss_epsilon','temperature','type_n']
            manually set bond_length_max = 3*lcr0*1.2, set 1.2 is to avoid 1,414 and 1.73 bond.
        exp:
            import symmetry_transformation_v4_3.list_code_analysis as lca
            agf = lca.analyze_a_series_of_gsd_file()
            agp = agf.get_bonds_from_simu_indices_type_n_from_csv()
        """
        import pandas as pd
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_csv = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/'#type_n_depin/
        filename_csv = prefix_csv + 'pin_hex_to_type_n_part_klt_2m_gauss_6613.csv'
        #depin_type_n_from_type_n_part_klt_2m_gauss_6053,6293
        record_csv = pd.read_csv(filename_csv)
        list_type_n_to_watch = [7,10,11]
        ggsd = ac.get_a_gsd_from_setup()
        for type_n in list_type_n_to_watch:
            record_type_n = record_csv[record_csv['type_n'] == type_n]
            list_index = record_type_n['simu_index'].values
            lcr0 = record_type_n['lcr'].values[0]
            list_trap_gauss_epsilon = record_type_n['trap_gauss_epsilon'].values
            for i in range(len(list_index)):
                #gsd_file = prefix_write + 'trajectory_auto'+str(int(index1))+'_9.gsd'
                ggsd.set_file_parameters(list_index[i],9)
                #ggsd.get_gsd_data_from_file()
                frame = ggsd.gsd_data[-1]
                gf = ac.get_data_from_a_gsd_frame(frame)#error:missing last frame
                gf.get_given_bonds_png_from_a_gsd_frame('pin_hex_to_type_type_'+str(type_n)+'_part_'+str(-int(list_trap_gauss_epsilon[i]))+'_.png',3*lcr0*1.2)
                #gf.get_given_bonds_png_from_a_gsd_frame('depin_from_type_'+str(type_n)+'_part_'+str(-int(list_trap_gauss_epsilon[i]))+'_.png',3*lcr0*1.2)
                del gf

    def analyze_gsd_files_and_record_as_csv(self):
        import pandas as pd
        prefix_write = "/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/"
        
        output_file_csv = prefix_write+'honeycomb_pin_scan_k100_1000'+'.csv'
        n_file = 10
        list_index = np.linspace(5793,5802,n_file,dtype=int)
        list_lcr = [0.81]*n_file
        list_trap_gauss_epsilon = np.linspace(-100,-1000,n_file,dtype=int)
        list_cn3 = np.zeros((n_file,))
        #df_all = pd.DataFrame()
        for seed in [0,1,2]:#,3,4,7,8,9
            #seed = 7
            list_seed_record = [seed]*n_file
            for i in range(n_file):
                simu_index = list_index[i]
                s2g = ac.get_a_gsd_from_setup()
                s2g.set_file_parameters(simu_index,seed)
                
                g2d = ac.get_data_from_a_gsd_frame()
                list_cn3[i] = g2d.get_cn3_from_a_gsd_frame(s2g.gsd_data[-1])
            
            
            df = pd.DataFrame(list_index,columns=['simu_index'])
            df['seed'] = list_seed_record
            df['lcr'] = list_lcr
            df['trap_gauss_epsilon'] = list_trap_gauss_epsilon
            df['cn3'] = list_cn3
            if seed == 0:
                df_all = df
            elif seed > 0:
                df_all = pd.concat([df_all,df])#.append(df)
        pd.DataFrame.to_csv(df_all,output_file_csv)

def check_lcr():
    import workflow_analysis as wa
    import numpy as np
    record_lcr0 = np.zeros((11,))
    for i in range(11):
        at = wa.archimedean_tilings()
        at.generate_type_n(i+1)
        cross_lattice = np.cross(at.a1,at.a2)
        area_per_particle = cross_lattice[2]/len(at.position)
        area_hex = np.sqrt(3)/2.0
        lcr0 = np.sqrt(area_hex/area_per_particle)
        record_lcr0[i] = lcr0
        #print('type'+str(i+1)+': '+str(np.round(lcr0,4) ))
        del at
    return record_lcr0
        
def show_type_3_11_dual():
    import workflow_analysis as wa
    import numpy as np
    sdl = wa.show_dual_lattice()
    for i in np.linspace(3,11,9,dtype=int):
        sdl.show_dual_type_n_part(i,xylim=5)

def show_type_3_11_polygon_dye():
    import workflow_analysis as wa
    import numpy as np
    sdl = wa.archimedean_tilings_polygon_dye()
    for i in np.linspace(3,11,9,dtype=int):
        sdl.workflow_type_n(i,xylim=5,n_plus=3)

def compare_differen_gauss():
    pass