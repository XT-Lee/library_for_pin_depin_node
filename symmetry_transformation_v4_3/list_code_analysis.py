import symmetry_transformation_v4_3.analysis_controller as ac
import numpy as np
class analyze_a_gsd_file:
    def __init__(self):
        pass
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
            ss.get_gsd_data_from_filename()
            lf = ac.get_data_from_a_gsd_frame(ss.gsd_data[-1])
            lf.get_bonds_from_a_gsd_frame()
    
    def get_bond_plot_with_time_from_a_gsd(self,index1):
        ss = ac.get_a_gsd_from_setup()
        ss.set_file_parameters(index1,9)
        ss.get_gsd_data_from_filename()
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
        ggfs = ac.get_gsds_from_mysql()
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
        ggfs = ac.get_gsds_from_mysql()
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
    
    def get_cn3s_from_csv_honeycomb_part(self,output_file_csv):
        R"""
            lcr0.81 honeycomb pin
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        """prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
        index1 = 6013
        output_file_csv = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_'+str(int(index1))+'.csv'"""
        ggfs = ac.get_gsds_from_mysql()
        gsds = ggfs.get_record_from_csv(output_file_csv)
        n_simu = len(gsds)
        list_cn3 = np.zeros((n_simu,))
        #list_cn3[:,:2] =  record_np[:,2:4]

        for i in range(n_simu):
            frame = ggfs.get_last_frame_from_gsd(gsds[i])
            lcr1 = ggfs.record['lcr'].values[i]
            tune_dis = lcr1*3#list_cn3[i,0]*3
            gf = ac.get_data_from_a_gsd_frame(frame)#error:missing last frame
            cn3 = gf.get_cn3_from_a_gsd_frame(frame,tune_dis)
            list_cn3[i] = cn3
            del gf
        
        #save as csv
        import pandas as pd
        df = pd.DataFrame(ggfs.record,columns=[ 'simu_index','seed',"lcr",'trap_gauss_epsilon','temperature'])
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
        csv_filename = prefix_write + 'honeycomb_part_pin_lcr81_scan_k3_1000_fill_brownian.csv'
        pd.DataFrame.to_csv(df,csv_filename)
        print(df.head(35))
        return list_cn3

    def plot_from_csv(self):
        import pandas as pd
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/'
        csv_filename = prefix_write + 'honeycomb_pin_scan_k3_1000.csv'
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

        
            