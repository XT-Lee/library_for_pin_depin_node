import symmetry_transformation_v4_3.analysis_controller as ac
import numpy as np
import os
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
            
            g2d = ac.get_data_from_a_gsd_frame_with_traps(s2g.gsd_data[-1])
            list_cn3[i] = g2d.get_cn_k_from_a_gsd_frame()
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
            
            g2d = ac.get_data_from_a_gsd_frame_with_traps(s2g.gsd_data[-1])
            list_cn3[i] = g2d.get_cn_k_from_a_gsd_frame()
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
                
                g2d = ac.get_data_from_a_gsd_frame_with_traps()
                list_cn3[i] = g2d.get_cn_k_from_a_gsd_frame(s2g.gsd_data[-1])
            
            
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
            
            g2d = ac.get_data_from_a_gsd_frame_with_traps()
            list_cn3[i] = g2d.get_cn_k_from_a_gsd_frame(s2g.gsd_data[-1])
            print(list_cn3[i])

    def get_bond_plot_from_a_gsd(self):
        ss = ac.get_a_gsd_from_setup()
        for index1 in np.linspace(5971,5974,4):#5950,5959
            ss.set_file_parameters(index1,9)
            ss.get_gsd_data_from_file()
            lf = ac.get_data_from_a_gsd_frame_with_traps(ss.gsd_data[-1])
            lf.get_bonds_from_a_gsd_frame()
    
    def get_bond_plot_with_time_from_a_gsd(self,index1):
        ss = ac.get_a_gsd_from_setup()
        ss.set_file_parameters(index1,9)
        ss.get_gsd_data_from_file()
        num_list = [100,500,1000,1500]#np.linspace(10,1999,10,dtype=int)
        for i in num_list:#range(2000-1):
            print(ss.gsd_data[i].configuration.box)
            #print(ss.gsd_data[i].configuration.box.dimensions)
            lf = ac.get_data_from_a_gsd_frame_with_traps(ss.gsd_data[i])
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
        ggfs = ac.get_gsds_from_csv()
        record = ggfs.get_record_from_sql_by_lcr()
        #record_np = np.array(record)
        gsds = ggfs.get_gsds_from_mysql_record(record)
        n_simu = len(gsds)
        list_cn3 = np.zeros((n_simu,))
        #list_cn3[:,:2] =  record_np[:,2:4]
        gf = ac.get_data_from_a_gsd_frame_with_traps()
        for i in range(n_simu):
            frame = ggfs.get_last_frame_from_gsd(gsds[i])
            tune_dis = record[i][2]*3#list_cn3[i,0]*3
            cn3 = gf.get_cn_k_from_a_gsd_frame(frame,tune_dis)
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
        ggfs = ac.get_gsds_from_csv()
        record = ggfs.get_record_from_sql_by_lcr(lcr1=0.80,table_name='pin_hex_to_honeycomb_part_klt_2m_gauss')
        #record_np = np.array(record)
        gsds = ggfs.get_gsds_from_mysql_record(record)
        n_simu = len(gsds)
        list_cn3 = np.zeros((n_simu,))
        #list_cn3[:,:2] =  record_np[:,2:4]

        for i in range(n_simu):
            frame = ggfs.get_last_frame_from_gsd(gsds[i])
            tune_dis = record[i][2]*3#list_cn3[i,0]*3
            gf = ac.get_data_from_a_gsd_frame_with_traps(frame)#error:missing last frame
            cn3 = gf.get_cn_k_from_a_gsd_frame(frame,tune_dis)
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
        ggfs = ac.get_gsds_from_csv()
        gsds = ggfs.get_record_from_csv(output_file_csv)
        n_simu = len(gsds)
        list_cn3 = np.zeros((n_simu,))
        #list_cn3[:,:2] =  record_np[:,2:4]

        for i in range(n_simu):
            frame = ggfs.get_last_frame_from_gsd(gsds[i])
            lcr1 = ggfs.record['lcr'].values[i]
            tune_dis = lcr1*3#list_cn3[i,0]*3
            gf = ac.get_data_from_a_gsd_frame_with_traps(frame)#error:missing last frame
            cn3 = gf.get_cn_k_from_a_gsd_frame(tune_dis)
            list_cn3[i] = cn3
            del gf
        
        #save as csv
        import pandas as pd
        df = pd.DataFrame(ggfs.record,columns=[ 'simu_index','seed',"lcr",'trap_gauss_epsilon','temperature'])
        df['cn3'] = list_cn3
        df['U_eq'] = df['trap_gauss_epsilon'].values*0.86466#*0.99613
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
    
    def get_cnks_from_csv_type_4569_part(self):
        R"""
        import symmetry_transformation_v4_3.list_code_analysis as lca
        asg = lca.analyze_a_series_of_gsd_file()
        asg.get_cnks_from_csv_type_4569_part()
        """
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/'

        #generate_simu_index_csv_type_4_pin_3_30()
        index1 = 6963#6823#[x]
        list_type_n = 4#[x]
        output_file_csv = prefix_write + 'pin_hex_to_type_'+str(int(list_type_n))+'_part_klt_2m_gauss_'+str(int(index1))+'.csv'#[x]
        self.get_cnks_from_csv_type_n_part_phi_core(3,output_file_csv)

        #generate_simu_index_csv_type_5_pin_3_30()
        index1 = 7013#6823#[x]
        list_type_n = 5#[x]
        output_file_csv = prefix_write + 'pin_hex_to_type_'+str(int(list_type_n))+'_part_klt_2m_gauss_'+str(int(index1))+'.csv'#[x]
        self.get_cnks_from_csv_type_n_part_phi_core(3,output_file_csv)

        #generate_simu_index_csv_type_6_pin_3_30()
        index1 = 7063#6823#[x]
        list_type_n = 6#[x]
        output_file_csv = prefix_write + 'pin_hex_to_type_'+str(int(list_type_n))+'_part_klt_2m_gauss_'+str(int(index1))+'.csv'#[x]
        self.get_cnks_from_csv_type_n_part_phi_core(3,output_file_csv)

        #generate_simu_index_csv_type_9_pin_3_30()
        index1 = 7113#6823#[x]
        list_type_n = 9#[x]
        output_file_csv = prefix_write + 'pin_hex_to_type_'+str(int(list_type_n))+'_part_klt_2m_gauss_'+str(int(index1))+'.csv'#[x]
        self.get_cnks_from_csv_type_n_part_phi_core(5,output_file_csv)

    def get_cnks_from_csv_type_n_part_phi_core(self,coord_num_k,output_csv_filename):
        R"""
        parameter:
            cn_k:(int) the k to calculate cn_k
        example:
            import symmetry_transformation_v4_3.list_code_analysis as lca
            asg = lca.analyze_a_series_of_gsd_file()
            asg.get_cn5s_from_csv_files()
            #["simu_index","seed","a","phi","u_dipole","temperature","type_n"]
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_csv_filename)
        lis_simu = gg.record['simu_index'].values
        list_as = gg.record['a'].values
        list_phis = gg.record['phi'].values#phi = pi/(2*sqrt(3)*a^2),a^2=pi/(2*sqrt(3)*phi)
        list_a_hexs = np.sqrt(np.pi/(list_phis*2*np.sqrt(3)))
        #lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        cnks = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            isExists=os.path.exists(gsds_filename[i])
            if isExists:
                file_size_b = os.path.getsize(gsds_filename[i])
                file_size_kb = file_size_b/1024
                if file_size_kb>100:
                    pdata.get_gsd_data_from_filename(gsds_filename[i])
                    gdata = ac.get_data_from_a_gsd_frame(pdata.gsd_data[-1])#get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
                    cnk = gdata.get_cn_k_from_a_gsd_frame(tune_dis=list_a_hexs[i],k=coord_num_k)#3*lcrs[i]
                    cnks[i] = cnk
                else:
                    cnks[i] = -1    
            else:
                cnks[i] = -1
        col_cnk = 'cn'+str(int(coord_num_k))
        gg.record[col_cnk] = cnks
        gg.record['a_hex'] = list_a_hexs
        r1 = np.exp(-0.25)
        gg.record['u_yukawa_r1'] =  gg.record['u_yukawa'].values*r1
        #gg.record['U_eq'] =  gg.record['trap_gauss_epsilon'].values*0.86466#*0.99613
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_csv_filename)
    
    def get_cnks_from_csv_type_mix_core(self,coord_num_k,output_csv_filename):
        R"""
        parameter:
            cn_k:(int) the k to calculate cn_k
        example:
            import symmetry_transformation_v4_3.list_code_analysis as lca
            asg = lca.analyze_a_series_of_gsd_file()
            asg.get_cn5s_from_csv_files()
            #["simu_index","seed","a","phi","u_dipole","temperature","type_n"]
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_csv_filename)
        if coord_num_k is None:
            list_type_n = gg.record['type_n'].values
            list_type_n = np.array(list_type_n,dtype=int)
        
        #lis_simu = gg.record['simu_index'].values
        #list_as = gg.record['a'].values
        #list_phis = gg.record['phi'].values#phi = pi/(2*sqrt(3)*a^2),a^2=pi/(2*sqrt(3)*phi)
        #list_a_hexs = np.sqrt(np.pi/(list_phis*2*np.sqrt(3)))
        lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        cnks = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            isExists=os.path.exists(gsds_filename[i])
            if isExists:
                file_size_b = os.path.getsize(gsds_filename[i])
                file_size_kb = file_size_b/1024
                if file_size_kb>1000:
                    pdata.get_gsd_data_from_filename(gsds_filename[i])
                    gdata = ac.get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])#get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
                    cnk = gdata.get_cn_k_from_a_gsd_frame(tune_dis=3*lcrs[i],k=coord_num_k)#
                    cnks[i] = cnk
                else:
                    cnks[i] = -1    
            else:
                cnks[i] = -1
        col_cnk = 'cnk'#+str(int(coord_num_k))
        gg.record[col_cnk] = cnks
        #gg.record['a_hex'] = list_a_hexs
        #r1 = np.exp(-0.25)
        #gg.record['u_yukawa_r1'] =  gg.record['u_yukawa'].values*r1
        gg.record['U_eq'] =  gg.record['trap_gauss_epsilon'].values*0.86466#*0.99613
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_csv_filename)

    def get_cnks_from_csv_dipole_core_traps(self,coord_num_k,output_csv_filename):
        R"""
        parameter:
            cn_k:(int) the k to calculate cn_k
        example:
            import symmetry_transformation_v4_3.list_code_analysis as lca
            asg = lca.analyze_a_series_of_gsd_file()
            asg.get_cn5s_from_csv_files()
            #["simu_index","seed","a","phi","u_dipole","temperature","type_n"]
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_csv_filename)
        #lis_simu = gg.record['simu_index'].values
        #list_as = gg.record['a'].values
        #list_phis = gg.record['phi'].values#phi = pi/(2*sqrt(3)*a^2),a^2=pi/(2*sqrt(3)*phi)
        #list_a_hexs = np.sqrt(np.pi/(list_phis*2*np.sqrt(3)))
        lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        cnks = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            isExists=os.path.exists(gsds_filename[i])
            if isExists:
                file_size_b = os.path.getsize(gsds_filename[i])
                file_size_kb = file_size_b/1024
                if file_size_kb>1000:
                    pdata.get_gsd_data_from_filename(gsds_filename[i])
                    gdata = ac.get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])#get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
                    cnk = gdata.get_cn_k_from_a_gsd_frame(tune_dis=3*lcrs[i],k=coord_num_k)#
                    cnks[i] = cnk
                else:
                    cnks[i] = -1    
            else:
                cnks[i] = -1
        col_cnk = 'cn'+str(int(coord_num_k))
        gg.record[col_cnk] = cnks
        #gg.record['a_hex'] = list_a_hexs
        #r1 = np.exp(-0.25)
        #gg.record['u_yukawa_r1'] =  gg.record['u_yukawa'].values*r1
        gg.record['U_eq'] =  gg.record['trap_gauss_epsilon'].values*0.86466#*0.99613
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_csv_filename)

    def get_cnks_from_csv_dipole_core(self,coord_num_k,output_csv_filename):
        R"""
        parameter:
            cn_k:(int) the k to calculate cn_k
        example:
            import symmetry_transformation_v4_3.list_code_analysis as lca
            asg = lca.analyze_a_series_of_gsd_file()
            asg.get_cn5s_from_csv_files()
            #["simu_index","seed","a","phi","u_dipole","temperature","type_n"]
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_csv_filename)
        #lis_simu = gg.record['simu_index'].values
        #list_as = gg.record['a'].values
        #list_phis = gg.record['phi'].values#phi = pi/(2*sqrt(3)*a^2),a^2=pi/(2*sqrt(3)*phi)
        #list_a_hexs = np.sqrt(np.pi/(list_phis*2*np.sqrt(3)))
        lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        cnks = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            isExists=os.path.exists(gsds_filename[i])
            if isExists:
                file_size_b = os.path.getsize(gsds_filename[i])
                file_size_kb = file_size_b/1024
                if file_size_kb>1000:
                    pdata.get_gsd_data_from_filename(gsds_filename[i])
                    gdata = ac.get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])#get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
                    cnk = gdata.get_cn_k_from_a_gsd_frame(tune_dis=3*lcrs[i],k=coord_num_k)#
                    cnks[i] = cnk
                else:
                    cnks[i] = -1    
            else:
                cnks[i] = -1
        col_cnk = 'cn'+str(int(coord_num_k))
        gg.record[col_cnk] = cnks
        #gg.record['a_hex'] = list_a_hexs
        #r1 = np.exp(-0.25)
        #gg.record['u_yukawa_r1'] =  gg.record['u_yukawa'].values*r1
        gg.record['U_eq'] =  gg.record['trap_gauss_epsilon'].values*0.86466#*0.99613
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_csv_filename)

    def extract_data_from_csv_u_dipole_phi_cn6(self,input_csv_filename,output_csv_filename):
        R"""
        import symmetry_transformation_v4_3.list_code_analysis as lca
        agf = lca.analyze_a_series_of_gsd_file()
        agp = agf.extract_data_from_csv(output_file_csv,"wca_dipole_"+str(int(index1))+"_p.csv")
        """
        import pandas as pd
        record = pd.read_csv(input_csv_filename)
        record_seed0 = record[record['seed']==0]
        sub_record = record_seed0[['u_dipole','phi','cn6']]
        pd.DataFrame.to_csv(sub_record,output_csv_filename)
    
    def extract_data_from_csv_u_yukawa_r1_phi_cn6(self,input_csv_filename,output_csv_filename):
        R"""
        import symmetry_transformation_v4_3.list_code_analysis as lca
        agf = lca.analyze_a_series_of_gsd_file()
        agp = agf.extract_data_from_csv(output_file_csv,"wca_dipole_"+str(int(index1))+"_p.csv")
        """
        import pandas as pd
        record = pd.read_csv(input_csv_filename)
        record_seed0 = record[record['seed']==0]
        sub_record = record_seed0[['u_yukawa_r1','phi','cn6']]
        pd.DataFrame.to_csv(sub_record,output_csv_filename)
    
    def extract_data_from_csv_u_yukawa_r1_U_eq_cnk(self,input_csv_filename,output_csv_filename,coord_num_k):
        R"""
        import symmetry_transformation_v4_3.list_code_analysis as lca
        agf = lca.analyze_a_series_of_gsd_file()
        agp = agf.extract_data_from_csv(output_file_csv,"wca_dipole_"+str(int(index1))+"_p.csv")
        """
        import pandas as pd
        record = pd.read_csv(input_csv_filename)
        record_seed0 = record[record['seed']==0]
        sub_record = record_seed0[['u_yukawa_r1','U_eq','cn'+str(coord_num_k)]]
        pd.DataFrame.to_csv(sub_record,output_csv_filename)
    
    def get_grs_from_csv_type_n_part(self,input_csv_filename):
        import symmetry_transformation_v4_3.analysis_controller as ac
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(input_csv_filename)
        lis_simu = gg.record['simu_index'].values
        list_seeds = gg.record['seed'].values
        list_as = gg.record['a'].values
        list_phis = gg.record['phi'].values#phi = pi/(2*sqrt(3)*a^2),a^2=pi/(2*sqrt(3)*phi)
        list_a_hexs = gg.record['a_hex'].values
        #lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        for i in range(len(gsds_filename)):
            if (list_seeds[i]==0) and (lis_simu[i]>=7753):
                pdata.get_gsd_data_from_filename(gsds_filename[i])
                gdata = ac.get_data_from_a_gsd_frame(pdata.gsd_data[-1])#get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
                gdata.get_grs_from_a_gsd_frame(list_a_hexs[i],str(lis_simu[i])+"_"+str(list_a_hexs[i])+".png",list_a_hexs[i])

   
    def get_cnks_from_csv_type_n_part_wca_yukawa_depin(self,coord_num_k,output_file_csv):
        R"""
        parameter:
            cn_k:(int) the k to calculate cn_k
        example:
            import symmetry_transformation_v4_3.list_code_analysis as lca
            asg = lca.analyze_a_series_of_gsd_file()
            asg.get_cn5s_from_csv_files()
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_file_csv)
        lis_simu = gg.record['simu_index'].values
        lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        cnks = np.zeros((len(gsds_filename),))
        cn6s = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            isExists=os.path.exists(gsds_filename[i])
            if isExists:
                #analyze gsd file only which have been written completely 
                file_size_b = os.path.getsize(gsds_filename[i])
                file_size_kb = file_size_b/1024
                if file_size_kb>1000:
                    pdata.get_gsd_data_from_filename(gsds_filename[i])
                    gdata = ac.get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
                    cnk = gdata.get_cn_k_from_a_gsd_frame(tune_dis=3*lcrs[i])#,k=coord_num_k
                    cnks[i] = cnk[coord_num_k]
                    cn6s[i] = cnk[6]
                else:
                    #print(gsds_filename[i])
                    #print(file_size_kb)
                    cnks[i] = -1   
                    cn6s[i] = -1
            else:
                cnks[i] = -1
                cn6s[i] = -1
        col_cnk = 'cn'+str(int(coord_num_k))
        gg.record[col_cnk] = cnks
        gg.record['cn6'] = cn6s
        gg.record['U_eq'] =  gg.record['gauss_epsilon'].values*(-0.99613)
        r1 = np.exp(-0.25)
        gg.record['u_yukawa_r1'] =  gg.record['u_yukawa'].values*r1
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_file_csv)
    
    def get_coordination_number_k_by_given_type_n(self,type_n):
        import workflow_analysis as wa
        at = wa.archimedean_tilings()
        coord_num_k = at.get_coordination_number_k_for_type_n(type_n)
        return coord_num_k
        
    def get_cn5s_from_csv_files_type_11_part(self):
        R"""
            import symmetry_transformation_v4_3.list_code_analysis as lca
            asg = lca.analyze_a_series_of_gsd_file()
            asg.get_cn5s_from_csv_files()
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/'
        index1 = 6733#index1 = 6643
        output_file_csv = prefix_write + 'pin_hex_to_type_11_part_klt_2m_gauss_'+str(int(index1))+'.csv'
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_file_csv)
        lis_simu = gg.record['simu_index'].values
        lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        cn5s = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            if lis_simu[i] > 6734:#lis_simu[i] < 6798 and 
                pdata.get_gsd_data_from_filename(gsds_filename[i])
                gdata = ac.get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
                cn5 = gdata.get_cn_k_from_a_gsd_frame(tune_dis=3*lcrs[i],k=5)
                cn5s[i] = cn5
            else:
                cn5s[i] = -1
        gg.record['cn5'] = cn5s
        gg.record['U_eq'] =  gg.record['trap_gauss_epsilon'].values*0.86466#*0.99613
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_file_csv)
    
    def get_cn4s_from_csv_files_type_7_part(self):
        R"""
            import symmetry_transformation_v4_3.list_code_analysis as lca
            asg = lca.analyze_a_series_of_gsd_file()
            asg.get_cn5s_from_csv_files()
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/'
        index1 = 6913#6823#[x]
        output_file_csv = prefix_write + 'pin_hex_to_type_7_part_klt_2m_gauss_'+str(int(index1))+'.csv'#[x]
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_file_csv)
        lis_simu = gg.record['simu_index'].values
        lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        cn4s = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            pdata.get_gsd_data_from_filename(gsds_filename[i])
            gdata = ac.get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
            cn4 = gdata.get_cn_k_from_a_gsd_frame(tune_dis=3*lcrs[i],k=4)
            cn4s[i] = cn4

        gg.record['cn4'] = cn4s
        gg.record['U_eq'] =  gg.record['trap_gauss_epsilon'].values*0.86466#*0.99613
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_file_csv)

    def get_cn4s_from_csv_files_type_8(self):
        R"""
            import symmetry_transformation_v4_3.list_code_analysis as lca
            asg = lca.analyze_a_series_of_gsd_file()
            asg.get_cn4s_from_csv_files_type_8()
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/'
        index1 = 243#[x]
        output_file_csv = prefix_write + 'pin_hex_to_type_8_klt_2m_gauss_'+str(int(index1))+'.csv'#[x]
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_file_csv)
        lis_simu = gg.record['simu_index'].values
        list_seed = gg.record['seed'].values
        lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        cn4s = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            if list_seed[i]==0:
                pdata.get_gsd_data_from_filename(gsds_filename[i])
                gdata = ac.get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
                cn4 = gdata.get_cn_k_from_a_gsd_frame(tune_dis=3*lcrs[i],k=4)
                cn4s[i] = cn4
            else:
                cn4s[i] = -1
        gg.record['cn4'] = cn4s
        gg.record['U_eq'] =  gg.record['trap_gauss_epsilon'].values*0.86466#*0.99613
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_file_csv)
    
    def get_cn4s_from_csv_files_type_8_part(self):
        R"""
            import symmetry_transformation_v4_3.list_code_analysis as lca
            asg = lca.analyze_a_series_of_gsd_file()
            asg.get_cn4s_from_csv_files_type_8()
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/'
        index1 = 513#[x]
        output_file_csv = prefix_write + 'pin_hex_to_type_8_klt_2m_gauss_'+str(int(index1))+'.csv'#[x]
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_file_csv)
        lis_simu = gg.record['simu_index'].values
        list_seed = gg.record['seed'].values
        lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        cn4s = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            if list_seed[i]==0:
                pdata.get_gsd_data_from_filename(gsds_filename[i])
                gdata = ac.get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
                cn4 = gdata.get_cn_k_from_a_gsd_frame(tune_dis=3*lcrs[i],k=4)
                cn4s[i] = cn4
            else:
                cn4s[i] = -1
        gg.record['cn4'] = cn4s
        gg.record['U_eq'] =  gg.record['trap_gauss_epsilon'].values*0.86466#*0.99613
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_file_csv)

    def get_cn3s_from_csv_files_type_3(self):
        R"""
        import symmetry_transformation_v4_3.list_code_analysis as lca
        agf = lca.analyze_a_series_of_gsd_file()
        agp = agf.get_cn3s_from_csv_files_type_3()
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_write = '/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/'
        output_file_csv = prefix_write + 'pin_hex_to_honeycomb_klt_2m_gauss_3_242.csv'
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_file_csv)
        list_seed = gg.record['seed'].values
        lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        cn3s = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            if list_seed[i]==9:
                pdata.get_gsd_data_from_filename(gsds_filename[i])
                gdata = ac.get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
                cn3 = gdata.get_cn_k_from_a_gsd_frame(tune_dis=3*lcrs[i])
                cn3s[i] = cn3
            else:
                cn3s[i] = -1
        gg.record['cn3'] = cn3s
        gg.record['U_eq'] =  gg.record['trap_gauss_epsilon'].values*0.86466#*0.99613
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_file_csv)
    
    def get_cn6s_from_csv_files_type_2(self):
        R"""
        introduction:
            column = ["simu_index","seed","a","phi","u_yukawa","temperature","type_n"]
        example:
            import symmetry_transformation_v4_3.list_code_analysis as lca
            agf = lca.analyze_a_series_of_gsd_file()
            agp = agf.get_cn6s_from_csv_files_type_2()
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_write = '/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/wca_yukawa/'
        index1 = 7223
        output_file_csv = prefix_write + "wca_yukawa_"+str(int(index1))+".csv"
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_file_csv)
        list_seed = gg.record['seed'].values
        list_simu_index = gg.record['simu_index'].values
        list_as = gg.record['a'].values
        pdata = ac.get_a_gsd_from_setup()
        cn6s = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            if (list_seed[i]==0)and(list_simu_index[i]>=7234) :#7278
                pdata.get_gsd_data_from_filename(gsds_filename[i])
                gdata = ac.get_data_from_a_gsd_frame(pdata.gsd_data[-1])
                cn6 = gdata.get_cn_k_from_a_gsd_frame(tune_dis=list_as[i],k=6)
                cn6s[i] = cn6
            else:
                cn6s[i] = -1
        gg.record['cn6'] = cn6s
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_file_csv)

    def get_cn3s_from_csv_files_type_3_part(self):
        R"""
            import symmetry_transformation_v4_3.list_code_analysis as lca
            asg = lca.analyze_a_series_of_gsd_file()
            asg.get_cn3s_from_csv_files_type_3_part()
        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
        output_file_csv = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_6373_6612.csv'
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_file_csv)
        lis_simu = gg.record['simu_index'].values
        list_seed = gg.record['seed'].values
        lcrs = gg.record['lcr'].values
        pdata = ac.get_a_gsd_from_setup()
        cn3s = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            if list_seed[i]<1 and lis_simu[i]>6373:
                pdata.get_gsd_data_from_filename(gsds_filename[i])
                gdata = ac.get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
                cn3 = gdata.get_cn_k_from_a_gsd_frame(tune_dis=3*lcrs[i],k=3)
                cn3s[i] = cn3
            else:
                cn3s[i] = -1
        gg.record['cn3'] = cn3s
        gg.record['U_eq'] =  gg.record['trap_gauss_epsilon'].values*0.86466#*0.99613
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_file_csv)
        
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
        gg = ac.get_gsds_from_csv()
        gsds_filename = gg.get_record_from_csv(output_file_csv3)
        """lcrs = gg.record['lcr'].values
        ks = gg.record['trap_gauss_epsilon'].values"""
        pdata = ac.get_a_gsd_from_setup()
        cn3s = np.zeros((len(gsds_filename),))
        for i in range(len(gsds_filename)):
            pdata.get_gsd_data_from_file(gsds_filename[i])
            gdata = ac.get_data_from_a_gsd_frame_with_traps(pdata.gsd_data[-1])
            cn3 = gdata.get_cn_k_from_a_gsd_frame()
            cn3s[i] = cn3
        gg.record['cn3'] = cn3s
        import pandas as pd
        list_col = gg.record.columns.values[1:]#to strip 'unnamed:0' column which contains [1,2,3...]
        pd.DataFrame.to_csv(gg.record[list_col],output_file_csv3)

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
    
    def add_cn3_csv_from_csv(self,csv_filename):
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
        #prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/'
        index0 = 3#10362#7223#6302
        for i in range(10):
            index1 = index0+i#658#
            #gsd_file = prefix_write + 'trajectory_auto'+str(int(index1))+'_9.gsd'
            ggsd = ac.get_a_gsd_from_setup()
            ggsd.set_file_parameters(index1,7)
            ggsd.get_gsd_data_from_file()
            frame = ggsd.gsd_data[-1]
            
            gf = ac.get_data_from_a_gsd_frame_with_traps(frame)#error:missing last frame
            """if i<5:
                gf.get_bonds_png_from_a_gsd_frame('depin_from_type_'+str(6)+'_part_special_'+str(index1)+'.png')#4+i
            else:
                gf.get_bonds_png_from_a_gsd_frame('depin_from_type_'+str(7)+'_part_special_'+str(index1)+'.png')#4+i
            del gf"""
            gf.get_bonds_png_from_a_gsd_frame(str(index1)+'_0_final.png')
            
    def get_bonds_from_simu_indices_list_type_n(self, list_index, list_seed=9, list_lcra=2.4):
        R"""

        """
        import symmetry_transformation_v4_3.analysis_controller as ac
        # list_index = [7001, 7052, 7112, 7137, 7141]
        for i in range(len(list_index)):  # range(5):
            index1 = list_index[i]
            seed1 = list_seed[i]
            lcra1 = list_lcra[i]
            # gsd_file = prefix_write + 'trajectory_auto'+str(int(index1))+'_9.gsd'
            ggsd = ac.get_a_gsd_from_setup()
            ggsd.set_file_parameters(index1, seed1)
            ggsd.get_gsd_data_from_file()
            frame = ggsd.gsd_data[-1]

            gf = ac.get_data_from_a_gsd_frame_with_traps(frame)  # error:missing last frame
            gf.get_bonds_png_from_a_gsd_frame('bond_'+str(index1)+'_'+str(seed1)+'.png', lcra1)
            del gf
    
    def get_bonds_from_simu_indices_type_n_from_csv(self,filename_csv):
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
        #prefix_csv = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_depin/'#type_n_depin/
        #filename_csv = prefix_csv + 'depin_type_n_from_type_n_part_klt_2m_gauss_7163.csv'
        #'pin_hex_to_type_n_part_klt_2m_gauss_6613.csv'#
        #depin_type_n_from_type_n_part_klt_2m_gauss_7163,6053,6293 
        record_csv = pd.read_csv(filename_csv)
        list_type_n_to_watch = [4,5,6,9]#[7,10,11]
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
                gf = ac.get_data_from_a_gsd_frame_with_traps(frame)#error:missing last frame
                #gf.get_given_bonds_png_from_a_gsd_frame('pin_hex_to_type_type_'+str(type_n)+'_part_'+str(-int(list_trap_gauss_epsilon[i]))+'_.png',3*lcr0*1.2)
                gf.get_given_bonds_png_from_a_gsd_frame('depin_from_type_'+str(type_n)+'_part_'+str(-int(list_trap_gauss_epsilon[i]))+'_.png',3*lcr0*1.2)
                del gf
    
    def get_bonds_from_simu_indices_type_n_mix_from_csv(self,filename_csv):
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
        record_csv = pd.read_csv(filename_csv)
        list_type_n = record_csv['type_n'].values
        list_seed = record_csv['seed'].values
        list_index = record_csv['simu_index'].values
        list_lcr0 = record_csv['lcr'].values
        list_trap_gauss_epsilon = record_csv['trap_gauss_epsilon'].values
        ggsd = ac.get_a_gsd_from_setup()
        for i in range(len(list_index)):
            ggsd.set_file_parameters(list_index[i],list_seed[i])
            frame = ggsd.gsd_data[-1]
            gf = ac.get_data_from_a_gsd_frame_with_traps(frame)#error:missing last frame
            gf.get_given_bonds_png_from_a_gsd_frame('pin_hex_to_type_type_'+str(list_type_n[i])+'_part_'+str(-int(list_trap_gauss_epsilon[i]))+'_.png',3*list_lcr0[i]*1.3)
            #gf.get_given_bonds_png_from_a_gsd_frame('depin_from_type_'+str(list_type_n[i])+'_part_'+str(-int(list_trap_gauss_epsilon[i]))+'_.png',3*list_lcr0[i]*1.2)
            del gf

    def get_bonds_from_simu_indices_type_11_from_csv(self):
        R"""
            column = ['simu_index','seed','lcr','trap_gauss_epsilon','temperature','type_n']
            manually set bond_length_max = 3*lcr0*1.2, set 1.2 is to avoid 1,414 and 1.73 bond.
        exp:
            import symmetry_transformation_v4_3.list_code_analysis as lca
            agf = lca.analyze_a_series_of_gsd_file()
            agp = agf.get_bonds_from_simu_indices_type_10_from_csv()
        """
        import pandas as pd
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/'
        index1 = 6733#index1 = 6643
        output_file_csv = prefix_write + 'pin_hex_to_type_11_part_klt_2m_gauss_'+str(int(index1))+'.csv'
        #depin_type_n_from_type_n_part_klt_2m_gauss_6053,6293
        record_csv = pd.read_csv(output_file_csv)
        ggsd = ac.get_a_gsd_from_setup()
        record_type_n = record_csv[record_csv['cn5'] > 0.2]
        list_index = record_type_n['simu_index'].values
        lcr0 = record_type_n['lcr'].values#[0]
        list_trap_gauss_epsilon = record_type_n['trap_gauss_epsilon'].values
        for i in range(len(list_index)):
            #gsd_file = prefix_write + 'trajectory_auto'+str(int(index1))+'_9.gsd'
            ggsd.set_file_parameters(list_index[i],9)
            #ggsd.get_gsd_data_from_file()
            frame = ggsd.gsd_data[-1]
            gf = ac.get_data_from_a_gsd_frame_with_traps(frame)#error:missing last frame
            gf.get_given_bonds_png_from_a_gsd_frame('pin_hex_to_type_type_11_part_'+str(int(list_index[i]))+'_9_.png',3*lcr0[i]*1.2)
            #gf.get_given_bonds_png_from_a_gsd_frame('depin_from_type_'+str(type_n)+'_part_'+str(-int(list_trap_gauss_epsilon[i]))+'_.png',3*lcr0*1.2)
            del gf

    def get_bonds_from_simu_indices_type_10_from_csv(self):
        R"""
            column = ['simu_index','seed','lcr','trap_gauss_epsilon','temperature','type_n']
            manually set bond_length_max = 3*lcr0*1.2, set 1.2 is to avoid 1,414 and 1.73 bond.
        exp:
            import symmetry_transformation_v4_3.list_code_analysis as lca
            agf = lca.analyze_a_series_of_gsd_file()
            agp = agf.get_bonds_from_simu_indices_type_10_from_csv()
        """
        import pandas as pd
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/'
        index1 = 6643#index1 = 6733#
        output_file_csv = prefix_write + 'pin_hex_to_type_10_part_klt_2m_gauss_'+str(int(index1))+'.csv'
        #depin_type_n_from_type_n_part_klt_2m_gauss_6053,6293
        record_csv = pd.read_csv(output_file_csv)
        ggsd = ac.get_a_gsd_from_setup()
        record_type_n = record_csv[record_csv['cn5'] > 0.8]
        list_index = record_type_n['simu_index'].values
        lcr0 = record_type_n['lcr'].values#[0]
        list_trap_gauss_epsilon = record_type_n['trap_gauss_epsilon'].values
        for i in range(len(list_index)):
            #gsd_file = prefix_write + 'trajectory_auto'+str(int(index1))+'_9.gsd'
            ggsd.set_file_parameters(list_index[i],9)
            #ggsd.get_gsd_data_from_file()
            frame = ggsd.gsd_data[-1]
            gf = ac.get_data_from_a_gsd_frame_with_traps(frame)#error:missing last frame
            gf.get_given_bonds_png_from_a_gsd_frame('pin_hex_to_type_type_10_part_'+str(int(list_index[i]))+'_9_.png',3*lcr0[i]*1.2)
            #gf.get_given_bonds_png_from_a_gsd_frame('depin_from_type_'+str(type_n)+'_part_'+str(-int(list_trap_gauss_epsilon[i]))+'_.png',3*lcr0*1.2)
            del gf

    def get_bonds_from_simu_indices_type_7_from_csv(self):
        R"""
            column = ['simu_index','seed','lcr','trap_gauss_epsilon','temperature','type_n']
            manually set bond_length_max = 3*lcr0*1.2, set 1.2 is to avoid 1,414 and 1.73 bond.
        exp:
            import symmetry_transformation_v4_3.list_code_analysis as lca
            agf = lca.analyze_a_series_of_gsd_file()
            agp = agf.get_bonds_from_simu_indices_type_7_from_csv()
        """
        import pandas as pd
        import symmetry_transformation_v4_3.analysis_controller as ac
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/'
        #index1 = 6643#index1 = 6733#
        output_file_csv = prefix_write + 'pin_hex_to_type_7_part_klt_2m_gauss_merge.csv'#+str(int(index1))+
        #depin_type_n_from_type_n_part_klt_2m_gauss_6053,6293
        record_csv = pd.read_csv(output_file_csv)
        ggsd = ac.get_a_gsd_from_setup()
        record_type_n = record_csv[record_csv['cn4'] > 0.74]
        list_index = record_type_n['simu_index'].values
        lcr0 = record_type_n['lcr'].values#[0]
        list_trap_gauss_epsilon = record_type_n['trap_gauss_epsilon'].values
        for i in range(len(list_index)):
            #gsd_file = prefix_write + 'trajectory_auto'+str(int(index1))+'_9.gsd'
            ggsd.set_file_parameters(list_index[i],9)
            #ggsd.get_gsd_data_from_file()
            frame = ggsd.gsd_data[-1]
            gf = ac.get_data_from_a_gsd_frame_with_traps(frame)#error:missing last frame
            gf.get_given_bonds_png_from_a_gsd_frame('pin_hex_to_type_7_part_'+str(int(list_index[i]))+'_9_.png',3*lcr0[i]*1.2)
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
                
                g2d = ac.get_data_from_a_gsd_frame_with_traps()
                list_cn3[i] = g2d.get_cn_k_from_a_gsd_frame(s2g.gsd_data[-1])
            
            
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

    #def get_cn3s_diagram_from_csv_files

class analyze_a_csv_file:
    def __init__(self):
        pass

    def analyze_pin_node_20240213(self):
        R"""
        | simu_index | seed | lcr  | trap_gauss_epsilon | temperature | type_n | cnk | U_eq |
        """
        prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/"#
        output_file_csv = prefix_write + "pin_hex_to_type_3_klt_2m_gauss_10333_10572.csv"
        import getDataAndDiagramCsv as gddc
        csvd = gddc.csv_data_processor(output_file_csv)
        csvd.select_single_seed(0)
        lcrs = csvd.sub_record['lcr'].values
        us = csvd.sub_record['U_eq'].values
        cn4s = csvd.sub_record['cn4s'].values
        diagramp = gddc.diagram_processor()
        diagramp.draw_diagram_scatter_oop()
        pass

def check_lcr():
    import workflow_analysis as wa
    at = wa.archimedean_tilings()
    record_lcr0 = at.get_type_n_lcr0()
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