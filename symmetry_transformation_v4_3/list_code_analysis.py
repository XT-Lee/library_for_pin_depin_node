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
            s2g = ac.from_setup_to_a_gsd()
            s2g.set_file_parameters(simu_index,seed)
            
            g2d = ac.from_a_gsd_to_data()
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
                s2g = ac.from_setup_to_a_gsd()
                s2g.set_file_parameters(simu_index,seed)
                
                g2d = ac.from_a_gsd_to_data()
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

        
            