##!/usr/bin/env python
import time
import pandas as pd
tm1=time.localtime(time.time())


prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/"#
output_file_csv = prefix_write + "pin_hex_to_type_3_klt_2m_gauss_11513_11752.csv"


import symmetry_transformation_v4_3.list_code_analysis as lca
agf = lca.analyze_a_series_of_gsd_file()
#["simu_index","seed","lcr","trap_gauss_epsilon","temperature"]+['U_eq','cn3']
csv_filename = output_file_csv
df = pd.read_csv(csv_filename)
#df_sub = df[df['trap_gauss_epsilon'] == -30]
df_sub = df[df['seed'] == 0]
agf.get_bonds_from_simu_indices_list_type_n(
    df_sub['simu_index'].values,
    df_sub['seed'].values,
    df_sub['lcr'].values*3)

#agp = agf.extract_data_from_csv(output_file_csv,"dipole_"+str(int(index1))+"_p.csv")
import getDataAndDiagramCsv as gddc
cp = gddc.csv_data_processor(output_file_csv)
cp.get_relative_rho(3)
cp.save_csv(output_file_csv)
gddc.get_diagram_from_csv_type3(output_file_csv)

tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)
