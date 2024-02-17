import time
tm1=time.localtime(time.time())

prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/wca_yukawa/"
index1 = 7223
output_file_csv = prefix_write + "wca_yukawa_"+str(int(index1))+".csv"
"""import symmetry_transformation_v4_3.simulation_controller as sc
wy = sc.simulation_controller_wca_yukawa()
wy.generate_simu_index_csv_scan_u_phi()
wy.generate_initial_state_type2_particle_scan_csv(output_file_csv)"""

import symmetry_transformation_v4_3.list_code_analysis as lca
agf = lca.analyze_a_series_of_gsd_file()
coord_num_k = agf.get_coordination_number_k_by_given_type_n(1)
agf.get_cnks_from_csv_type_n_part_phi_core(coord_num_k,output_file_csv)
agf.extract_data_from_csv_u_yukawa_r1_phi_cn6(output_file_csv,"wca_yukawa_"+str(int(index1))+"_1_p.csv")

tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)