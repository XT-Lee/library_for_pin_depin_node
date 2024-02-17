import time
tm1=time.localtime(time.time())

#read all the codes written in lxt_code_py16_node,
#write codes in bachsub_wca_dipole_pin.py,
#to operate the simulation similar to the codes in past_script/batchsub_wca_yukawa_depin_fix_phi.py

import symmetry_transformation_v4_3.simulation_controller as sc
wy = sc.simulation_controller_wca_dipole()
wy.generate_simu_index_csv_scan_u_phi()
prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/wca_dipole/"#
index1 = 7683
output_file_csv = prefix_write + "wca_dipole_"+str(int(index1))+".csv"
wy.generate_initial_state_type2_particle_scan_csv(output_file_csv)#n_size=16, temporarily

"""import symmetry_transformation_v4_3.list_code_analysis as lca
agf = lca.analyze_a_series_of_gsd_file()
#["simu_index","seed","a","phi","u_dipole","temperature","type_n"]
agp = agf.get_cnks_from_csv_type_n_part_core(6,output_file_csv)
agp = agf.extract_data_from_csv(output_file_csv,"wca_dipole_"+str(int(index1))+"_p.csv")
"""
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)