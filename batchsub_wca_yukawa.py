import time
tm1=time.localtime(time.time())


import symmetry_transformation_v4_3.simulation_controller as sc
wy = sc.simulation_controller_wca_yukawa()
#wy.generate_simu_index_csv_scan_u_phi()
prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/wca_yukawa/"
index1 = 7223
output_file_csv = prefix_write + "wca_yukawa_"+str(int(index1))+".csv"
wy.generate_initial_state_type2_particle_scan_csv(output_file_csv)

tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)