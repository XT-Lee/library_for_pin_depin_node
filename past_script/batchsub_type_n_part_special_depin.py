##!/usr/bin/env python
import time
tm1=time.localtime(time.time())

#import symmetry_transformation_v4_3.list_code_simulation as lca
import symmetry_transformation_v4_3.simulation_controller as sc
hpt = sc.simulation_controller_type_n_part_traps()
#hpt.generate_simu_index_csv_depin_special_30()#
prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_depin/"
index1 = 12533
output_file_csv = prefix_write + "depin_type_n_from_type_n_part_special_klt_2m_gauss_"+str(int(index1))+".csv"
#hpt.generate_initial_state_type_n_particle_type_n_part_trap_depin_scan_csv_node(output_file_csv)#not checked

import symmetry_transformation_v4_3.list_code_analysis as lca
agf = lca.analyze_a_series_of_gsd_file()
#["simu_index","seed","lcr","trap_gauss_epsilon","temperature"]+['U_eq','cn3']
agf.get_bonds_from_simu_indices_type_n_mix_from_csv(output_file_csv)

tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)
