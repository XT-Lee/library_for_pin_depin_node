##!/usr/bin/env python
import time
tm1=time.localtime(time.time())

#import symmetry_transformation_v4_3.list_code_simulation as lca
import symmetry_transformation_v4_3.simulation_controller as sc
hpt = sc.simulation_controller_kagome_traps()
#hpt.generate_simu_index_csv_kagome_pin_10813_11082()
prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/"#
output_file_csv = prefix_write + "pin_hex_to_type_8_klt_2m_gauss_10813.csv"
hpt.generate_initial_state_hex_particle_kagome_trap_scan_csv(output_file_csv)

tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)
