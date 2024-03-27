##!/usr/bin/env python
import time
tm1=time.localtime(time.time())

#import symmetry_transformation_v4_3.list_code_simulation as lca
import symmetry_transformation_v4_3.simulation_controller as sc
hpt = sc.simulation_controller_kagome_part_traps()
hpt.generate_simu_index_csv_kagome_part_pin_12263_12532()#
prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/"#
output_file_csv = prefix_write + "pin_hex_to_type_8_part_klt_2m_gauss_12263.csv"
hpt.generate_initial_state_hex_particle_kagome_part_trap_scan_csv(output_file_csv)

"""import symmetry_transformation_v4_3.list_code_analysis as lca
agf = lca.analyze_a_series_of_gsd_file()
#["simu_index","seed","lcr","trap_gauss_epsilon","temperature"]+['U_eq','cn3']
agp = agf.get_cnks_from_csv_dipole_core(4,output_file_csv)
#agp = agf.extract_data_from_csv(output_file_csv,"dipole_"+str(int(index1))+"_p.csv")
import getDataAndDiagramCsv as gddc
cp = gddc.csv_data_processor(output_file_csv)
cp.get_relative_rho(8)
cp.save_csv(output_file_csv)
gddc.get_diagram_from_csv_type_n(output_file_csv,8,True)"""

tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)
