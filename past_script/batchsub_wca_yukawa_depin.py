import time
tm1=time.localtime(time.time())

import symmetry_transformation_v4_3.list_code_analysis as lca
agf = lca.analyze_a_series_of_gsd_file()

type_ns = [8,10,11]#3,72,
index1s = [7333+70*2,7333+70*3,7333+70*4]#7333,7333+70,
"""import symmetry_transformation_v4_3.simulation_controller as sc
wy = sc.simulation_controller_wca_yukawa_traps()
wy.get_lcr0()
print(wy.lcr0)
for i in range(5):
    type_n = type_ns[i]
    index1 = index1s[i]
    wy.generate_simu_index_csv_depin_near_melt(index1,type_n)"""
for i in len(type_ns):
    type_n = type_ns[i]
    index1 = index1s[i]
    prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/wca_yukawa/"
    output_file_csv = prefix_write + "wca_yukawa_depin_from_type"+str(type_n)+"_"+str(int(index1))+".csv"
    #wy.generate_initial_state_type_n_particle_type_n_part_trap_scan_csv(output_file_csv)
    coord_num_k = agf.get_coordination_number_k_by_given_type_n(type_n)
    agf.get_cnks_from_csv_type_n_part_wca_yukawa_depin(coord_num_k,output_file_csv)
    agf.extract_data_from_csv_u_yukawa_r1_U_eq_cnk(output_file_csv,"wca_yukawa_depin_from_type"+str(type_n)+"_"+str(int(index1))+"_p.csv",coord_num_k)

tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)