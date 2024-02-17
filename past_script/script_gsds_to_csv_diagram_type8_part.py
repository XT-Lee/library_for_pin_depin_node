import time
tm1=time.localtime(time.time())

#after finishing simulation, we should calculate cnk,rho_trap_relative
prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/type_n_pin/"#
output_file_csv = prefix_write + "pin_hex_to_type_8_part_klt_2m_gauss_11083.csv"
cnk = 4
type_n = 8
part = True
import symmetry_transformation_v4_3.list_code_analysis as lca
asg = lca.analyze_a_series_of_gsd_file()
asg.get_cnks_from_csv_dipole_core_traps(cnk,output_file_csv)

import getDataAndDiagramCsv as gdadc
cdp = gdadc.csv_data_processor(csv_filename=output_file_csv)
cdp.get_relative_rho(type_n)#cdp.get_data_diagram(cnk)
cdp.save_csv(output_file_csv)
gdadc.get_diagram_from_csv_type_n(output_file_csv,type_n,part)

tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)

