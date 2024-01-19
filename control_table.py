#import threading
import time
import computeTime as ct
tm1=time.localtime(time.time())

output_file_csv = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/pin_hex_to_honeycomb_klt_2m_gauss_3_242.csv"#0-6
            
import symmetry_transformation_v4_3.list_code_analysis as lca
asg = lca.analyze_a_series_of_gsd_file()

asg.get_cnks_from_csv_type_n_part_core(3,output_file_csv)


tm2=time.localtime(time.time())
#calculate the time cost
ct.getTimeCost(tm1,tm2)

