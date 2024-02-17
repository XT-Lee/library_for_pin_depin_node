import time
tm1=time.localtime(time.time())

import symmetry_transformation_v4_3.list_code_analysis as ac
aa = ac.analyze_a_series_of_gsd_file()
aa.get_bonds_from_simu_indices_type_n()
    
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)

