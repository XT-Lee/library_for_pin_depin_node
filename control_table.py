#import threading
import time

tm1=time.localtime(time.time())

import getDataAndScatter as scatt
scatt.workflow_mysql_to_data_pin_hex_to_honeycomb_part_klt_2m(account='remote')
"""
import data_analysis_cycle as da
da.save_from_gsd(simu_index=4302,seed=9,
                    final_cut=True,
                    psik_plot=3,
                    account='remote')

"""



#time.sleep(1)
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)

