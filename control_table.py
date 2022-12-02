#import threading
import time
import numpy as np
tm1=time.localtime(time.time())
import workflow_part as tt
index1=1183
lists = np.linspace(0.77,0.84,8)


for lcr in lists:
    tt.workflow_simu_to_mysql_depin_from_honeycomb_part(index1,lcr)
    #print(index1,lcr)
    index1=index1+10


import getDataAndScatter as scatt
scatt.workflow_mysql_to_data_depin_from_honeycomb_part1()

#time.sleep(1)
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)

