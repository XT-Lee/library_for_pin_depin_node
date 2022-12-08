#import threading
import time

tm1=time.localtime(time.time())

#pin sequence-GPU
import workflow_part as tt
import numpy
lcr_list = numpy.linspace(0.77,0.85,9)
lcr_list[-1] = 0.816
print(lcr_list)
index1=5299
for lcr1 in lcr_list:
    end_index = tt.workflow_simu_to_mysql_pin_hex_to_honeycomb_part_oop_klt_2m(index1=index1,lcr=lcr1,account='remote')
    #print(index1,lcr1)
    index1=end_index+1 #index1=index1+10 #



#time.sleep(1)
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)

