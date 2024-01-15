#import threading
import time
import computeTime as ct
tm1=time.localtime(time.time())

import os
prefix_write = "/home/lixt/home/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_1/"
gsds_filename = prefix_write+"trajectory_auto7607_0.gsd"
file_size_b = os.path.getsize(gsds_filename)
file_size_kb = file_size_b/1024
print(file_size_b)
print(file_size_kb)

tm2=time.localtime(time.time())
#calculate the time cost
ct.getTimeCost(tm1,tm2)

