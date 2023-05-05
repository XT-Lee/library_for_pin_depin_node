#import threading
import time
import computeTime as ct
tm1=time.localtime(time.time())





import file_for_CUDA.test_cuda as tc
tcc = tc.comapre_speed_cpu_vs_gpu()
t1,t2,r2 = tcc.calculate_cpu(tcc.rm1)
t1,t2,r2c = tcc.calculate_gpu(tcc.rm1)


#time.sleep(1)
tm2=time.localtime(time.time())
#calculate the time cost
ct.getTimeCost(tm1,tm2)

