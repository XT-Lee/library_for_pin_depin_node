from sys import prefix
import time
tm1=time.localtime(time.time())

import function_plot.example_plot as ep
fpm = ep.functions_plot_module()

#fpm.generate_compare()
t,y = fpm.generate_glass_list()#(b=1,m=0.99)
fpm.plot_function_glass(t,y)
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)