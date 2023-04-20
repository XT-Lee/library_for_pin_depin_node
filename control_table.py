#import threading
import time
tm1=time.localtime(time.time())

#get reference positions
import numpy as np
import workflow_analysis as wa
spd = wa.show_disp_field()
#spd.get_points_plot()
pd = wa.show_polygon_dye()
pd.get_points_plot()#_xylim()






#time.sleep(1)
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)

