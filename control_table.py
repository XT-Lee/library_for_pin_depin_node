#import threading
import time

tm1=time.localtime(time.time())

import workflow_part as tt
import numpy
"""
x1=[True,False,False]
x2=[False,False]
x=numpy.max(x2)
print(x)
"""
#print(int(2999.7))
"""
import symmetry_transformation.simple_simulation as st 
hex = st.workflow_uniform(5208,'remote',kT=0.1,seed_set=5,mode='--mode=cpu')
hex.workflow()


"""
import points_analysis_2D as pa
gsd_data = pa.proceed_gsd_file(simu_index=5208,seed=9,account='remote')
gsd_data.get_trajectory_data()
#txyz=numpy.load('hex_txyz.npy') 
#box=numpy.load('hex_box.npy')
msd_class = pa.msd(gsd_data.txyz,gsd_data.box,account='remote')#pa.msd(txyz,box,account='remote')
msd_class.compute_t_chips()
msd_class.plot()


"""
import freud
msds = freud.msd.MSD(gsd_data.box)#the class is fault,,'direct'
msds.compute(positions=msd_class.txyz_stable)
import matplotlib.pyplot as plt 
plt.figure()
plt.plot(msds.msd)
plt.title("Mean Squared Displacement")
plt.xlabel("$t$")
plt.ylabel("MSD$(t)$")
png_filename = 'msd_'+'index5208_9'+'.png'
plt.savefig(png_filename)#png_filename
plt.close()
"""
#print(tr.shape)



#time.sleep(1)
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)

