#import threading
import time
import numpy 
tm1=time.localtime(time.time())

"""
x1=[True,False,False]
x2=[False,False]
x=numpy.max(x2)
print(x)
"""
import points_analysis_2D as pa
gsd_data = pa.proceed_gsd_file(simu_index=5208,seed=9)
gsd_data.get_trajectory_data()

msd_class = pa.msd(gsd_data.txyz,gsd_data.box)

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
#print(tr.shape)



"""
import gsd.hoomd
import freud
traj = gsd.hoomd.open('/home/tplab/hoomd-examples_0/trajectory_auto5208_9.gsd')
snap=traj.read_frame(-1)
msd = freud.msd.MSD( snap.configuration.box)#the class is fault,,'direct'
msd.compute(positions=traj[:,:,:3])#,images=pos_list
import matplotlib.pyplot as plt

plt.figure()
plt.plot(msd.msd)
plt.title("Mean Squared Displacement")
plt.xlabel("$t$")
plt.ylabel("MSD$(t)$")

plt.savefig('/home/tplab/Downloads/msd.png')

"""

#import data_analysis_cycle as da
#da.save_from_gsd(simu_index=5208,seed=9)



#time.sleep(1)
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)

