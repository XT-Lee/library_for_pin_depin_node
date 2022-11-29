#import threading
import time

tm1=time.localtime(time.time())

import workflow_part as tt
seed=9
#pin sequence-GPU
index1=4918
lcr1=0.80#less than 0.75 is dangerous! some particles may not effected by trap!
while lcr1<0.901:
    print(index1,lcr1,seed)
    end_index=tt.workflow_simu_to_mysql_pin_hex_to_kagome_oop_klt_2m(index1=index1,lcr=lcr1,kT=0.1,seed=seed,account='remote')
    index1=end_index+1#index1=index1+10#
    lcr1=lcr1+0.01
"""
import particle_tracking as pt
frame=pt.particle_track()
filename='/home/remote/xiaotian_file/20210417/DefaultImage.jpg'
Dia=19#25-95%
minmass=1000
frame.single_frame_particle_tracking(filename,Dia,minmass)
"""


  


"""
import data_analysis_cycle as da
i=4448
seed=9
account='remote'
da.save_from_gsd(simu_index=i,seed=seed,final_cut=True,
                                bond_plot =True,
                                show_traps=True,
                                trap_filename="/home/"+account+"/hoomd-examples_0/testkagome_part3-11-6",
                                trap_lcr=0.89,
                                account=account)
"""


"""
seed=9
index1=4176
lcr1=0.89
kT=0.1
#while lcr1 < 0.895:
index_end=tt.workflow_simu_to_mysql_pin_hex_to_kagome_cycle_oop_kT(index1=index1,lcr=lcr1,kT=kT,seed=seed)
print(index1)
print(kT)
index1 = index_end + 1
"""

#time.sleep(1)
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)

