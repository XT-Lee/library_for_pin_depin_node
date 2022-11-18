#import threading
import time
import workflow_part as tt

tm1=time.localtime(time.time())

seed=0
index1=4709
lcr1=0.86
while seed<8.5:
    tt.workflow_simu_to_mysql_pin_hex_to_kagome_oop_klt_2m_kt25(index1=index1,lcr=lcr1,seed=seed,account='remote')
    #print(index1,lcr1,seed)
    seed=seed+1

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

#time.sleep(1)
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


tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)

