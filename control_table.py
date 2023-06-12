#import threading
import time
import computeTime as ct
tm1=time.localtime(time.time())

"""import workflow_analysis as wa
prefix_image='image_to_proceed/'
wti = wa.show_waiting_time_interstitial_motion()
wti.trajectory_coarse_grain_single_particle()"""

"""
import data_analysis_cycle as dac
dac.saveIndexCN346PCairoSeed(2563,2572,0,10,0.681,9)
dac.saveIndexCN346PCairoSeed(2573,2582,0,1,0.681,9)
"""
import data_analysis_cycle as dac
daw = dac.data_analysis_workflow()
daw.get_info_from_mysql_bond()


"""import opertateOnMysql as osql
osql.loadDataToMysql()
osql.createTableInMysql('depin_from_cairo_egct2lcra',None,'pin_hex_to_cairo_egct2lcra')"""
"""
import workflow_part as tt
index1=2153#2583#
lcr1=0.60#less than 0.60 is dangerous! some particles may not effected by trap!
while lcr1<0.805:
    tt.workflow_simu_to_mysql_pin_hex_to_cairo_egct(index1=index1,lcr=lcr1,seed=9)
    index1=index1+10
    lcr1=lcr1+0.01
"""

#import getDataAndDiagram  as scatt
#scatt.workflow_mysql_to_data_pin_hex_to_cairo_egct()
"""spd.draw_bonds_conditional_ridge_oop(prefix_write,8,True)
spd.draw_bonds_conditional_ridge_oop(prefix_write,104,True)
spd.draw_bonds_conditional_ridge_oop(prefix_write,1427,True)"""


"""trm = pa.trajectory_module()
prefix = '/home/tplab/xiaotian_file/lxt_code_py/4302_9/'
fig,axs = plt.subplots(1,3,sharey=True)
ax0 = trm.trajectory_coarse_grain_general(None,2,[0,8],prefix,'trajectory_displacement',axs[0],fig)
ax1 = trm.trajectory_coarse_grain_general(None,2,[8,104],prefix,'trajectory_displacement',axs[1],fig)
ax2 = trm.trajectory_coarse_grain_general(None,2,[104,1427],prefix,'trajectory_displacement',axs[2],fig)
plt.show()
"""


tm2=time.localtime(time.time())
#calculate the time cost
ct.getTimeCost(tm1,tm2)

