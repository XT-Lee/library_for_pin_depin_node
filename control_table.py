#import threading
import time
import computeTime as ct
tm1=time.localtime(time.time())

#draw a diagram of CN3 4 6
#draw a series of plot about CN3 4 6 vs t
import workflow_analysis as wa
import matplotlib.pyplot as plt
at = wa.archimedean_tilings()
at.generate_type11()
points = at.generate_lattices([18,12])
dula = at.get_dual_lattice(points)
fig,ax = plt.subplots()
ax.scatter(points[:,0],points[:,1])
ax.scatter(dula[:,0],dula[:,1],marker = 'x')#
ax.set_xlabel('x label')  # Add an x-label to the axes.
ax.set_ylabel('y label')  # Add a y-label to the axes.
ax.set_title("Simple Plot")  # Add a title to the axes
ax.set_aspect('equal','box')
plt.show()
"""getDataAndDiagram.get_dual_lattice()
dl = wa.show_dual_lattice()
dl.go()"""
"""import points_analysis_2D as pa
import numpy as np
list_index = np.linspace(4436,4445,10,dtype=int)
prefix = '/home/remote/Downloads/index'
postfix = '_9'
for i in range(10):
    filename = prefix+str(list_index[i])+postfix
    data = np.loadtxt(filename)
    spa = pa.static_points_analysis_2d(points=data[:,:2])
    spa.bond_first_minima_left()
    bond_cut_off=6
    trap_lcr=lcr
    trap_filename='/home/remote/hoomd-examples_0/testhoneycomb3-8-12-part1'
    self.x_unit = '($\sigma$)'
    png_filename1 = folder_name +'bond_hist_index'+str_index+'_'+str(int(i))+'.png'#+"/"
    png_filename2 = folder_name +'bond_plot_1st_minima_index'+str_index+'_'+str(int(i))+'.png'#+"/"
    
    a_frame.get_first_minima_bond_length_distribution(lattice_constant=1,hist_cutoff=bond_cut_off,png_filename=png_filename1)#,png_filename=png_filename1
    #a_frame.draw_bonds_conditional_bond(check=[0.4, a_frame.bond_first_minima_left], png_filename=png_filename2,
    #                               show_traps=show_traps,LinearCompressionRatio=trap_lcr,trap_filename=trap_filename,
    #                                nb_change=ids,x_unit=self.x_unit)
    a_frame.draw_bonds_conditional_bond_oop(check=[0.4, a_frame.bond_first_minima_left], png_filename=png_filename2,
                                            xy_stable=xy,x_unit=self.x_unit,
                                    LinearCompressionRatio=trap_lcr, trap_filename=trap_filename)"""

print('s')
"""
import tensorflow as tf
import jax
import cmake_example as ce
res = ce.add(1,2)
print(res)"""

"""import workflow_analysis as wa
prefix_image='image_to_proceed/'
wti = wa.show_waiting_time_interstitial_motion()
wti.trajectory_coarse_grain_single_particle()"""

#pcairo evolution

"""
import data_analysis_cycle as dac
dac.saveIndexCN346PCairoSeed(2563,2572,0,10,0.681,9)
dac.saveIndexCN346PCairoSeed(2573,2582,0,1,0.681,9)
"""
"""import opertateOnMysql as osql
res = osql.showTables('like \'pin_hex_to_kago%\'')"""


"""
import data_analysis_cycle as dac
simu_indexs = [5244,5245,5247]#[5236,5237,5238]
for simu_index in simu_indexs:
    dac.save_from_gsd_to_cn3(simu_index,9,coordination_number=True,account='remote')"""
"""import opertateOnMysql as osql
#res = osql.showTables('like \'pin%\'')
con = "where HarmonicK<70"
res = osql.getDataFromMysql(table_name='pin_hex_to_honeycomb_klt_2m',search_condition=con)
import workflow_part as tt
import numpy
list_lcr = numpy.linspace(0.77,0.85,9)
list_lcr[-1] = 0.816
list_seed = [6,7,8]
#print(lcr_list)
for seed in list_seed:
    index1=5209#~5298
    for lcr1 in list_lcr:
        end_index = tt.workflow_simu_to_mysql_pin_hex_to_honeycomb_oop_klt_2m(
            index1=index1,lcr=lcr1,seed=seed,account='remote',check=True)
        print(index1,end_index,lcr1,seed)
        index1=end_index+1 """
"""import data_analysis_cycle as dac
daw = dac.data_analysis_workflow()
daw.get_info_from_mysql_bond()
"""

import workflow_analysis as wa
#cn3 = wa.controller_get_honey_part_cn3_vs_u_sub()
#cn3.control_tb()
ccn3 = wa.controller_get_honey_part_cn3_vs_u_sub()
ccn3.control_tb()
"""import matplotlib.pyplot as plt
plt.figure()
plt.scatter(extended_positions[:,0],extended_positions[:,1])
plt.axis('equal')
plt.show()



import workflow_analysis as wa
spio = wa.show_pin_interstitial_order_parameter()
#spio.workflow_data_honey_part(io_only=True)#0.79,5238,9,io_only=False
spio.workflow_data_honey(io_only=True)#0.79,5238,9,io_only=False
import pandas as pd
simu_index = '5238_9'
csv_filename = '/home/remote/Downloads/'+simu_index+'/pin_check/time_vs_activation_ratio.csv'
data = pd.read_csv(csv_filename)#/home/remote/Downloads/5238_9/pin_check/time_vs_activation_ratio.csv
print(data.columns)
t = data['t(step)'].values
y1 = data['pin_ratio(1)'].values
y2 = data['pin_act_ratio(1)'].values
y3 = data['depin_act_ratio(1)'].values
spio.plot_to_check_act_ratio
spio.plot_to_check_act_share_ratio(t,y1,y2,'5238_9')#'t(step)', 'pin_ratio(1)', 'pin_act_ratio(1)'"""
"""
import particle_tracking as pt
pa = pt.particle_track()
prefix = '/home/remote/Downloads/image_proceed/honey_part/'
image_filename = 'DefaultImage_10.jpg'
x0 = 141
y0 = 141
limit = [x0,1023-x0,y0,1023-y0]
diameter = 17
f1 = pa.single_frame_particle_tracking(prefix+image_filename,D=diameter,minmass=1400,calibration=False,axis_limit=limit)#D = 16.7 pixel
import particle_decorate as pd 
pdr = pd.particle_decorator((1023-2*x0+2,1023-2*x0+2),pa.xy,diameter) #((100,100),[[50,50],[0,0]],100)#  
save_filename = '/home/remote/Downloads/image_proceed/honey_part/img10.png'
pdr.draw_raw_image(f1,save_filename)

"""
#import trackpy.preprocessing as pp
#pp.bandpass()

"""import numpy as np
list_index = np.linspace(2675,2679,5)
import data_analysis_cycle as dac
for index1 in list_index:
    #cairo simulation 2675-2679 are saved in account = 'tplab'.
    dac.save_from_gsd(simu_index=index1,seed=9,coordination_number=True,p_cairo=True)#,account='remote'"""

"""
import opertateOnMysql as osql
osql.loadDataToMysql()
osql.createTableInMysql('depin_from_cairo_egct2lcra',None,'pin_hex_to_cairo_egct2lcra')
import workflow_part as tt
index1=2153#2583#
lcr1=0.60#less than 0.60 is dangerous! some particles may not effected by trap!
while lcr1<0.805:
    tt.workflow_simu_to_mysql_pin_hex_to_cairo_egct(index1=index1,lcr=lcr1,seed=9)
    index1=index1+10
    lcr1=lcr1+0.01

"""



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

