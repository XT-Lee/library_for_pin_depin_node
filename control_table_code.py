import time
tm1=time.localtime(time.time())


import symmetry_transformation_v4_3.list_code_simulation as lcs
lcs.pin_hex_to_hp_fill_brownian_workflow()
"""_series_of_gsd_file()
cn3s = sgf.get_cn3s_from_csv_honeycomb_part()
print(cn3s)"""
"""import matplotlib.pyplot as plt
import symmetry_transformation_v4_3.system_parameters_generators as pg
import numpy as np"""
"""import symmetry_transformation_v4_3.list_code_analysis as lca
agf = lca.analyze_a_gsd_file()
agf.analyze_gsd_files_and_record_as_csv()"""
"""import symmetry_transformation_v4_3.list_code_analysis as lca
agf = lca.analyze_a_gsd_file()
agf.get_bond_plot_from_a_gsd()"""
#dac.save_from_gsd(5793,0,-1,bond_plot=True,)
"""import workflow_analysis as wa
atp = wa.archimedean_tilings_polygon_dye()

atp.workflow_type8()"""

#sct.generate_simu_index_csv()



#sco.simulation_core_traps()

#pg.initial_state_generator()

"""particles = wa.archimedean_tilings()
particles.generate_type9(a=3)
n_size = [10,10]
particle_points = particles.generate_lattices(n_size)

traps = wa.archimedean_tilings()
traps.generate_type3(a=3)
isg = pg.initial_state_generator()
prefix = "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
output_gsd_filename = prefix+'particle_type9_and_trap_square.gsd'
isg.set_new_gsd_file_2types(particles,n_size,particle_points,traps,output_gsd_filename)
#isg.set_new_gsd_file(at,n_size,points)#get_gsd_sample()
isg = pg.initial_state_generator()
isg.read_gsd_file(output_gsd_filename)
points = isg.particles.position
import numpy as np
ids = np.array(isg.snap.particles.typeid)
list_p = ids == 0
list_t = ids == 1

isg.snap.particles.types
fig,ax = plt.subplots()
ax.scatter(points[list_p,0],points[list_p,1],color='k')#
#ax.scatter(points[list_t,0],points[list_t,1],color='r')#
#ax.scatter(dula[:,0],dula[:,1],facecolors='none',edgecolors='k')#,marker = 'x'
ax.set_xlabel('x label')  # Add an x-label to the axes.
ax.set_ylabel('y label')  # Add a y-label to the axes.
ax.set_title("Simple Plot")  # Add a title to the axes
ax.set_aspect('equal','box')
plt.show()"""

"""import gsd.hoomd as gh
gsd_prefix = '/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/'
gsd_filename = gsd_prefix+'trajectory_auto9902_9.gsd'
gsd_file = gh.open(gsd_filename)
print(len(gsd_file))
frame = gsd_file.read_frame(0)
pos = frame.particles.position
print('123')
#simulation
import symmetry_transformation_v4_3.simulation_core as ts
un = ts.simulation_core_traps()
un.operate_simulation()#simulation()


#read the file
import proceed_file as pf
gsd_prefix = '/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/'
gsd_filename = gsd_prefix+'trajectory_auto9902_9.gsd'
gf = pf.proceed_gsd_file(gsd_filename)
save_prefix= '/home/remote/Downloads/'
gf.get_trajectory_data(save_prefix,9902,9)
npy = save_prefix+'txyz_9902_9.npy'
txyz = np.load(npy)
snap = gf.trajectory.read_frame(0)
ids = snap.particles.typeid#configuration.box
list_p = ids == 0
list_t = ids == 1

sz = np.shape(txyz)
for i in range(sz[0]):
    fig,ax = plt.subplots()
    ax.scatter(txyz[i,list_p,0],txyz[i,list_p,1],color='k')#
    ax.scatter(txyz[i,list_t,0],txyz[i,list_t,1],color='r',marker = 'x')#
    #ax.scatter(dula[:,0],dula[:,1],facecolors='none',edgecolors='k')#,marker = 'x'
    ax.set_xlabel('x label')  # Add an x-label to the axes.
    ax.set_ylabel('y label')  # Add a y-label to the axes.
    ax.set_title("Simple Plot")  # Add a title to the axes
    ax.set_xlim([-8,8])
    ax.set_ylim([-10,10])
    ax.set_aspect('equal','box')
    
    fn = save_prefix + 'index9902_9_'+str(i)+'.png'
    fig.savefig(fn)
    plt.close()"""


"""import function_plot.example_plot as ep
fpm = ep.functions_plot_module()

#fpm.generate_compare()
t,y = fpm.generate_glass_list()#(b=1,m=0.99)
fpm.plot_function_glass(t,y)"""
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)