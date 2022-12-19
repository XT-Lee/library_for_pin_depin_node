#import threading
import time


tm1=time.localtime(time.time())


import pandas as pd
path_to_results = '/home/remote/xiaotian_file/data/20221129/video_7'
features_file = path_to_results+'/'+'feature.csv'
features = pd.read_csv(features_file)
print(features['y'].head())
features['y'] = 1420-1 - features['y']
print(features['y'].head())

"""
import data_analysis_cycle as da
import points_analysis_2D as pa
import numpy as np
path_to_results = '/home/remote/xiaotian_file/data/20221129/video_7'
txyz_npy_filename = path_to_results+'/'+'txyz_stable.npy'
txyz_stable = np.load(txyz_npy_filename)
msds = pa.dynamic_points_analysis_2d(txyz_stable,mode='exp')
msds.plot_trajectory_single_particle(path_to_results+'/traj_stable/')#trajectory_stable.png
#txyz_stable[:] = txyz_stable[:]/2.0
"""
"""
msds.compute_atmsd_t_chips(0.95)
time_log_file = path_to_results+'/'+'DefaultVideo_7.txt'
time_log = np.loadtxt(time_log_file)
png_filename=path_to_results+'/'+'msd_scantchips_loglog_um_95%.png'
msds.plot_msd(time_log,png_filename)

"""
"""



"""
#print(txyz_stable[0,:5])

#a0 = {('A'+unit):result.bond_length_median*lc}
#record_msd_id = pa.compute_atmsd_t_chips(0.9,True)

"""
ts_id_dxy['particle']
#pa.compute_atmsd_t_chips(interval_max=0.9,msd_per_particle=True)
path_to_folder = '/home/remote/xiaotian_file/data/20221129/video5/'
time_log = path_to_folder+'DefaultVideo_5.txt'
png_filename = path_to_folder+'particle_'+str(id)+'_'+'msds_loglog.png'
pa.plot_msd_uniform(m_msd=1,time_log=time_log,png_filename=png_filename)
"""

#time.sleep(1)
tm2=time.localtime(time.time())
#calculate the time cost
import computeTime as ct
ct.getTimeCost(tm1,tm2)

