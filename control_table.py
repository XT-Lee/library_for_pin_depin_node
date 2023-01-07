#import threading
import time
tm1=time.localtime(time.time())


R"""
select * from pin_hex_to_honeycomb_part_klt_2m where HarmonicK = 700;
| SimuIndex | HarmonicK | LinearCompressionRatio | kT   | Psi3     | Psi6     | RandomSeed |
|      4302 |       700 |                   0.81 |    1 | 0.927068 | 0.123686 |          9 |

select * from pin_hex_to_honeycomb_klt_2m where HarmonicK = 900;
+-----------+-----------+------------------------+------+----------+----------+------------+
| SimuIndex | HarmonicK | LinearCompressionRatio | kT   | Psi3     | Psi6     | RandomSeed |
|      4634 |       900 |                   0.79 |    1 | 0.862018 | 0.159095 |          9 |

select * from pin_hex_to_honeycomb_klt_2m where SimuIndex = 5238;
+-----------+-----------+------------------------+------+----------+----------+------------+
| SimuIndex | HarmonicK | LinearCompressionRatio | kT   | Psi3     | Psi6     | RandomSeed |
+-----------+-----------+------------------------+------+----------+----------+------------+
|      5238 |        60 |                   0.79 |    1 | 0.885731 | 0.196146 |          9 |
+-----------+-----------+------------------------+------+----------+----------+------------+
"""
#import getDataAndScatter as scatt
#scatt.workflow_mysql_to_data_pin_hex_to_honeycomb_klt_2m(account='remote')

import data_analysis_cycle as da
get_traj = da.data_analysis()
directory,data_name = get_traj.gsd_to_txyz('remote',5238,9)
get_traj.txyz_to_bond_plot(directory,data_name)
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

