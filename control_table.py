#import threading
import time
tm1=time.localtime(time.time())



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

