import numpy as np
import cupy as cp#import numpy 
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from scipy.spatial import  Delaunay
from cupyx.scipy.spatial import distance#from scipy.spatial import  distance
import cudf as pd#import pandas as pd
#import scikit-cuda 
import matplotlib 
import os
import datetime as dt
import time



class comapre_speed_cpu_vs_gpu:
    R"""
    for RTX 2080,
    operate n row random matrixmultiply(matmul), O(n)~n^3,
    when row<<3k, cpu is faster;row>>3k, gpu faster 
    """
    def __init__(self):
        lengt = 3000#
        self.rm1 = np.random.rand(lengt,lengt)


    def calculate_cpu(self,rm1):
        tm1=dt.datetime.now()
        rm2 = np.matmul(rm1,rm1) 
        tm2=dt.datetime.now()
        self.get_time_cost(tm1,tm2)
        return tm1,tm2,rm2
    
    def calculate_gpu(self,rm1):
        tm1=dt.datetime.now()
        rm1c = cp.array(rm1)#cost 2s

        #tm1=dt.datetime.now()
        rm2c = cp.matmul(rm1c,rm1c)
        tm2=dt.datetime.now()
        self.get_time_cost(tm1,tm2)
        #rm2c = np.array(rm2c)
        return tm1,tm2,rm2c.get()
    
    def get_time_cost(self,t1,t2):
        R"""
        parameters"
            t1:start time as time.localtime()
            t2:end time as time.localtime()
            dt: time cost
        example:
            t1=time.localtime()
            time.sleep(63)
            t2=time.localtime()
        """
        #calculate dt
        dt_s=t2.second-t1.second
        dt_us=t2.microsecond - t1.microsecond

        if dt_us<0:
            dt_us=dt_us+1e6
            dt_s=dt_s-1
            

        if dt_s>0:
            print("cost ",dt_s,"(s):",dt_us/1000,"(ms)")
        else:
            print("cost ",dt_us/1000,"(ms)")