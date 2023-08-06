#import time

def getTimeCost(t1,t2):
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
    dt_d=t2.tm_mday-t1.tm_mday
    dt_h=t2.tm_hour-t1.tm_hour
    dt_m=t2.tm_min-t1.tm_min
    dt_s=t2.tm_sec-t1.tm_sec
    if dt_s<0:
        dt_s=dt_s+60
        dt_m=dt_m-1
    if dt_m<0:
        dt_m=dt_m+60
        dt_h=dt_h-1
    if dt_h<0:
        dt_h=dt_h+24
        dt_d=dt_d-1

    if dt_d>0.9:
        print("start",t1.tm_mday,":",t1.tm_hour,":",t1.tm_min,":",t1.tm_sec)
        print("end  ",t2.tm_mday,":",t2.tm_hour,":",t2.tm_min,":",t2.tm_sec)
        print("cost ",dt_d,":",dt_h,":",dt_m,":",dt_s)
    elif dt_d<-0.1:
        print("start",t1.tm_hour,":",t1.tm_min,":",t1.tm_sec)
        print("end  ",t2.tm_hour,":",t2.tm_min,":",t2.tm_sec)
        print("cost ",dt_h+30,":",dt_m,":",dt_s)
    elif dt_d<0.1:
        print("start",t1.tm_hour,":",t1.tm_min,":",t1.tm_sec)
        print("end  ",t2.tm_hour,":",t2.tm_min,":",t2.tm_sec)
        print("cost ",dt_h,":",dt_m,":",dt_s)
        