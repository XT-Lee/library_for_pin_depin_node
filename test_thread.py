#test_thread

import threading
import time

def run(n):
    print('task',n)
    time.sleep(1)
    print('2s')
    time.sleep(1)
    print('1s')
    time.sleep(1)
    print('0s')
    time.sleep(1)

'''
import threading
import test_thread as tt

t1 = threading.Thread(target=tt.run,args=('t1',))
t2 = threading.Thread(target=tt.run,args=('t2',))
t3 = threading.Thread(target=tt.run,args=('t3',))
t1.start()
t2.start()
t3.start()
'''

