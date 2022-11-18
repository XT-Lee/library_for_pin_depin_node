import hoomd
import hoomd.md
import hoomd.azplugins.sequence_generator as sg
import numpy
import hoomd.azplugins.multi_positions as mp
from matplotlib import pyplot

class simulation_unit:
    R"""
    not finished

    """
    def __init__(self):
        self.workflow()
        
    def workflow(self):
        #set parameters
        self.kset=500.0
        self.rcut=1.0
        self.linear_compression_ratio=0.816
        
        self.trap_prefix='/home/tplab/hoomd-examples_0/'
        self.trap_filename=self.trap_prefix+"testhoneycomb3-8-12"#data of trap position

        self.index1=194
        #tune k
        self.step=100.0
        self.k1 = 100.0
        self.k_end=400.0
        self.num=(self.k_end-self.k1)/self.step+1
        self.num=round(self.num)#get the num of repeat times

        #simulation setup
        for i in numpy.linspace(1,self.num,self.num):
            #set parameters
            #linear_compression_ratio=linear_compression_ratio1+step*(i-1.0)
            self.kset=self.k1+self.step*(i-1.0)
            self.index=self.index1+(i-1.0)
            self.str_index=str(int(self.index))
            self.log_filename='log-output_auto'+self.str_index+'.log'
            self.file_gsd='trajectory_auto'+self.str_index+'.gsd'
            self.prefix='/home/tplab/Downloads/'
            self.result_filename=self.prefix+'index'+self.str_index
            self.png_filename=self.prefix+'index'+self.str_index+'.png'

            self.simulation()
            print('this is the '+str(i.astype(int))+'th cycle\n')


    def simulation(self):
        R"""
        parameters:
            index1: the index where to start from
            linear_compression_ratio: the ratio to which 
                the trapping postions will be compressed linearly.
            trap_filename: the file which gives trapping positions.
            i: this is the i-th cycle time
        """
        #initialization
        #hoomd.context.initialize('')
        
        a_set=3*0.816
        seq=sg.sequence()
        seq.generate_honeycomb(a=a_set,n=[8,12])
        self.sys=seq.system

        #set particles
        self.set_particles()
        
        #set traps
        self.set_traps()

        #set system and log
        self.set_system()
        self.set_log()
        #run the system
        hoomd.run(2e4)
        #save simulated points as 'index123.txt'
        self.save_points()
    
    def set_init(self):
        pass

    def set_particles(self):
        nl = hoomd.md.nlist.cell()
        Yukawa = hoomd.md.pair.yukawa(r_cut=15, nlist=nl)
        Yukawa.pair_coeff.set('A', 'A', epsilon=300.0,kappa=0.25)

    def set_traps(self):
        traps_pos=numpy.loadtxt(self.trap_filename)
        traps_pos=numpy.dot(self.linear_compression_ratio,traps_pos)
        mp.set_multi_positions(traps_pos,k=self.kset,r_cut=self.rcut,system=self.sys)

    def set_system_and_log(self,dt_set=0.002,seed_set=9):
        hoomd.md.integrate.mode_standard(dt=dt_set)
        all = hoomd.group.all()
        hoomd.md.integrate.langevin(group=all, kT=1.0, seed=seed_set)

    def set_log(self):
        hoomd.analyze.log(filename=self.log_filename,
                        quantities=['potential_energy', 'temperature'],
                        period=100,
                        overwrite=True)
        
        hoomd.dump.gsd(self.file_gsd, period=2e3, group=all, overwrite=True)

    def save_points(self):
        #save simulated points as 'index123.txt'
        snap=self.sys.take_snapshot()
        points=snap.particles.position[:]
        numpy.savetxt(self.result_filename,points)