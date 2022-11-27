import hoomd
import hoomd.md
import hoomd.azplugins.sequence_generator as sg
import numpy
import hoomd.azplugins.multi_positions as mp
from matplotlib import pyplot
import gsd.hoomd
def workflow(index_old,index1,kset,linear_compression_ratio):
    #set default parameters
    rcut=1.0
    
    trap_prefix='/home/tplab/hoomd-examples_0/'
    trap_filename=trap_prefix+"testhoneycomb3-8-12-part1"#data of trap position


    #set parameters
    #linear_compression_ratio=linear_compression_ratio1+step*(i-1.0)
    
    str_index=str(int(index1))
    log_prefix='/home/tplab/hoomd-examples_0/'
    file_gsd_old=log_prefix+'trajectory_auto'+str(index_old)+'.gsd'#[x]cation: index new/old !!
    file_log=log_prefix+'log-output_auto'+str_index+'.log'
    file_gsd=log_prefix+'trajectory_auto'+str_index+'.gsd'
    prefix='/home/tplab/Downloads/'
    result_filename=prefix+'index'+str_index
    png_filename=prefix+'index'+str_index+'.png'

    simulation1_pin_step2(file_gsd_old,linear_compression_ratio,trap_filename,kset,rcut,file_log,file_gsd,result_filename,png_filename)


def simulation1_pin_step2(file_gsd_old,linear_compression_ratio,trap_filename,kset,rcut,file_log,file_gsd,result_filename,png_filename):
    R"""
    introduction:
        the simulation shows the pinning process from hex toward honeycomb 
        through increasing trap potential of one half of honeycomb traps. 
    parameters:
        index1: the index where to start from
        linear_compression_ratio: the ratio to which 
            the trapping postions will be compressed linearly.
        trap_filename: the file which gives trapping positions.
        i: this is the i-th cycle time
    """
    #initialization
    hoomd.context.initialize("--mode=gpu");
    sys=hoomd.init.read_gsd(file_gsd_old,frame=-1)
    
    #file_gsd="trajectory_auto1502.gsd"
    #data_gsd=gsd.hoomd.open(file_gsd)
    #snap=data_gsd.read_frame(1)
    #snap.constraints.group.all

    #set particles
    nl = hoomd.md.nlist.cell();
    Yukawa = hoomd.md.pair.yukawa(r_cut=15, nlist=nl);
    Yukawa.pair_coeff.set('A', 'A', epsilon=300.0,kappa=0.25);
    #set traps
    traps_pos=numpy.loadtxt(trap_filename);
    traps_pos=numpy.dot(linear_compression_ratio,traps_pos)

    mp.set_multi_positions(traps_pos,k=kset,r_cut=rcut,system=sys)

    #set system
    hoomd.md.integrate.mode_standard(dt=0.002);
    all = hoomd.group.all();
    hoomd.md.integrate.langevin(group=all, kT=1.0, seed=9);
    
    
    hoomd.analyze.log(filename=file_log,
                      quantities=['potential_energy', 'temperature'],
                      period=100,
                      overwrite=True);
    
    hoomd.dump.gsd(file_gsd, period=1e2, group=all, overwrite=True);#2e3 default

    #run the system
    hoomd.run(2e4);

    #save simulated points as 'index123.txt'
    snap=sys.take_snapshot()
    #box=snap.box
    points=snap.particles.position[:]
    
    numpy.savetxt(result_filename,points)
    #save simulated points as 'index123.png'
