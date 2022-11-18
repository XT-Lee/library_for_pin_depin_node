import hoomd
import hoomd.md
import hoomd.azplugins.sequence_generator as sg
from hoomd import azplugins
import numpy
import hoomd.azplugins.multi_positions as mp
from matplotlib import pyplot

def workflow(index1,k1,step,k_end,linear_compression_ratio,seed_set):
    R"""
    Parameters:
        initial state: testsquare3-20-20
        traps:testsquare3-20-20
        safe lcr: [0.85,0.99]
        critical lcr: 0.931
    """
    #set parameters
    rcut=1.0
    trap_prefix='/home/tplab/hoomd-examples_0/'
    trap_filename=trap_prefix+"testcairo3-6-6"#data of trap position

    """
    #generate a trap sequence
    pp = sg.sequence()
    pp.generate_honeycomb(a=3,n=[8,12])
    pp.save(trap_filename)
    """
    num=(k_end-k1)/step+1
    num=round(num)#get the num of repeat times

    #simulation setup
    for i in numpy.linspace(1,num,num):
        #set parameters
        #linear_compression_ratio=linear_compression_ratio1+step*(i-1.0)
        kset=k1+step*(i-1.0)
        index=index1+(i-1.0)
        str_index=str(int(index))
        log_prefix='/home/tplab/hoomd-examples_0/'
        file_log=log_prefix+'log-output_auto'+str_index+'.log'
        file_gsd=log_prefix+'trajectory_auto'+str_index+'.gsd'
        prefix='/home/tplab/Downloads/'
        result_filename=prefix+'index'+str_index
        png_filename=prefix+'index'+str_index+'.png'

        simulation1(linear_compression_ratio,trap_filename,kset,rcut,seed_set,file_log,file_gsd,result_filename,png_filename)
        #print('this is the '+str(i.astype(int))+'th cycle\n')

    return index.astype(int) #as index_end


def simulation1(linear_compression_ratio,trap_filename,kset,rcut,seed_set,file_log,file_gsd,result_filename,png_filename):
    R"""
    introduction:
        the simulation shows the symmetry transition process 
        from square to hexagonal through decreasing trap potential.
    
    parameters:
        index1: the index where to start from
        linear_compression_ratio: the ratio to which 
            the trapping postions will be compressed linearly.
        trap_filename: the file which gives trapping positions.
        i: this is the i-th cycle time
    """
    #initialization
    a_set=3*linear_compression_ratio
    seq=sg.sequence()
    seq.generate_sq(a=a_set,n=20)
    sys=seq.system

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
    hoomd.md.integrate.langevin(group=all, kT=0.5, seed=seed_set);
    
    
    hoomd.analyze.log(filename=file_log,
                      quantities=['potential_energy', 'temperature'],
                      period=100,
                      overwrite=True);
    
    hoomd.dump.gsd(file_gsd, period=2e3, group=all, overwrite=True);

    #run the system
    hoomd.run(2e4);

    #save simulated points as 'index123.txt'
    snap=sys.take_snapshot()
    box=snap.box
    points=snap.particles.position[:]
    
    numpy.savetxt(result_filename,points)

    #save simulated points as 'index123.png'
    savefig=False
    if savefig:
        pyplot.figure();
        #pyplot.cla()
        pyplot.scatter(points[:,0], points[:,1]);
        pyplot.xlabel('x');
        pyplot.ylabel('y');
        pyplot.axis('equal');
        rate=0.5
        pyplot.xlim(box.Lx*(-rate),box.Lx*rate)
        pyplot.ylim(box.Ly*(-rate),box.Ly*rate)
        pyplot.savefig(png_filename)