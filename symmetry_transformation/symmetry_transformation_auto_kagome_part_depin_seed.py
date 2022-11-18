import hoomd
import hoomd.md
import hoomd.azplugins.sequence_generator as sg
import numpy
import hoomd.azplugins.multi_positions as mp
from matplotlib import pyplot

def workflow(index1,k1,step,k_end,linear_compression_ratio,seed_set=9):
    R"""
    Introduction:
        The critical Linear Compression Ratio for kagome lattice is 0.866.

    """
    #set parameters
    rcut=1.0
    trap_prefix='/home/tplab/hoomd-examples_0/'
    trap_filename=trap_prefix+"testkagome_part3-11-6"#data of trap position
    #caution: LCR > 0.80 !

    num=(k_end-k1)/step+1
    num=round(num)#get the num of repeat times

    #simulation setup
    for i in numpy.linspace(1,num,num):
        #set parameters
        #linear_compression_ratio=linear_compression_ratio1+step*(i-1.0)
        kset=k1+step*(i-1.0)
        index=index1+(i-1.0)
        str_index=str(int(index))
        str_seed=str(int(seed_set))
        str_index=str_index+"_"+str_seed
        log_prefix='/home/tplab/hoomd-examples_0/'
        file_log=log_prefix+'log-output_auto'+str_index+'.log'
        file_gsd=log_prefix+'trajectory_auto'+str_index+'.gsd'
        prefix='/home/tplab/Downloads/'
        result_filename=prefix+'index'+str_index
        png_filename=prefix+'index'+str_index+'.png'

        #run a case with 10 different radom seeds.
        simulation1(linear_compression_ratio,trap_filename,kset,rcut,seed_set,file_log,file_gsd,result_filename,png_filename)
        print('this is the '+str(i.astype(int))+'th cycle\n')

    return index.astype(int) #as index_end


def simulation1(linear_compression_ratio,trap_filename,kset,rcut,seed_set,file_log,file_gsd,result_filename,png_filename=""):
    R"""
    introduction:
        the simulation shows the depinning process 
        from kagome through decreasing trap potential.

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
    seq.generate_kagome(a=a_set,n=[10,6])#initial state! Not traps!
    sys=seq.system#sys=hoomd.init.create_lattice(unitcell=hoomd.lattice.hex(a=3), n=[12,6])
        #ex_render.render_disk_frame(sys.take_snapshot())

    #set particles
    nl = hoomd.md.nlist.cell()
    Yukawa = hoomd.md.pair.yukawa(r_cut=15, nlist=nl)
    Yukawa.pair_coeff.set('A', 'A', epsilon=300.0,kappa=0.25)
    #set traps
    traps_pos=numpy.loadtxt(trap_filename)
    traps_pos=numpy.dot(linear_compression_ratio,traps_pos)

    mp.set_multi_positions(traps_pos,k=kset,r_cut=rcut,system=sys)

    #set system
    hoomd.md.integrate.mode_standard(dt=0.002)
    all = hoomd.group.all()
    hoomd.md.integrate.langevin(group=all, kT=1.0, seed=seed_set)
    
    period=2e3#to analyze data, the period should be the same!
    #given period, the last frame in .log and .gsd is not the real last frame!
    #if the fluctuation in system is large, period must be set small!
    hoomd.analyze.log(filename=file_log,
                      quantities=['potential_energy', 'temperature'],
                      period=period,
                      overwrite=True)
    
    hoomd.dump.gsd(file_gsd, period=period, group=all, overwrite=True)

    #run the system
    hoomd.run(2e4)

    #save simulated points as 'index123_seed.txt'
    #and it is the real last frame!
    snap=sys.take_snapshot()
    points=snap.particles.position[:]
    numpy.savetxt(result_filename,points)
    """
    #save simulated points as 'index123.png'
    box=snap.box
    pyplot.figure()
    #pyplot.cla()
    pyplot.scatter(points[:,0], points[:,1])
    pyplot.xlabel('x')
    pyplot.ylabel('y')
    pyplot.axis('equal')
    rate=0.5
    pyplot.xlim(box.Lx*(-rate),box.Lx*rate)
    pyplot.ylim(box.Ly*(-rate),box.Ly*rate)

    pyplot.savefig(png_filename)
    """
    