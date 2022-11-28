import hoomd
import hoomd.md
import hoomd.azplugins.sequence_generator as sg
import numpy
import hoomd.azplugins.multi_positions as mp
from matplotlib import pyplot

def workflow(index1,k1,step,k_end,linear_compression_ratio,seed_set):
    #set parameters
    #kset=500.0
    #rcut=1.0
    #linear_compression_ratio=0.816
    
    trap_prefix='/home/tplab/hoomd-examples_0/'
    trap_filename=trap_prefix+"testhoneycomb3-8-12"#data of trap position

    #index1=206
    #tune k
    #step=10.0
    #k1 = 0.0
    #k_end=90.0
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

        simulation1(linear_compression_ratio,seed_set,0.01,kset,5e4,file_log,file_gsd,result_filename,png_filename)
        print('this is the '+str(i.astype(int))+'th cycle\n')

    return index.astype(int) #as index_end


def simulation1(linear_compression_ratio,seed_set,dt_set,kT_set,step_set,file_log,file_gsd,result_filename,png_filename,plot=False):
    R"""
    parameters:
        index1: the index where to start from
        linear_compression_ratio: the ratio to which 
            the trapping postions will be compressed linearly.
        trap_filename: the file which gives trapping positions.
        i: this is the i-th cycle time
    """
    
    a_set=3*linear_compression_ratio
    seq=sg.sequence()
    seq.generate_honeycomb(a=a_set,n=[8,12])
    sys=seq.system#sys=hoomd.init.create_lattice(unitcell=hoomd.lattice.hex(a=3), n=[12,6]);
        #ex_render.render_disk_frame(sys.take_snapshot())

    #set particles
    nl = hoomd.md.nlist.cell();
    Yukawa = hoomd.md.pair.yukawa(r_cut=15, nlist=nl);
    Yukawa.pair_coeff.set('A', 'A', epsilon=300.0,kappa=0.25);


    #set system
    hoomd.md.integrate.mode_standard(dt=dt_set);
    all = hoomd.group.all();
    hoomd.md.integrate.langevin(group=all, kT=kT_set, seed=seed_set);
    
    
    hoomd.analyze.log(filename=file_log,
                      quantities=['potential_energy', 'temperature','pressure'],
                      period=100,
                      overwrite=True);
    
    hoomd.dump.gsd(file_gsd, period=2e3, group=all, overwrite=True);

    #run the system
    hoomd.run(step_set);

    #save simulated points as 'index123.txt'
    snap=sys.take_snapshot()
    box=snap.box
    points=snap.particles.position[:]
    
    numpy.savetxt(result_filename,points)

    #save simulated points as 'index123.png'
    if plot:
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