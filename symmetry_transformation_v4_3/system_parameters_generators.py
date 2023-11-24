import numpy as np
import proceed_file as pf#gsd.hoomd
import gsd.hoomd
#import workflow_analysis as wa

class initial_state_generator:
    def __init__(self):
        R"""
        import symmetry_transformation_v4_3.system_parameters_generators as spg
        isg = spg.initial_state_generator()
        isg.get_gsd_sample()
        """
        pass

    def get_gsd_sample(self,simu_index=4302,seed=9,account='remote'):
        prefix = "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
        str_index=str(int(simu_index))+'_'+str(seed)
        file_gsd = prefix+'trajectory_auto'+str_index+'.gsd'#+'_'+str(seed)
        self.gsd_data = pf.proceed_gsd_file(filename_gsd_seed=file_gsd,account=account)
        self.snap = self.gsd_data.trajectory._read_frame(0)
        
        self.angles = self.snap.angles
        self.bonds = self.snap.bonds
        self.box = self.snap.configuration.box#[LX,LY,LZ,0,0,0]
        self.constraints = self.snap.constraints
        self.dihedrals = self.snap.dihedrals
        self.pairs = self.snap.pairs
        self.particles = self.snap.particles
        print(self.snap.particles.typeid)#[0]*N_particles, the value of the i-th element means the type('A'or'B') of the i-th particle.
        #print(self.snap.particles.types)#['A'] mark the particles
        
        #print(self.snap.particles.orientation)#for spherical particle, list([1,0,0,0])
        #print(self.snap.particles.types)#types of particles,list('A','B')
    
    def read_gsd_file(self,file_gsd):
        #prefix = "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
        #file_gsd = prefix+'particle_and_trap.gsd'#'lattice.gsd'#'trajectory_auto'+str_index+'.gsd'#+'_'+str(seed)
        self.trajectory = gsd.hoomd.open(file_gsd)
        self.snap = self.trajectory.read_frame(0)# _read_frame
        self.particles = self.snap.particles
        print('ok')

    def set_new_gsd_file(self,lattice_pointer,n_size,positions):
        R"""
            lattice_pointer: 
                lp = workflow_analysis.archimedean_tilings(), after lp.generate_typex(),
                lp.a1,lp.a2,lp.a3 are lattice constant;
                lp.position record the positions of points in a single lattice
            n_size:
                [nx,ny] lattices the box records.

        """#https://gsd.readthedocs.io/en/v3.2.0/python-module-gsd.hoomd.html#gsd.hoomd.open
        snap = gsd.hoomd.Snapshot()#.Frame()
        nps = np.shape(positions)[0]
        snap.particles.N = nps
        snap.particles.position = positions
        snap.particles.orientation = [(1,0,0,0)]*nps
        snap.particles.typeid = [0]*nps#record the rank in types, 0 for 'A', 1 for 'B'.
        snap.particles.types = ['A']
        lp = lattice_pointer
        snap.configuration.box = [lp.a1[0]*n_size[0],lp.a2[1]*n_size[1],lp.a3[2],0,0,0]

        prefix = "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
        with gsd.hoomd.open(name=prefix+'lattice.gsd', mode='xb') as f:
            f.append(snap)
    
    def set_new_gsd_file_2types(self,particle_pointer,n_size,positions,trap_pointer,output_gsd_filename):
        R"""
        introduction:
            Generate a gsd file containing an array of particles and an array of traps
        parameters:
            lattice_pointer(substrate_pointer): 
                lp = workflow_analysis.archimedean_tilings(), after lp.generate_typex(),
                lp.a1,lp.a2,lp.a3 are lattice constant;
                lp.position record the positions of points in a single lattice
            n_size:
                [nx,ny] lattices the box records.
        example:
            particles = wa.archimedean_tilings()
            particles.generate_type8_part(a=3)
            n_size = [3,2]
            particle_points = particles.generate_lattices(n_size)

            traps = wa.archimedean_tilings()
            traps.generate_type10(a=0.8*3)# a*lcr !
            isg = pg.initial_state_generator()
            isg.set_new_gsd_file_2types(particles,n_size,particle_points,traps)

        example_show_result:
            import matplotlib.pyplot as plt
            import symmetry_transformation_v4_3.system_parameters_generators as pg
            isg = pg.initial_state_generator()
            isg.read_gsd_file()
            points = isg.particles.position
            import numpy as np
            ids = np.array(isg.snap.particles.typeid)
            list_p = ids == 0
            list_t = ids == 1

            isg.snap.particles.types
            fig,ax = plt.subplots()
            ax.scatter(points[list_p,0],points[list_p,1],color='k')#
            ax.scatter(points[list_t,0],points[list_t,1],color='r')#
            #ax.scatter(dula[:,0],dula[:,1],facecolors='none',edgecolors='k')#,marker = 'x'
            ax.set_xlabel('x label')  # Add an x-label to the axes.
            ax.set_ylabel('y label')  # Add a y-label to the axes.
            ax.set_title("Simple Plot")  # Add a title to the axes
            ax.set_aspect('equal','box')
            plt.show()
            ids = np.array(isg.snap.particles.typeid)

        """#https://gsd.readthedocs.io/en/v3.2.0/python-module-gsd.hoomd.html#gsd.hoomd.open
        snap = gsd.hoomd.Snapshot()#.Frame()

        #generate type 'trap' particles
        pp = particle_pointer
        snap.configuration.box = [pp.a1[0]*n_size[0],pp.a2[1]*n_size[1],pp.a3[2],0,0,0]
        sp = trap_pointer
        bx = snap.configuration.box
        trap_size = [int(bx[0]/sp.a1[0]),int(bx[1]/sp.a2[1])]
        trap_positions = sp.generate_lattices(trap_size)
        nps_trap = np.shape(trap_positions)[0]

        nps = np.shape(positions)[0]
        snap.particles.N = nps + nps_trap
        snap.particles.position = np.concatenate((positions,trap_positions),axis=0) 
        snap.particles.orientation = [(1,0,0,0)]*(nps+nps_trap)
        list_typeid = [0]*nps+[1]*nps_trap
        snap.particles.typeid = list_typeid
        snap.particles.types = ['particle','trap']
        

        #prefix = "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
        #prefix+'particle_and_trap.gsd'
        # wb: overwrite, xb: not overwrite, ab: append after a existing file
        with gsd.hoomd.open(name=output_gsd_filename, mode='wb') as f:
            f.append(snap)

    def generate_yukawa(self):
        #frame = gsd.hoomd.Frame()
        """frame.particles.N = N_particles
        frame.particles.position = position
        frame.particles.orientation = orientation
        frame.particles.typeid = [0] * N_particles
        frame.particles.types = ['octahedron']
        frame.configuration.box = [L, L, L, 0, 0, 0]
        with gsd.hoomd.open(name='lattice.gsd', mode='x') as f:
            f.append(frame)
        """
        pass
        

class potential_generator:
    def __init__(self):
        pass

    def generate_yukawa(self,a,kappa,r_min,r_max,width=1000):
        R"""
            u = a*exp(-kappa*r)/r
            f = -a*(kappa*r^-1 + r^-2)*exp(-kappa*r)
        """
        if r_min>0:
            r = np.linspace(r_min,r_max,width)
            list_index = -kappa*r
            u_unit = np.divide(np.exp(list_index),r)
            u = a*u_unit

            interval = (r_max-r_min)/float(width)
            f = np.zeros((width,))
            f[1:-1] = (u[2:] - u[:-2])/(2*interval)
            f[0] = (u[1] - u[0]) / interval
            f[-1 ]= (u[-1] - u[-2]) / interval

            list_u_and_f = np.zeros((width,2))
            list_u_and_f[:,0] = u
            list_u_and_f[:,1] = f
        else:
            print('Error: r_min should be no less than zero!')
            print('Error!')

        return u,f
    
    def save_txt(self,filename,table):
        np.savetxt(filename,table)