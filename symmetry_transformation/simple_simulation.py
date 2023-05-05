import hoomd
import hoomd.md
import numpy
import hoomd.azplugins.sequence_generator as sg

class workflow_uniform:
    R"""
        Introduction
        exp:
            import symmetry_transformation.simple_simulation as pin
            wk = pin.workflow_uniform(index1=5208,account='tplab',mode="--mode=cpu")
            end_index = wk.workflow()
    """
    def __init__(self,index1,account='tplab',linear_compression_ratio=1.0,kT=1.0,seed_set=9,mode="--mode=cpu"):
        R"""
            mode: '--mode=cpu' or  '--mode=gpu'
        """
        #set parameters
        self.account = account
        self.index_start = index1
        self.linear_compression_ratio = linear_compression_ratio
        self.kT = kT
        self.seed_set = seed_set
        self.mode=mode
        #this should be set independently!
        self.set_init_state_parameters(pin=True)

    def workflow(self):
        R"""
        Introduction:
            The critical Linear Compression Ratio for kagome lattice is 0.866.
        """
        #set parameters

        #simulation setup

        #set parameters
        index=self.index_start
        str_index=str(int(index))
        str_seed=str(int(self.seed_set))
        str_index=str_index+"_"+str_seed
        log_prefix='/home/'+self.account+'/hoomd-examples_0/'
        file_log=log_prefix+'log-output_auto'+str_index+'.log'#log-output_auto3232_8.log
        file_gsd=log_prefix+'trajectory_auto'+str_index+'.gsd'#log-output_auto3232_8.gsd
        prefix='/home/'+self.account+'/Downloads/'
        result_filename=prefix+'index'+str_index

        #run a case with 10 different radom seeds.
        self.__simulation1(self.seed_set,file_log,file_gsd,result_filename)

        return index #as index_end

    def __simulation1(self,seed_set,file_log,file_gsd,result_filename):
        R"""
        introduction:
            the simulation shows the symmetry transition process 
            from hexagonal to kagome through improving trap potential.
        
        parameters:
            self.index_start: the index where to start from
            linear_compression_ratio: the ratio to which 
                the trapping postions will be compressed linearly.
            trap_filename: the file which gives trapping positions.
            i: this is the i-th cycle time
        """
        #initialization
        hoomd.context.initialize(self.mode);#--mode=cpu
        self.__init_state_launch()
            #ex_render.render_disk_frame(self.sys.take_snapshot())

        #set particles
        nl = hoomd.md.nlist.cell();
        Yukawa = hoomd.md.pair.yukawa(r_cut=15, nlist=nl);
        Yukawa.pair_coeff.set('A', 'A', epsilon=300.0,kappa=0.25);

        #set system
        hoomd.md.integrate.mode_standard(dt=0.002);
        all = hoomd.group.all();
        hoomd.md.integrate.langevin(group=all, kT=self.kT, seed=seed_set);
        
        period=1000#to analyze data, the period should be the same!
        #given period, the last frame in .log and .gsd is not the real last frame!
        #if the fluctuation in system is large, period must be set small!
        hoomd.analyze.log(filename=file_log,
                        quantities=['potential_energy', 'temperature'],
                        period=period,
                        overwrite=True)
        
        hoomd.dump.gsd(file_gsd, period=period, group=all, overwrite=True)

        #run the system
        hoomd.run(1e6+1)#+1 to let the last period be record.
        #eg. 0-100 is for 101~200 steps while 0-100-200 is for 201 steps

        #save simulated points as 'index123_seed.txt'
        #and it is the real last frame!
        snap=self.sys.take_snapshot()
        points=snap.particles.position[:]
        numpy.savetxt(result_filename,points)

    def set_init_state_parameters(self,a=3,nx=16,ny=8,pin=False,depin_from_kagome=False,depin_from_honeycomb=False):
        R"""
        here pin = true just means that the init_state is hexagonal lattice.
        """
        self.a = a
        self.nx= nx
        self.ny = ny
        #choose one init state to launch(pin by default)
        self.pin = pin
        self.depin_from_kagome = depin_from_kagome
        self.depin_from_honeycomb = depin_from_honeycomb

    def __init_state_launch(self):
        if self.pin:
            self.__init_state_pin()
        if self.depin_from_kagome:
            self.__init_state_depin_from_kagome()
        if self.depin_from_honeycomb:    
            self.__init_state_depin_from_honeycomb()

    def __init_state_pin(self):
        self.sys=hoomd.init.create_lattice(unitcell=hoomd.lattice.hex(a=self.a), n=[self.nx,self.ny]);

    def __init_state_depin_from_kagome(self):
        a_set=3*self.linear_compression_ratio
        seq=sg.sequence()
        seq.generate_kagome(a=a_set,n=[self.nx,self.ny])#12,6 -> 11,6 ;20220830,20:04
        self.sys=seq.system
        #pin = hex + traps
        #depin = traps + traps
        #heat = gsd + traps; kT 
        #cool = gsd + traps
    def __init_state_depin_from_honeycomb(self):
        a_set=3*self.linear_compression_ratio
        seq=sg.sequence()
        seq.generate_honeycomb(a=a_set,n=[self.nx,self.ny])
        self.sys=seq.system
        "testhoneycomb3-8-12-part1"


#continue

class workflow_uniform_continue:
    R"""
        Introduction
        exp:
            import symmetry_transformation.simple_simulation as pin
            wk = pin.workflow_uniform(index1=5208,account='tplab',mode="--mode=cpu")
            end_index = wk.workflow()
    """
    def __init__(self,index1,index_old,account='tplab',kT=1.0,seed_set=9,mode="--mode=cpu"):
        R"""
            mode: '--mode=cpu' or  '--mode=gpu'
        """
        #set parameters
        self.account = account
        self.index_start = index1
        self.index_old = index_old
        self.kT = kT
        self.seed_set = seed_set
        self.mode=mode
        #this should be set independently!
        self.set_init_state_parameters(pin=True)

    def workflow(self):
        R"""
        Introduction:
            The critical Linear Compression Ratio for kagome lattice is 0.866.
        """
        #set parameters

        #simulation setup
        
        #set parameters
        index=self.index_start
        str_index=str(int(index))
        str_seed=str(int(self.seed_set))
        str_index=str_index+"_"+str_seed
        log_prefix='/home/'+self.account+'/hoomd-examples_0/'
        file_log=log_prefix+'log-output_auto'+str_index+'.log'#log-output_auto3232_8.log
        file_gsd=log_prefix+'trajectory_auto'+str_index+'.gsd'#log-output_auto3232_8.gsd
        prefix='/home/'+self.account+'/Downloads/'
        result_filename=prefix+'index'+str_index

        str_index_old = str(int(self.index_old))
        str_index_old=str_index_old+"_"+str_seed
        self.file_gsd_old = log_prefix+'trajectory_auto'+str_index_old+'.gsd'
        #run a case with 10 different radom seeds.
        self.__simulation1(self.seed_set,file_log,file_gsd,result_filename)

        return index #as index_end

    def __simulation1(self,seed_set,file_log,file_gsd,result_filename):
        R"""
        introduction:
            the simulation shows the symmetry transition process 
            from hexagonal to kagome through improving trap potential.
        
        parameters:
            self.index_start: the index where to start from
            linear_compression_ratio: the ratio to which 
                the trapping postions will be compressed linearly.
            trap_filename: the file which gives trapping positions.
            i: this is the i-th cycle time
        """
        #initialization
        hoomd.context.initialize(self.mode);#--mode=cpu
        self.sys=hoomd.init.read_gsd(self.file_gsd_old,frame=-1)
            #ex_render.render_disk_frame(self.sys.take_snapshot())

        #set particles
        nl = hoomd.md.nlist.cell();
        Yukawa = hoomd.md.pair.yukawa(r_cut=15, nlist=nl);
        Yukawa.pair_coeff.set('A', 'A', epsilon=300.0,kappa=0.25);

        #set system
        hoomd.md.integrate.mode_standard(dt=0.002);
        all = hoomd.group.all();
        hoomd.md.integrate.langevin(group=all, kT=self.kT, seed=seed_set);
        
        period=10#to analyze data, the period should be the same!
        #given period, the last frame in .log and .gsd is not the real last frame!
        #if the fluctuation in system is large, period must be set small!
        hoomd.analyze.log(filename=file_log,
                        quantities=['potential_energy', 'temperature'],
                        period=period,
                        overwrite=True)
        
        hoomd.dump.gsd(file_gsd, period=period, group=all, overwrite=True)

        #run the system
        hoomd.run(1e4+1)#+1 to let the last period be record.
        #eg. 0-100 is for 101~200 steps while 0-100-200 is for 201 steps

        #save simulated points as 'index123_seed.txt'
        #and it is the real last frame!
        snap=self.sys.take_snapshot()
        points=snap.particles.position[:]
        numpy.savetxt(result_filename,points)

    def set_init_state_parameters(self,a=3,nx=16,ny=8,pin=False,depin_from_kagome=False,depin_from_honeycomb=False):
        R"""
        here pin = true just means that the init_state is hexagonal lattice.
        """
        self.a = a
        self.nx= nx
        self.ny = ny
        #choose one init state to launch(pin by default)
        self.pin = pin
        self.depin_from_kagome = depin_from_kagome
        self.depin_from_honeycomb = depin_from_honeycomb

