import hoomd
import hoomd.md
import numpy
import hoomd.azplugins.multi_positions as mp
import hoomd.azplugins.sequence_generator as sg

class workflow_uniform:
    R"""
        Introduction
    """
    def __init__(self,index1,account='tplab',k1=100.0,step=100.0,k_end=1000.0,linear_compression_ratio=1.0,kT=1.0,seed_set=9,
                trap_name="testkagome_cycle3-4-6",mode="--mode=cpu",period=1e3,steps=2e6+1):
        #set parameters
        self.account = account
        self.index_start = index1
        self.k1 =k1 #spring constant
        self.k_step = step
        self.k_end = k_end
        self.linear_compression_ratio = linear_compression_ratio
        self.kT = kT
        self.seed_set = seed_set
        self.trap_name = trap_name
        self.mode=mode
        self.period = period
        self.steps = steps
        #this should be set independently!
        self.set_init_state_parameters(pin_from_hex=True)

    def workflow(self):
        R"""
        Introduction:
            The critical Linear Compression Ratio for kagome lattice is 0.866.
        """
        #set parameters
        rcut=1.0
        trap_prefix='/home/'+self.account+'/hoomd-examples_0/'
        self.trap_filename=trap_prefix+self.trap_name#"testkagome_part3-11-6"#data of trap position
        #caution: LCR > 0.80 !

        num=(self.k_end-self.k1)/self.k_step+1
        num=round(num)#get the num of repeat times

        #simulation setup
        for i in numpy.linspace(1,num,num):
            #set parameters
            kset=self.k1+self.k_step*(i-1.0)
            index=self.index_start+(i-1.0)
            str_index=str(int(index))
            str_seed=str(int(self.seed_set))
            str_index=str_index+"_"+str_seed
            #log_prefix='/home/'+self.account+'/hoomd-examples_0/'
            log_prefix='/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/'
            file_log=log_prefix+'log-output_auto'+str_index+'.log'#log-output_auto3232_8.log
            file_gsd=log_prefix+'trajectory_auto'+str_index+'.gsd'#log-output_auto3232_8.gsd
            prefix='/home/'+self.account+'/Downloads/'
            result_filename=prefix+'index'+str_index
            png_filename=prefix+'index'+str_index+'.png'

            #run a case with 10 different radom seeds.
            self.__simulation1(self.linear_compression_ratio,self.trap_filename,kset,rcut,self.seed_set,file_log,file_gsd,result_filename,png_filename)
            print('this is the '+str(i.astype(int))+'th cycle\n')

        return index.astype(int) #as index_end

    def __simulation1(self,linear_compression_ratio,trap_filename,kset,rcut,seed_set,file_log,file_gsd,result_filename,png_filename):
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
        self.__launch_init_state()
            #ex_render.render_disk_frame(self.sys.take_snapshot())

        #set particles
        nl = hoomd.md.nlist.cell();
        Yukawa = hoomd.md.pair.yukawa(r_cut=15, nlist=nl);
        Yukawa.pair_coeff.set('A', 'A', epsilon=300.0,kappa=0.25);
        #set traps
        traps_pos=numpy.loadtxt(trap_filename);
        traps_pos=numpy.dot(linear_compression_ratio,traps_pos)

        mp.set_multi_positions(traps_pos,k=kset,r_cut=rcut,system=self.sys)

        #set system
        hoomd.md.integrate.mode_standard(dt=0.002);
        all = hoomd.group.all();
        hoomd.md.integrate.langevin(group=all, kT=self.kT, seed=seed_set);
        
        #period=1000#to analyze data, the period should be the same!
        #given period, the last frame in .log and .gsd is not the real last frame!
        #if the fluctuation in system is large, period must be set small!
        hoomd.analyze.log(filename=file_log,
                        quantities=['potential_energy', 'temperature'],
                        period=self.period,
                        overwrite=True)
        
        hoomd.dump.gsd(file_gsd, period=self.period, group=all, overwrite=True)

        #run the system
        hoomd.run(self.steps)#+1 to let the last period be record.
        #eg. 0-100 is for 101~200 steps while 0-100-200 is for 201 steps

        #save simulated points as 'index123_seed.txt'
        #and it is the real last frame!
        snap=self.sys.take_snapshot()
        points=snap.particles.position[:]
        numpy.savetxt(result_filename,points)

    def set_init_state_parameters(self,a=3,nx=16,ny=8,pin_from_hex=False,init_gsd=None,depin_from_kagome=False,depin_from_honeycomb=False):
        R"""
            init_gsd: .gsd filename
        """
        self.a = a
        self.nx= nx
        self.ny = ny
        #choose one init state to launch(pin by default)
        self.pin_from_hex = pin_from_hex
        self.depin_from_kagome = depin_from_kagome
        self.depin_from_honeycomb = depin_from_honeycomb
        if init_gsd is None:
            self.pin_from_gsd = False
        else:
            self.pin_from_gsd = True
            self.init_gsd = init_gsd

    def __launch_init_state(self):
        if self.pin_from_hex:
            self.__init_state_pin_from_hex()
        if self.pin_from_gsd:
            self.__init_state_pin_from_gsd()
        if self.depin_from_kagome:
            self.__init_state_depin_from_kagome()
        if self.depin_from_honeycomb:    
            self.__init_state_depin_from_honeycomb()

    def __init_state_pin_from_hex(self):
        self.sys=hoomd.init.create_lattice(unitcell=hoomd.lattice.hex(a=self.a), n=[self.nx,self.ny]);

    def __init_state_pin_from_gsd(self):
        self.sys = hoomd.init.read_gsd(filename=self.init_gsd,frame=-1)

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
        seq.generate_honeycomb(a=a_set,n=[8,12])
        self.sys=seq.system
        "testhoneycomb3-8-12-part1"

class workflow_pin:
    R"""
        Introduction
    """
    def __init__(self,index1,account='tplab',k1=100.0,step=100.0,k_end=1000.0,linear_compression_ratio=1.0,seed_set=9,trap_name="testkagome_part3-11-6"):
        #set parameters
        self.account = account
        self.index_start = index1
        self.k1 =k1 #spring constant
        self.k_step = step
        self.k_end = k_end
        self.linear_compression_ratio = linear_compression_ratio
        self.seed_set = seed_set
        self.trap_name = trap_name

        self.set_init_state_pin()

    def workflow(self):
        R"""
        Introduction:
            The critical Linear Compression Ratio for kagome lattice is 0.866.
        """
        #set parameters
        rcut=1.0
        trap_prefix='/home/remote/hoomd-examples_0/'
        self.trap_filename=trap_prefix+self.trap_name#"testkagome_part3-11-6"#data of trap position
        #caution: LCR > 0.80 !

        num=(self.k_end-self.k1)/self.k_step+1
        num=round(num)#get the num of repeat times

        #simulation setup
        for i in numpy.linspace(1,num,num):
            #set parameters
            kset=self.k1+self.k_step*(i-1.0)
            index=self.index_start+(i-1.0)
            str_index=str(int(index))
            str_seed=str(int(self.seed_set))
            str_index=str_index+"_"+str_seed
            log_prefix='/home/remote/hoomd-examples_0/'
            file_log=log_prefix+'log-output_auto'+str_index+'.log'#log-output_auto3232_8.log
            file_gsd=log_prefix+'trajectory_auto'+str_index+'.gsd'#log-output_auto3232_8.gsd
            prefix='/home/remote/Downloads/'
            result_filename=prefix+'index'+str_index
            png_filename=prefix+'index'+str_index+'.png'

            #run a case with 10 different radom seeds.
            self.simulation1(self.linear_compression_ratio,self.trap_filename,kset,rcut,self.seed_set,file_log,file_gsd,result_filename,png_filename)
            print('this is the '+str(i.astype(int))+'th cycle\n')

        return index.astype(int) #as index_end

    def simulation1(self,linear_compression_ratio,trap_filename,kset,rcut,seed_set,file_log,file_gsd,result_filename,png_filename):
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
        hoomd.context.initialize("");#--mode=cpu
        self.__init_state_pin()
            #ex_render.render_disk_frame(self.sys.take_snapshot())

        #set particles
        nl = hoomd.md.nlist.cell();
        Yukawa = hoomd.md.pair.yukawa(r_cut=15, nlist=nl);
        Yukawa.pair_coeff.set('A', 'A', epsilon=300.0,kappa=0.25);
        #set traps
        traps_pos=numpy.loadtxt(trap_filename);
        traps_pos=numpy.dot(linear_compression_ratio,traps_pos)

        mp.set_multi_positions(traps_pos,k=kset,r_cut=rcut,system=self.sys)

        #set system
        hoomd.md.integrate.mode_standard(dt=0.002);
        all = hoomd.group.all();
        hoomd.md.integrate.langevin(group=all, kT=1.0, seed=seed_set);
        
        period=100#to analyze data, the period should be the same!
        #given period, the last frame in .log and .gsd is not the real last frame!
        #if the fluctuation in system is large, period must be set small!
        hoomd.analyze.log(filename=file_log,
                        quantities=['potential_energy', 'temperature'],
                        period=period,
                        overwrite=True)
        
        hoomd.dump.gsd(file_gsd, period=period, group=all, overwrite=True)

        #run the system
        hoomd.run(2e4+1)#+1 to let the last period be record.
        #eg. 0-100 is for 101~200 steps while 0-100-200 is for 201 steps

        #save simulated points as 'index123_seed.txt'
        #and it is the real last frame!
        snap=self.sys.take_snapshot()
        points=snap.particles.position[:]
        numpy.savetxt(result_filename,points)

    def set_init_state_pin(self,a=3,nx=16,ny=8):
        self.a = a
        self.nx= nx
        self.ny = ny

    def __init_state_pin(self):
        self.sys=hoomd.init.create_lattice(unitcell=hoomd.lattice.hex(a=self.a), n=[self.nx,self.ny]);
        #pin = hex + traps
        #depin = traps + traps
        #heat = gsd + traps; kT 
        #cool = gsd + traps


class workflow_depin:#[x]
    R"""
        Introduction
    """
    def __init__(self,index1,account='tplab',k1=100.0,step=100.0,k_end=1000.0,linear_compression_ratio=1.0,kT=1.0,seed_set=9,trap_name="testkagome_cycle3-4-6"):
        #set parameters
        self.account = account
        self.index_start = index1
        self.k1 =k1 #spring constant
        self.k_step = step
        self.k_end = k_end
        self.linear_compression_ratio = linear_compression_ratio
        self.kT = kT
        self.seed_set = seed_set
        self.trap_name = trap_name

        self.set_init_state_depin()

    def workflow(self):
        R"""
        Introduction:
            The critical Linear Compression Ratio for kagome lattice is 0.866.
        """
        #set parameters
        rcut=1.0
        trap_prefix='/home/'+self.account+'/hoomd-examples_0/'
        self.trap_filename=trap_prefix+self.trap_name#"testkagome_part3-11-6"#data of trap position
        #caution: LCR > 0.80 !

        num=(self.k_end-self.k1)/self.k_step+1
        num=round(num)#get the num of repeat times

        #simulation setup
        for i in numpy.linspace(1,num,num):
            #set parameters
            kset=self.k1+self.k_step*(i-1.0)
            index=self.index_start+(i-1.0)
            str_index=str(int(index))
            str_seed=str(int(self.seed_set))
            str_index=str_index+"_"+str_seed
            log_prefix='/home/'+self.account+'/hoomd-examples_0/'
            file_log=log_prefix+'log-output_auto'+str_index+'.log'#log-output_auto3232_8.log
            file_gsd=log_prefix+'trajectory_auto'+str_index+'.gsd'#log-output_auto3232_8.gsd
            prefix='/home/'+self.account+'/Downloads/'
            result_filename=prefix+'index'+str_index
            png_filename=prefix+'index'+str_index+'.png'

            #run a case with 10 different radom seeds.
            self.simulation1(self.linear_compression_ratio,self.trap_filename,kset,rcut,self.seed_set,file_log,file_gsd,result_filename,png_filename)
            print('this is the '+str(i.astype(int))+'th cycle\n')

        return index.astype(int) #as index_end

    def simulation1(self,linear_compression_ratio,trap_filename,kset,rcut,seed_set,file_log,file_gsd,result_filename,png_filename):
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
        hoomd.context.initialize("");#--mode=cpu
        self.__init_state_depin_from_kagome()
            #ex_render.render_disk_frame(self.sys.take_snapshot())

        #set particles
        nl = hoomd.md.nlist.cell();
        Yukawa = hoomd.md.pair.yukawa(r_cut=15, nlist=nl);
        Yukawa.pair_coeff.set('A', 'A', epsilon=300.0,kappa=0.25);
        #set traps
        traps_pos=numpy.loadtxt(trap_filename);
        traps_pos=numpy.dot(linear_compression_ratio,traps_pos)

        mp.set_multi_positions(traps_pos,k=kset,r_cut=rcut,system=self.sys)

        #set system
        hoomd.md.integrate.mode_standard(dt=0.002);
        all = hoomd.group.all();
        hoomd.md.integrate.langevin(group=all, kT=self.kT, seed=seed_set);
        
        period=100#to analyze data, the period should be the same!
        #given period, the last frame in .log and .gsd is not the real last frame!
        #if the fluctuation in system is large, period must be set small!
        hoomd.analyze.log(filename=file_log,
                        quantities=['potential_energy', 'temperature'],
                        period=period,
                        overwrite=True)
        
        hoomd.dump.gsd(file_gsd, period=period, group=all, overwrite=True)

        #run the system
        hoomd.run(2e4+1)#+1 to let the last period be record.
        #eg. 0-100 is for 101~200 steps while 0-100-200 is for 201 steps

        #save simulated points as 'index123_seed.txt'
        #and it is the real last frame!
        snap=self.sys.take_snapshot()
        points=snap.particles.position[:]
        numpy.savetxt(result_filename,points)

    def set_init_state_depin(self,a=3,nx=16,ny=8):
        self.a = a
        self.nx= nx
        self.ny = ny

    def __init_state_depin_from_kagome(self):
        a_set=3*self.linear_compression_ratio
        seq=sg.sequence()
        seq.generate_kagome(a=a_set,n=[12,6])
        self.sys=seq.system
        #pin = hex + traps
        #depin = traps + traps
        #heat = gsd + traps; kT 
        #cool = gsd + traps

