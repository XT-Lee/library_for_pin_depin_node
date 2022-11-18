import hoomd
import numpy
import hoomd.azplugins as azplugins

class multi_positions():
    R""" 
    v2.0beta
    Add an array of harmonic restraining potentials based on specified positions.[x][should _del_ be defined manually?]
    [x]In reset_multi_positions method, here exist num_traps_use and num_traps_record two types. Hence you'd better write __del__ code in c/cpp class.
    
    Args:
        traps_pos(float array_like): an array of n*3 which defines n potentials you will apply to particles.
            The 3rd column should be set as 0.0 to ensure that the traps are 2D traps.
        k (float or array_like): Force constant, isotropic or in each of *x*-, *y*-, and *z*-directions
        r_cut (non-negative float): Force constant, a scalar
        system(object): the object created and initialized by hoomd.init

    The Hamiltonian is augmented with a harmonic potential calculated based on the distance between the current and
    initial positions of particles. This effectively allows particles to have nearly fixed position while still retaining
    their degrees of freedom. The potential has the following form:
    .. math::
        :nowrap:
        \begin{equation*}
        V(\mathbf{r}) = \frac{1}{2} \mathbf{k} \mathbf{\Delta r} \mathbf{\Delta r}^T  where \Delta r < r_cut, 
        V(\mathbf{r}) = \frac{1}{2} \mathbf{k} \mathbf{r_cut} \mathbf{r_cut}^T  where \Delta r >= r_cut, 
        \end{equation*}
    The strength of the potential depends on a spring constant :math:`\mathbf{k}` so that the particle position can be
    restrained in any of the three coordinate directions. If :math:`\mathbf{k}` is very large, the particle position is
    essentially constrained. However, because the particles retain their degrees of freedom, shorter integration timesteps
    must be taken for large :math:`\mathbf{k}` to maintain stability.
    .. note::
        The displacement is calculated using the minimum image convention.
    Examples::
        import hoomd.azplugins.multi_positions as mp
        sys=hoomd.init.create_lattice(unitcell=hoomd.lattice.sq(a=5), n=9);
        traps_pos=[[10,10,0],[-10,-10,0]];
        mp.set_multi_positions(traps_pos,k=100.0,r_cut=6.0,system=sys)
    .. warning::
        Virial calculation is not implemented because the particles are tethered to fixed positions. A warning will be raised
        if any calls to :py:class:`hoomd.analyze.log` are made because the logger always requests the virial flags. However,
        this warning can be safely ignored if the pressure (tensor) is not being logged or the pressure is not of interest.
    Parameters:
        group: select all the particles
        trap_unit: a trap whose k=0.0,r_cut=0.0
        num_traps: the number of traps activated(k>0,r_cut>0) now
        num_traps_max: the maximum number of traps(no matter activated or not) having been created.
        traps_list: a list of trap-positions
        num_particles: the number of particles to apply traps on
        
    Methods:
        get_particle_number: get particle number in the given system.
        init_trap_unit: create a trap as a seed to generate an array of traps, whose k=0.0 and r_cut=0.0 
        set_multi_positions: create an array of traps.
        alter_multi_positions: if given new trap positions, this method will alter traps as inputting data.
    
    """
    def __init__(self,traps_pos,k,r_cut,system):
        self.get_particle_number(system)
        self.group=hoomd.group.all()
        self.init_trap_unit()
        self.set_multi_positions(traps_pos,k,r_cut)
        
    def init_trap_unit(self,kset=0.0,rcut=0.0):
        self.trap_unit = azplugins.restrain.position(group=self.group,k=kset,r_cut=rcut)
    
    def get_particle_number(self,system):
        snap = system.take_snapshot(all=True)
        self.num_particles=len(snap.particles.position[:,0])
        
    def set_multi_positions(self,traps_pos,k,r_cut):
        R"""
        Description:
            this method supports only multi-traps_pos and single k & r_cut.[x]
        Parameters:
            system: used to get the num of particles existing in the system. 
            traps_pos: giving the positions where traps should be set and the numbers of traps.
        """
        #get information of the traps to create.
        self.num_traps=self.check_the_num_of_inputted_traps(traps_pos)
        arr_traps = range(self.num_traps)
        #generate an array for following for loop
        arr_particles = range(self.num_particles)

        #generate an array of traps for object 'position'
        self.traps_list=[self.trap_unit]*self.num_traps
        
        #add_new_traps, alter_traps_positions?
        for (i) in  arr_traps :
            self.traps_list[i]=azplugins.restrain.position(group=self.group,k=k,r_cut=r_cut)
            for (j) in  arr_particles :
                self.traps_list[i].set_position(int(j),traps_pos[i][:])
        #update global data in the class
        self.num_traps_max=self.num_traps
    
    def alter_multi_positions(self,traps_pos_new,k,r_cut):
        #get information of the traps to alter as.
        self._num_traps_new=self.check_the_num_of_inputted_traps(traps_pos_new)
        self._num_traps_old=self.num_traps
        #reset the olds array and generate a new array of traps for object 'position'
        if self._num_traps_new<self._num_traps_old:
            #init redundant traps
            small=self._num_traps_new
            large=self._num_traps_old
            self._arr_to_zero=numpy.linspace(small,large-1,large-small)
            for i in  self._arr_to_zero.astype(int) :
                self.traps_list[i].set_params(k=0.0, r_cut=0.0)
            #alter traps positions
            self.alter_traps_positions(self._num_traps_new,traps_pos_new,k,r_cut)
            #update global data in the class
            self.num_traps=self._num_traps_new
            #delete used data
            del self._arr_to_zero,self._num_traps_new,self._num_traps_old

        else: #self._num_traps_new>self._num_traps_old:
            if self._num_traps_new > self.num_traps_max:
                #add new traps
                self.add_new_traps(self._num_traps_new)
                #update global data in the class
                self.num_traps_max=self._num_traps_new
                self.num_traps=self._num_traps_new
            #alter traps positions
            self.alter_traps_positions(self._num_traps_new,traps_pos_new,k,r_cut)
            

    def add_new_traps(self,num_traps_new):
        self._num_traps_old=self.num_traps
        self._num_traps_new=num_traps_new
        self._traps_list_new=[self.trap_unit]*self._num_traps_new   #init
        self._traps_list_new[0:self._num_traps_old]=self.traps_list    #succeed old traps
        #self.traps_list=self._traps_list_new
        #self._arr_to_add=self.alter_list(self._num_traps_old,self._num_traps_new)
        small=self._num_traps_old
        large=self._num_traps_new
        self._arr_to_add=numpy.linspace(small,large-1,large-small)
        for (i) in  self._arr_to_add.astype(int) :
            self._traps_list_new[i]=azplugins.restrain.position(group=self.group,k=0.0,r_cut=0.0)
        self.traps_list=self._traps_list_new
        #delete used data
        del self._num_traps_old,self._traps_list_new

    def alter_traps_positions(self,num_traps_new,traps_pos_new,k,r_cut):
        arr_particles = range(self.num_particles)
        arr_traps = range(num_traps_new)
        for i in  arr_traps :
            self.traps_list[i].set_params(k=k, r_cut=r_cut)
            for j in  arr_particles :
                self.traps_list[i].set_position(int(j),traps_pos_new[i][:]) 
    '''
    what the fuck bug? takes 2 arguments but 3 were given

    def alter_list(small,large):
        array_to_alter=numpy.linspace(small,large-1,large-small)
        return array_to_alter
    '''
    
    
    def check_the_num_of_inputted_traps(self,traps_pos):
        traps_shape = numpy.asarray(traps_pos)
        num_traps=traps_shape.shape[0]

        #k_shape = numpy.asarray(k)
        #num_k=k_shape[0]

        #rcut_shape = numpy.asarray(r_cut)
        #num_rcut=rcut_shape[0]
        return num_traps#,num_k,num_rcut
        '''
        if not (num_rcut==1):
            print('Error: you should input only one r_cut!')

        if num_k==1:
            #k and rcut should be set together
            pass
        elif not (num_k == 1):
            #k and rcut should be set separately
            pass
        else:
            pass
        '''
        
    
                