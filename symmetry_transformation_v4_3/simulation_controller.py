import workflow_analysis as wa
import symmetry_transformation_v4_3.system_parameters_generators as pg
import symmetry_transformation_v4_3.simulation_core as sco
import numpy as np
import os

class simulation_controller_traps:
    def __init__(self):
        pass

    def generate_initial_state_hexagonal_particle_honeycomb_trap_single(self):
        R"""
        Introduction:
                a comprehensive workflow from creating initial state, 
            to record trajectory as a series of .png images, covering a 
            single simulation point (simu_index,seed,lcr,k)
        Example:
            import symmetry_transformation_v4_3.simulation_controller as sc
            import numpy as np
            sct = sc.simulation_controller_traps()
            sct.generate_initial_state_hexagonal_particle_honeycomb_trap()
        """
        

        prefix = "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
        filename_initial_state = 'particle_hex_and_trap_honeycomb'
        gsd_filename = prefix+filename_initial_state+'.gsd'

        
        a_particle = 3
        lcr = 0.81
        particles = wa.archimedean_tilings()
        particles.generate_type1(a=a_particle)
        n_size = [16,8]
        particle_points = particles.generate_lattices(n_size)

        traps = wa.archimedean_tilings()
        traps.generate_type3(a=a_particle*lcr)
        isg = pg.initial_state_generator()
        
        isg.set_new_gsd_file_2types(particles,n_size,particle_points,traps,gsd_filename)
        
        #check the initial configuration of points and traps
        isg = pg.initial_state_generator()
        isg.read_gsd_file( gsd_filename)
        points = isg.particles.position
        
        ids = np.array(isg.snap.particles.typeid)
        list_p = ids == 0
        list_t = ids == 1

        isg.snap.particles.types
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots()
        ax.scatter(points[list_p,0],points[list_p,1],color='k')#
        ax.scatter(points[list_t,0],points[list_t,1],color='r',marker = 'x')#
        #ax.scatter(dula[:,0],dula[:,1],facecolors='none',edgecolors='k')#,marker = 'x'
        ax.set_xlabel('x label')  # Add an x-label to the axes.
        ax.set_ylabel('y label')  # Add a y-label to the axes.
        ax.set_title("Simple Plot")  # Add a title to the axes
        ax.set_aspect('equal','box')
        plt.show()
        plt.close()
        print('particle_hex_and_trap_honeycomb')

        #operate simulation
        sim = sco.simulation_core_traps(9999,9)
        sim.seed=9
        sim.mode = 'cpu'
        sim.gauss_epsilon = -50
        sim.input_file_gsd = gsd_filename
        sim.operate_simulation()
       
        
        #save trap position txt
        isg = pg.initial_state_generator()
        isg.read_gsd_file( gsd_filename)
        points = isg.particles.position
        ids = np.array(isg.snap.particles.typeid)
        list_p = ids == 0
        list_t = ids == 1
        save_file_txt = prefix + filename_initial_state + '.txt'
        np.savetxt(save_file_txt,points[list_t]) 


        import data_analysis_cycle as dac
        dac.save_from_gsd(9999,9,bond_plot=True,show_traps=True,trap_filename=save_file_txt,trap_lcr=1.0,list_traps=list_t,account='remote')

    def generate_initial_state_hexagonal_particle_honeycomb_trap_scan_k(self):
        R"""
        Introduction:
                a comprehensive workflow from creating initial state, 
            to record trajectory as a series of .png images, covering a 
            group of scanning k simulation points (simu_index,seed,lcr,k)
        Example:
            import symmetry_transformation_v4_3.simulation_controller as sc
            import numpy as np
            sct = sc.simulation_controller_traps()
            sct.generate_initial_state_hexagonal_particle_honeycomb_trap()
        """
        #series setup
        list_simu_index = np.linspace(5773,5782,10,dtype=int)
        list_gauss_epsilon = np.linspace(-3,-30,10)

        prefix = "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
        filename_initial_state = 'particle_hex_and_trap_honeycomb'
        init_gsd_filename = prefix+filename_initial_state+'.gsd'#describe the lcr

        #set_initial_state_particles
        a_particle = 3
        lcr = 0.81
        particles = wa.archimedean_tilings()
        particles.generate_type1(a=a_particle)
        n_size = [16,8]
        particle_points = particles.generate_lattices(n_size)
        #set_initial_state_traps
        traps = wa.archimedean_tilings()
        traps.generate_type3(a=a_particle*lcr)
        isg = pg.initial_state_generator()
        

        isg.set_new_gsd_file_2types(particles,n_size,particle_points,traps,init_gsd_filename)
        
        #check the initial configuration of points and traps
        isg = pg.initial_state_generator()
        isg.read_gsd_file( init_gsd_filename)
        points = isg.particles.position
        
        ids = np.array(isg.snap.particles.typeid)
        list_p = ids == 0
        list_t = ids == 1

        import matplotlib.pyplot as plt
        fig,ax = plt.subplots()
        ax.scatter(points[list_p,0],points[list_p,1],color='k')#
        ax.scatter(points[list_t,0],points[list_t,1],color='r',marker = 'x')#
        #ax.scatter(dula[:,0],dula[:,1],facecolors='none',edgecolors='k')#,marker = 'x'
        ax.set_xlabel('x label')  # Add an x-label to the axes.
        ax.set_ylabel('y label')  # Add a y-label to the axes.
        ax.set_title("Simple Plot")  # Add a title to the axes
        ax.set_aspect('equal','box')
        plt.show()
        plt.close('all')
        print('particle_hex_and_trap_honeycomb')

        #for seed in range(10):
        seed = 9
        for i in range(len(list_simu_index)):
            #operate simulation
            simu_index = list_simu_index[i]

            sim = sco.simulation_core_traps(simu_index,seed)
            sim.seed=seed
            sim.mode = 'cpu'
            sim.gauss_epsilon = list_gauss_epsilon[i]#-50
            sim.input_file_gsd = init_gsd_filename
            sim.operate_simulation()
        
            
            """#save trap position txt
            isg = pg.initial_state_generator()
            isg.read_gsd_file( init_gsd_filename)
            points = isg.particles.position
            ids = np.array(isg.snap.particles.typeid)
            list_p = ids == 0
            list_t = ids == 1
            save_file_txt = prefix + filename_initial_state + '.txt'
            np.savetxt(save_file_txt,points[list_t]) 


            import data_analysis_cycle as dac
            dac.save_from_gsd(simu_index,seed,bond_plot=True,show_traps=True,trap_filename=save_file_txt,trap_lcr=1.0,list_traps=list_t,account='remote')
            """


    def generate_initial_state_hexagonal_particle_honeycomb_trap_scan_lcr(self):
        R"""
        Introduction:
                a comprehensive workflow from creating initial state, 
            to record trajectory as a series of .png images, covering a 
            group of scanning k simulation points (simu_index,seed,lcr,k)
        Example:
            import symmetry_transformation_v4_3.simulation_controller as sc
            import numpy as np
            sct = sc.simulation_controller_traps()
            sct.generate_initial_state_hexagonal_particle_honeycomb_trap()
        """
        #series setup
        list_simu_index = np.linspace(5773,5782,10,dtype=int)
        list_gauss_epsilon = np.linspace(-3,-30,10)
        a_particle = 3
        lcr = 0.81
        str_lcr = '081'#in case the string like '0.810000000032'
        n_size = [16,8]

        prefix = "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
        #filename_initial_state should be like 
        # 'init_particle_<lattice_type>_<a>_<nx>_<ny>_and_trap_<lattice_type>_<int_lcr>
        # to prevent the mixing of different initial state .gsd file
        filename_initial_state = 'init_particle_hex_'+str(a_particle)+'_'+str(n_size[0])+'_'+str(n_size[1])+'_and_trap_honeycomb_'+str_lcr
                
        init_gsd_filename = prefix+filename_initial_state+'.gsd'#describe the lcr

        #check if the folder exists
        isExists=os.path.exists(init_gsd_filename)
        if isExists:
            pass
        else:
            #set_initial_state_particles
            particles = wa.archimedean_tilings()
            particles.generate_type1(a=a_particle)
            
            particle_points = particles.generate_lattices(n_size)
            #set_initial_state_traps
            traps = wa.archimedean_tilings()
            traps.generate_type3(a=a_particle*lcr)
            isg = pg.initial_state_generator()
            
            isg.set_new_gsd_file_2types(particles,n_size,particle_points,traps,init_gsd_filename)
            
            #check the initial configuration of points and traps
            isg = pg.initial_state_generator()
            isg.read_gsd_file( init_gsd_filename)
            points = isg.particles.position
            
            ids = np.array(isg.snap.particles.typeid)
            list_p = ids == 0
            list_t = ids == 1
            #check the init state
            import matplotlib.pyplot as plt
            fig,ax = plt.subplots()
            ax.scatter(points[list_p,0],points[list_p,1],color='k')#
            ax.scatter(points[list_t,0],points[list_t,1],color='r',marker = 'x')#
            #ax.scatter(dula[:,0],dula[:,1],facecolors='none',edgecolors='k')#,marker = 'x'
            ax.set_xlabel('x label')  # Add an x-label to the axes.
            ax.set_ylabel('y label')  # Add a y-label to the axes.
            ax.set_title("Simple Plot")  # Add a title to the axes
            ax.set_aspect('equal','box')
            plt.show()
            plt.close('all')
            print('particle_hex_and_trap_honeycomb')

        #for seed in range(10):
        seed = 9
        for i in range(len(list_simu_index)):
            #operate simulation
            simu_index = list_simu_index[i]

            sim = sco.simulation_core_traps(simu_index,seed)
            sim.seed=seed
            sim.mode = 'cpu'
            sim.gauss_epsilon = list_gauss_epsilon[i]#-50
            sim.input_file_gsd = init_gsd_filename
            sim.operate_simulation()