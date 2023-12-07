import workflow_analysis as wa
import symmetry_transformation_v4_3.system_parameters_generators as pg
import symmetry_transformation_v4_3.simulation_core as sco
import numpy as np
import os
import pandas as pd

class simulation_controller_honeycomb_traps:
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
        sim.operate_simulation_langevin()
       
        
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
        n_files = 10
        list_simu_index = np.linspace(5773,5782,n_files,dtype=int)
        list_gauss_epsilon = np.linspace(-3,-30,n_files)

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
            sim.operate_simulation_langevin()
        
            
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


    def generate_initial_state_hexagonal_particle_honeycomb_trap_scan_csv(self,csv_filename):
        R"""
        Introduction:
                a comprehensive workflow from creating initial state, 
            to record trajectory as a series of .png images, covering a 
            group of scanning k simulation points (simu_index,seed,lcr,k)

            csv format: column = ['simu_index','seed','lcr','trap_gauss_epsilon','temperature']
        Example:
            import symmetry_transformation_v4_3.simulation_controller as sc
            import numpy as np
            sct = sc.simulation_controller_traps()
            sct.generate_initial_state_hexagonal_particle_honeycomb_trap()
        """
        #series setup
        df = pd.read_csv(csv_filename)
        
        col = df.columns
        print(df.head())
        
        n_size = [16,8]
        a_particle = 3
        list_simu_index = df[col[1]].values
        list_seed = df[col[2]].values
        list_lcr = df[col[3]].values
        list_trap_gauss_epsilon = df[col[4]].values
        list_temperature = df[col[5]].values

        for i in range(len(list_simu_index)):
            #operate simulation
            simu_index = list_simu_index[i]
            seed = list_seed[i]
            lcr = list_lcr[i]
            gauss_epsilon = list_trap_gauss_epsilon[i]
            temperature = list_temperature[i]
            init_gsd_filename = self.generate_initial_gsd(a_particle,n_size,lcr)
            """print(simu_index,seed,lcr,gauss_epsilon,temperature,init_gsd_filename)
            sim = sco.simulation_core_traps(simu_index,seed)
            sim.seed=seed
            sim.mode = 'gpu'
            sim.gauss_epsilon = gauss_epsilon#-50
            sim.kT = temperature
            sim.input_file_gsd = init_gsd_filename
            sim.operate_simulation()"""
            print(i+1,'/',len(list_simu_index))

    def generate_initial_gsd(sel,a_particle,n_size,lcr):
        prefix = "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
        #filename_initial_state should be like 
        # 'init_particle_<lattice_type>_<a>_<nx>_<ny>_and_trap_<lattice_type>_<int_lcr>
        # to prevent the mixing of different initial state .gsd file
        str_lcr = str(lcr)#'081'#in case the string like '0.810000000032'
        filename_initial_state = 'init_particle_hex_'+str(a_particle)+'_'+str(n_size[0])+'_'+str(n_size[1])+'_and_trap_honeycomb_'+str_lcr
                
        init_gsd_filename = prefix+filename_initial_state+'.gsd'#describe the lcr
         
        #check if the file exists
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
            
            """#check the initial configuration of points and traps
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
            print('particle_hex_and_trap_honeycomb')"""
            print('create:'+init_gsd_filename)
        return init_gsd_filename

    def generate_simu_index_csv(self):
        #series setup
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/'
        output_file_csv = prefix_write + 'pin_hex_to_honeycomb_klt_2m_gauss.csv'
        index1 = 5813
        #list_simu_index = np.linspace(index1,index1+n_files-1,n_files,dtype=int)
        gauss_epsilon1 = -3
        gauss_epsilon_step = -3
        n_gauss_epsilons = 30
        list_gauss_epsilon = np.linspace(gauss_epsilon1,gauss_epsilon1+gauss_epsilon_step*(n_gauss_epsilons-1),n_gauss_epsilons)
        list_lcr = [0.77,0.78,0.79,0.80]
        temperature = 1
        t = temperature

        
        for seed in range(10):
            simu_index = index1
            for lcr in list_lcr:
                for gauss_epsilon in list_gauss_epsilon:
                    if ((simu_index == index1) & (seed == 0)):
                        record1 = [[index1,seed,list_lcr[0],gauss_epsilon1,t]]
                        column = ['simu_index','seed','lcr','trap_gauss_epsilon','temperature']#,'cn3'
                        record_all = pd.DataFrame(record1)
                        record_all.columns = column
                        simu_index = simu_index + 1
                    else:
                        record1 = pd.DataFrame([[simu_index,seed,lcr,gauss_epsilon,t]],columns=column)
                        record_all = pd.concat([record_all,record1])
                        simu_index = simu_index + 1
        pd.DataFrame.to_csv(record_all,output_file_csv)
    
    def generate_simu_index_csv_5773_5812(self):
        #series setup
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/'
        output_file_csv = prefix_write + 'pin_hex_to_honeycomb_klt_2m_gauss_5773_5812.csv'
        index1 = 5773
        #list_simu_index = np.linspace(index1,index1+n_files-1,n_files,dtype=int)
        gauss_epsilon1 = -3
        gauss_epsilon_step = -3
        n_gauss_epsilons = 20
        list_gauss_epsilon1 = np.linspace(gauss_epsilon1,gauss_epsilon1+gauss_epsilon_step*(n_gauss_epsilons-1),n_gauss_epsilons)
        gauss_epsilon1 = -100
        gauss_epsilon_step = -100
        n_gauss_epsilons = 10
        list_gauss_epsilon2 = np.linspace(gauss_epsilon1,gauss_epsilon1+gauss_epsilon_step*(n_gauss_epsilons-1),n_gauss_epsilons)
        gauss_epsilon1 = -63
        gauss_epsilon_step = -3
        n_gauss_epsilons = 10
        list_gauss_epsilon3 = np.linspace(gauss_epsilon1,gauss_epsilon1+gauss_epsilon_step*(n_gauss_epsilons-1),n_gauss_epsilons)
        list_gauss_epsilon = np.concatenate((list_gauss_epsilon1,list_gauss_epsilon2,list_gauss_epsilon3))
        list_lcr = 0.81
        temperature = 1
        t = temperature

        
        for seed in range(10):
            simu_index = index1
            for gauss_epsilon in list_gauss_epsilon:
                if ((simu_index == index1) & (seed == 0)):
                    record1 = [[index1,seed,list_lcr,list_gauss_epsilon[0],t]]
                    column = ['simu_index','seed','lcr','trap_gauss_epsilon','temperature']#,'cn3'
                    record_all = pd.DataFrame(record1)
                    record_all.columns = column
                    simu_index = simu_index + 1
                else:
                    record1 = pd.DataFrame([[simu_index,seed,list_lcr,gauss_epsilon,t]],columns=column)
                    record_all = pd.concat([record_all,record1])
                    simu_index = simu_index + 1
        pd.DataFrame.to_csv(record_all,output_file_csv)
    
    def record_simu_csv_into_sql(self):
        import opertateOnMysql as osql
        filename_csv = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_pin/honeycomb_pin_scan_k3_60_simu.csv'
        #osql.createTableInMysql('pin_hex_to_honeycomb_klt_2m_gauss','simu_index integer unsigned not null,seed integer unsigned not null,lcr float, trap_gauss_epsilon float, temperature float')
        osql.loadCsvDataToMysql\
        (path_to_file_name=filename_csv,
        table_name='pin_hex_to_honeycomb_klt_2m_gauss')

class simulation_controller_honeycomb_part_traps:
    def __init__(self):
        pass    

    def generate_simu_index_csv_lcr080(self):
        #series setup
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
        output_file_csv = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_5973.csv'#5933,
        index1 = 5973
        #list_simu_index = np.linspace(index1,index1+n_files-1,n_files,dtype=int)
        gauss_epsilon1 = -3
        gauss_epsilon_step = -3
        n_gauss_epsilons = 30
        list_gauss_epsilon1 = np.linspace(gauss_epsilon1,gauss_epsilon1+gauss_epsilon_step*(n_gauss_epsilons-1),n_gauss_epsilons)
        gauss_epsilon1 = -100
        gauss_epsilon_step = -100
        n_gauss_epsilons = 10
        list_gauss_epsilon2 = np.linspace(gauss_epsilon1,gauss_epsilon1+gauss_epsilon_step*(n_gauss_epsilons-1),n_gauss_epsilons)
        list_gauss_epsilon = np.concatenate((list_gauss_epsilon1,list_gauss_epsilon2))
        list_lcr = 0.80
        temperature = 1
        t = temperature

        
        seed=9#for seed in range(10):
        simu_index = index1
        for gauss_epsilon in list_gauss_epsilon:
            if (simu_index == index1):
                record1 = [[index1,seed,list_lcr,list_gauss_epsilon[0],t]]
                column = ['simu_index','seed','lcr','trap_gauss_epsilon','temperature']#,'cn3'
                record_all = pd.DataFrame(record1)
                record_all.columns = column
                simu_index = simu_index + 1
            else:
                record1 = pd.DataFrame([[simu_index,seed,list_lcr,gauss_epsilon,t]],columns=column)
                record_all = pd.concat([record_all,record1])
                simu_index = simu_index + 1
        pd.DataFrame.to_csv(record_all,output_file_csv)
    
    def generate_simu_index_csv_lcr081(self):
        #series setup
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
        index1 = 6013
        output_file_csv = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_'+str(int(index1))+'.csv'#5933,
        #list_simu_index = np.linspace(index1,index1+n_files-1,n_files,dtype=int)
        gauss_epsilon1 = -3
        gauss_epsilon_step = -3
        n_gauss_epsilons = 30
        list_gauss_epsilon1 = np.linspace(gauss_epsilon1,gauss_epsilon1+gauss_epsilon_step*(n_gauss_epsilons-1),n_gauss_epsilons)
        gauss_epsilon1 = -100
        gauss_epsilon_step = -100
        n_gauss_epsilons = 10
        list_gauss_epsilon2 = np.linspace(gauss_epsilon1,gauss_epsilon1+gauss_epsilon_step*(n_gauss_epsilons-1),n_gauss_epsilons)
        list_gauss_epsilon = np.concatenate((list_gauss_epsilon1,list_gauss_epsilon2))
        list_lcr = 0.81
        temperature = 1
        t = temperature

        
        seed=9#for seed in range(10):
        simu_index = index1
        for gauss_epsilon in list_gauss_epsilon:
            if (simu_index == index1):
                record1 = [[index1,seed,list_lcr,list_gauss_epsilon[0],t]]
                column = ['simu_index','seed','lcr','trap_gauss_epsilon','temperature']#,'cn3'
                record_all = pd.DataFrame(record1)
                record_all.columns = column
                simu_index = simu_index + 1
            else:
                record1 = pd.DataFrame([[simu_index,seed,list_lcr,gauss_epsilon,t]],columns=column)
                record_all = pd.concat([record_all,record1])
                simu_index = simu_index + 1
        pd.DataFrame.to_csv(record_all,output_file_csv)
    
    def generate_simu_index_csv_lcr081_brownian(self):
        #series setup
        prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
        index1 = 6013
        output_file_csv = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_b_'+str(int(index1))+'.csv'#5933,
        #list_simu_index = np.linspace(index1,index1+n_files-1,n_files,dtype=int)
        gauss_epsilon1 = -3
        gauss_epsilon_step = -3
        n_gauss_epsilons = 30
        list_gauss_epsilon1 = np.linspace(gauss_epsilon1,gauss_epsilon1+gauss_epsilon_step*(n_gauss_epsilons-1),n_gauss_epsilons)
        gauss_epsilon1 = -100
        gauss_epsilon_step = -100
        n_gauss_epsilons = 10
        list_gauss_epsilon2 = np.linspace(gauss_epsilon1,gauss_epsilon1+gauss_epsilon_step*(n_gauss_epsilons-1),n_gauss_epsilons)
        list_gauss_epsilon = np.concatenate((list_gauss_epsilon1,list_gauss_epsilon2))
        list_lcr = 0.81
        temperature = 1
        t = temperature

        
        seed=8#for seed in range(10):
        simu_index = index1
        for gauss_epsilon in list_gauss_epsilon:
            if (simu_index == index1):
                record1 = [[index1,seed,list_lcr,list_gauss_epsilon[0],t]]
                column = ['simu_index','seed','lcr','trap_gauss_epsilon','temperature']#,'cn3'
                record_all = pd.DataFrame(record1)
                record_all.columns = column
                simu_index = simu_index + 1
            else:
                record1 = pd.DataFrame([[simu_index,seed,list_lcr,gauss_epsilon,t]],columns=column)
                record_all = pd.concat([record_all,record1])
                simu_index = simu_index + 1
        pd.DataFrame.to_csv(record_all,output_file_csv)

    def generate_initial_state_hexagonal_particle_honeycomb_part_trap_scan_csv(self,csv_filename):
        R"""
        Introduction:
                a comprehensive workflow from creating initial state, 
            to record trajectory as a series of .png images, covering a 
            group of scanning k simulation points (simu_index,seed,lcr,k)

            csv format: column = ['simu_index','seed','lcr','trap_gauss_epsilon','temperature']
        Example:
            import symmetry_transformation_v4_3.simulation_controller as sc
            #import symmetry_transformation_v4_3.simulation_core as sco
            sct = sc.simulation_controller_honeycomb_part_traps()
            prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
            file_csv = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_5933.csv'
            sct.generate_initial_state_hexagonal_particle_honeycomb_part_trap_scan_csv(file_csv)
        """
        #series setup
        df = pd.read_csv(csv_filename)
        
        col = df.columns
        print(df.head())
        
        n_size = [16,8]
        a_particle = 3
        list_simu_index = df[col[1]].values
        list_seed = df[col[2]].values
        list_lcr = df[col[3]].values
        list_trap_gauss_epsilon = df[col[4]].values
        list_temperature = df[col[5]].values

        for i in range(len(list_simu_index)):
            #operate simulation
            simu_index = list_simu_index[i]
            seed = list_seed[i]
            lcr = list_lcr[i]
            gauss_epsilon = list_trap_gauss_epsilon[i]
            temperature = list_temperature[i]
            init_gsd_filename = self.generate_initial_gsd(a_particle,n_size,lcr)
            print(simu_index,seed,lcr,gauss_epsilon,temperature,init_gsd_filename)
            sim = sco.simulation_core_traps(simu_index,seed)
            sim.seed=seed
            sim.mode = 'gpu'
            sim.gauss_epsilon = gauss_epsilon#-50
            sim.kT = temperature
            sim.input_file_gsd = init_gsd_filename
            sim.operate_simulation_brownian()
            print(i+1,'/',len(list_simu_index))
        
    def generate_initial_gsd(sel,a_particle,n_size,lcr):
        prefix = "/media/remote/32E2D4CCE2D49607/file_lxt/hoomd-examples_0/"
        #filename_initial_state should be like 
        # 'init_particle_<lattice_type>_<a>_<nx>_<ny>_and_trap_<lattice_type>_<int_lcr>
        # to prevent the mixing of different initial state .gsd file
        str_lcr = str(lcr)#'081'#in case the string like '0.810000000032'
        filename_initial_state = 'init_particle_hex_'+str(a_particle)+'_'+str(n_size[0])+'_'+str(n_size[1])+'_and_trap_honeycomb_part_'+str_lcr
                
        init_gsd_filename = prefix+filename_initial_state+'.gsd'#describe the lcr
         
        #check if the file exists
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
            traps.generate_type3_part(a=a_particle*lcr)
            isg = pg.initial_state_generator()
            isg.set_new_gsd_file_2types(particles,n_size,particle_points,traps,init_gsd_filename,trap_fill_box=True)
            
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
            print('create:'+init_gsd_filename)
        return init_gsd_filename