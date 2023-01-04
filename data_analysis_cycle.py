from genericpath import exists
import numpy 
import matplotlib.pyplot as plt
import points_analysis_2D
import os

def saveIndexPsi3Psi6(start_index,end_index,k1,step,linear_compression_ratio):
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data [| simu_index | HarmonicK | LinearCompressionRatio | Psi3Global | Psi6Global|RandomSeed]
    """
    #start_index=193
    #end_index=205
    diference_index=end_index-start_index+1
    prefix='/home/tplab/Downloads/'
    record=numpy.zeros((diference_index,6))

    for index in numpy.linspace(start_index,end_index,diference_index):
        #index=index0+(i-1.0)
        data_filename=prefix+'index'+str(index.astype(int))
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        png_filename=None#prefix+'index'+str(index)+'Psi3.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=3,plot=True,png_filename=png_filename)
        record[(index-start_index).astype(int),3]=obj_of_simu_index.Psi_k_global_cut_edge

        png_filename=None#prefix+'index'+str(index)+'Psi6.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=6,plot=True,png_filename=png_filename)#plot=True
        record[(index-start_index).astype(int),4]=obj_of_simu_index.Psi_k_global_cut_edge
        #print(obj_of_simu_index.Psi_k_global)

        record[(index-start_index).astype(int),5]=9

    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'kl'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexPsi4Psi6(start_index,end_index,k1,step,linear_compression_ratio):
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data [| simu_index | HarmonicK | LinearCompressionRatio | Psi3Global | Psi6Global]
    """
    #start_index=193
    #end_index=205
    diference_index=end_index-start_index+1
    prefix='/home/tplab/Downloads/'
    record=numpy.zeros((diference_index,5))

    for index in numpy.linspace(start_index,end_index,diference_index):
        #index=index0+(i-1.0)
        data_filename=prefix+'index'+str(index.astype(int))
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        png_filename=prefix+'index'+str(index)+'Psi4.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=4,plot=True,png_filename=png_filename)
        record[(index-start_index).astype(int),3]=obj_of_simu_index.Psi_k_global

        png_filename=prefix+'index'+str(index)+'Psi6.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=6,plot=True,png_filename=png_filename)#plot=True
        record[(index-start_index).astype(int),4]=obj_of_simu_index.Psi_k_global
        #print(obj_of_simu_index.Psi_k_global)
    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'kl4'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexPsi3Ratio(start_index,end_index,k1,step,linear_compression_ratio):#[x]
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data [| simu_index | HarmonicK | LinearCompressionRatio | Psi3Global | Psi6Global]
    """
    #start_index=193
    #end_index=205
    diference_index=end_index-start_index+1
    prefix='/home/tplab/Downloads/'
    record=numpy.zeros((diference_index,5))

    for index in numpy.linspace(start_index,end_index,diference_index):
        #index=index0+(i-1.0)
        data_filename=prefix+'index'+str(index.astype(int))
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        png_filename=prefix+'index'+str(index)+'Psi3.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=3,plot=True,png_filename=png_filename)
        record[(index-start_index).astype(int),3]=obj_of_simu_index.Psi_k_global

        #png_filename=prefix+'index'+str(index)+'Psi3Ratio.png'
        obj_of_simu_index.get_bond_orientational_order_selected(k_set=3)#plot=True
        record[(index-start_index).astype(int),4]=obj_of_simu_index.Psi_k_rate
        #print(obj_of_simu_index.Psi_k_global)
    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'kl'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexCairoCN3456(start_index,end_index,k1,step,linear_compression_ratio,randomseed):
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data 
    [ simu_index | HarmonicK | LinearCompressionRatio | 
        CoordinationNum3Rate | CoordinationNum4Rate | 
        CoordinationNum6Rate |RandomSeed ]
    """
    diference_index=end_index-start_index+1
    prefix='/home/tplab/Downloads/'
    record=numpy.zeros((diference_index,8))

    for index in numpy.linspace(start_index,end_index,diference_index):
        data_filename=prefix+'index'+str(index.astype(int))
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        #get 1st neighbor bonds
        s = "/home/tplab/Downloads/"
        index_name = "index"+str(int(index))+'.png'
        oname1 = s+"bond_hist_"+index_name
        obj_of_simu_index.draw_bond_length_distribution_and_first_minima(lattice_constant=3*linear_compression_ratio,png_filename=oname1)
        oname2 = s+"bond_plot_"+index_name
        obj_of_simu_index.draw_bonds_conditional_bond(check=[0.9, obj_of_simu_index.bond_first_minima_left], png_filename=oname2)

        obj_of_simu_index.get_coordination_number_conditional()
        ccn=obj_of_simu_index.count_coordination_ratio
        record[(index-start_index).astype(int),3]=ccn[3]

        record[(index-start_index).astype(int),4]=ccn[4]

        record[(index-start_index).astype(int),5]=ccn[5]

        record[(index-start_index).astype(int),6]=ccn[6]

        record[(index-start_index).astype(int),7]=randomseed

        
        
    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'klcn3456'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexCairoPsi345(start_index,end_index,k1,step,linear_compression_ratio):#[x]
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data [| simu_index | HarmonicK | LinearCompressionRatio | Psi3Global | Psi6Global]
    """
    #start_index=193
    #end_index=205
    diference_index=end_index-start_index+1
    prefix='/home/tplab/Downloads/'
    record=numpy.zeros((diference_index,6))

    for index in numpy.linspace(start_index,end_index,diference_index):
        #index=index0+(i-1.0)
        data_filename=prefix+'index'+str(index.astype(int))
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        png_filename=prefix+'index'+str(index)+'Psi3.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=3,plot=True,png_filename=png_filename)
        record[(index-start_index).astype(int),3]=obj_of_simu_index.Psi_k_global

        png_filename=prefix+'index'+str(index)+'Psi4.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=4,plot=True,png_filename=png_filename)
        record[(index-start_index).astype(int),4]=obj_of_simu_index.Psi_k_global

        png_filename=prefix+'index'+str(index)+'Psi5.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=5,plot=True,png_filename=png_filename)
        record[(index-start_index).astype(int),5]=obj_of_simu_index.Psi_k_global

        #print(obj_of_simu_index.Psi_k_global)
    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'klp345'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexCN346Seed(start_index,end_index,k1,step,linear_compression_ratio,randomseed):
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data 
    [ simu_index | HarmonicK | LinearCompressionRatio | 
        CoordinationNum3Rate | CoordinationNum4Rate | 
        CoordinationNum6Rate |RandomSeed ]
    """
    diference_index=end_index-start_index+1
    prefix='/home/tplab/Downloads/'
    record=numpy.zeros((diference_index,7))

    for index in numpy.linspace(start_index,end_index,diference_index):
        data_filename=prefix+'index'+str(index.astype(int))
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        
        obj_of_simu_index.get_coordination_number_conditional()
        ccn=obj_of_simu_index.count_coordination_ratio
        record[(index-start_index).astype(int),3]=ccn[3]

        record[(index-start_index).astype(int),4]=ccn[4]

        record[(index-start_index).astype(int),5]=ccn[6]

        record[(index-start_index).astype(int),6]=randomseed

        
        
    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'klcn346'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexCN4CN6SeedPsi6(start_index,end_index,k1,step,linear_compression_ratio,randomseed):
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data 
    [ simu_index | HarmonicK | LinearCompressionRatio | 
    CoordinationNum4Rate | CoordinationNum6Rate | RandomSeed |
      Psi6Global]
    """
    diference_index=end_index-start_index+1
    prefix='/home/tplab/Downloads/'
    record=numpy.zeros((diference_index,7))

    for index in numpy.linspace(start_index,end_index,diference_index):
        data_filename=prefix+'index'+str(index.astype(int))
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        
        obj_of_simu_index.get_coordination_number_conditional()
        ccn=obj_of_simu_index.count_coordination_ratio
        record[(index-start_index).astype(int),3]=ccn[4]

        record[(index-start_index).astype(int),4]=ccn[6]

        record[(index-start_index).astype(int),5]=randomseed

        #png_filename=prefix+'index'+str(index)+'Psi6.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=6)#plot=True
        record[(index-start_index).astype(int),6]=obj_of_simu_index.Psi_k_global_cut_edge
        
        
    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'kl4'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexCN4CN6SeedPsi6Seed(start_index,end_index,k1,step,linear_compression_ratio,randomseed):
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data 
    [ simu_index | HarmonicK | LinearCompressionRatio | 
    CoordinationNum4Rate | CoordinationNum6Rate | RandomSeed |
      Psi6Global]
    """
    diference_index=end_index-start_index+1
    prefix='/home/tplab/Downloads/'
    record=numpy.zeros((diference_index,7))

    for index in numpy.linspace(start_index,end_index,diference_index):
        data_filename=prefix+'index'+str(index.astype(int))+"_"+str(randomseed)
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        
        obj_of_simu_index.get_coordination_number_conditional()
        ccn=obj_of_simu_index.count_coordination_ratio
        record[(index-start_index).astype(int),3]=ccn[4]

        record[(index-start_index).astype(int),4]=ccn[6]

        record[(index-start_index).astype(int),5]=randomseed

        #png_filename=prefix+'index'+str(index)+'Psi6.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=6)#plot=True
        record[(index-start_index).astype(int),6]=obj_of_simu_index.Psi_k_global_cut_edge

        """
        log_prefix='/home/tplab/hoomd-examples_0/'
        str_index=str(int(index))
        str_seed=str(int(randomseed))
        str_index=str_index+"_"+str_seed
        record[(index-start_index).astype(int),7]=log_prefix+'trajectory_auto'+str_index+'.gsd'
        """
        
        
        
    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'kl4'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexCN4CN3SeedPsi6(start_index,end_index,k1,step,linear_compression_ratio,randomseed,account='tplab'):
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data 
    [ simu_index | HarmonicK | LinearCompressionRatio | 
    CoordinationNum4Rate | CoordinationNum3Rate | RandomSeed |
      Psi6Global]
    """
    diference_index=end_index-start_index+1
    prefix='/home/'+account+'/Downloads/'
    record=numpy.zeros((diference_index,7))

    for index in numpy.linspace(start_index,end_index,diference_index):
        data_filename=prefix+'index'+str(index.astype(int))+"_"+str(randomseed)
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        
        obj_of_simu_index.get_coordination_number_conditional()
        ccn=obj_of_simu_index.count_coordination_ratio
        record[(index-start_index).astype(int),3]=ccn[4]

        record[(index-start_index).astype(int),4]=ccn[3]

        record[(index-start_index).astype(int),5]=randomseed

        #png_filename=prefix+'index'+str(index)+'Psi6.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=6)#plot=True
        record[(index-start_index).astype(int),6]=obj_of_simu_index.Psi_k_global_cut_edge

        """
        log_prefix='/home/tplab/hoomd-examples_0/'
        str_index=str(int(index))
        str_seed=str(int(randomseed))
        str_index=str_index+"_"+str_seed
        record[(index-start_index).astype(int),7]=log_prefix+'trajectory_auto'+str_index+'.gsd'
        """

    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'kl4'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexklTCN3CN4Seed(start_index,end_index,k1,step,linear_compression_ratio,kT,randomseed,account='tplab'):
    R"""
    This function will save a txt file named 'start index - end index klt', which contains
    n rows of data 
     | SimuIndex | HarmonicK | LinearCompressionRatio | kT | 
    CoordinationNum3Rate | CoordinationNum4Rate | RandomSeed | 

    CHECK: [v]
    """
    diference_index=end_index-start_index+1
    prefix='/home/'+account+'/Downloads/'
    record=numpy.zeros((diference_index,7))

    for index in numpy.linspace(start_index,end_index,diference_index):
        data_filename=prefix+'index'+str(index.astype(int))+"_"+str(randomseed)
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        record[(index-start_index).astype(int),3]=kT

        obj_of_simu_index.get_coordination_number_conditional()
        ccn=obj_of_simu_index.count_coordination_ratio
        record[(index-start_index).astype(int),4]=ccn[3]

        record[(index-start_index).astype(int),5]=ccn[4]

        record[(index-start_index).astype(int),6]=randomseed

    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'klt'
    numpy.savetxt(save_file_name,record)
    return save_file_name

def saveIndexklTPsi36Seed(start_index,end_index,k1,step,linear_compression_ratio,kT,randomseed,account='tplab'):
    R"""
    This function will save a txt file named 'start index - end index klt', which contains
    n rows of data 
    |simu_index | HarmonicK | LinearCompressionRatio | kT |
      Psi3Global | Psi6Global | RandomSeed | 

    CHECK: [v]
    """
    diference_index=end_index-start_index+1
    prefix='/home/'+account+'/Downloads/'
    record=numpy.zeros((diference_index,7))

    for index in numpy.linspace(start_index,end_index,diference_index):
        data_filename=prefix+'index'+str(index.astype(int))+"_"+str(randomseed)
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        record[(index-start_index).astype(int),3]=kT

        obj_of_simu_index.get_bond_orientational_order(k_set=3)
        record[(index-start_index).astype(int),4]=obj_of_simu_index.Psi_k_global_cut_edge

        obj_of_simu_index.get_bond_orientational_order(k_set=6)
        record[(index-start_index).astype(int),5]=obj_of_simu_index.Psi_k_global_cut_edge

        record[(index-start_index).astype(int),6]=randomseed

    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'klt'
    numpy.savetxt(save_file_name,record)
    return save_file_name

def saveIndexkTCN4CN3SeedPsi6(start_index,end_index,kt1,step,linear_compression_ratio,kset,randomseed):
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data 
    [ simu_index | HarmonicK | LinearCompressionRatio | kT |
    CoordinationNum4Rate | CoordinationNum3Rate | RandomSeed |
      Psi6Global]
    """
    diference_index=end_index-start_index+1
    prefix='/home/tplab/Downloads/'
    record=numpy.zeros((diference_index,8))

    for index in numpy.linspace(start_index,end_index,diference_index):
        data_filename=prefix+'index'+str(index.astype(int))+"_"+str(randomseed)
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=kset

        record[(index-start_index).astype(int),2]=linear_compression_ratio
        
        record[(index-start_index).astype(int),3]=kt1+step*(index-start_index)
        
        obj_of_simu_index.get_coordination_number_conditional()
        ccn=obj_of_simu_index.count_coordination_ratio
        record[(index-start_index).astype(int),4]=ccn[4]

        record[(index-start_index).astype(int),5]=ccn[3]

        record[(index-start_index).astype(int),6]=randomseed

        #png_filename=prefix+'index'+str(index)+'Psi6.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=6)#plot=True
        record[(index-start_index).astype(int),7]=obj_of_simu_index.Psi_k_global_cut_edge

        """
        log_prefix='/home/tplab/hoomd-examples_0/'
        str_index=str(int(index))
        str_seed=str(int(randomseed))
        str_index=str_index+"_"+str_seed
        record[(index-start_index).astype(int),7]=log_prefix+'trajectory_auto'+str_index+'.gsd'
        """

    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'klt43'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexkTCN4CN3SeedPsi6_cool(start_index,end_index,kt1,linear_compression_ratio,kset,randomseed):
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data 
    [ simu_index | HarmonicK | LinearCompressionRatio | kT |
    CoordinationNum4Rate | CoordinationNum3Rate | RandomSeed |
      Psi6Global]
    """
    diference_index=end_index-start_index+1
    prefix='/home/tplab/Downloads/'
    record=numpy.zeros((diference_index,8))

    for index in numpy.linspace(start_index,end_index,diference_index):
        data_filename=prefix+'index'+str(index.astype(int))+"_"+str(randomseed)
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=kset

        record[(index-start_index).astype(int),2]=linear_compression_ratio
        
        record[(index-start_index).astype(int),3]=kt1#+step*(index-start_index)
        
        obj_of_simu_index.get_coordination_number_conditional()
        ccn=obj_of_simu_index.count_coordination_ratio
        record[(index-start_index).astype(int),4]=ccn[4]

        record[(index-start_index).astype(int),5]=ccn[3]

        record[(index-start_index).astype(int),6]=randomseed

        #png_filename=prefix+'index'+str(index)+'Psi6.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=6)#plot=True
        record[(index-start_index).astype(int),7]=obj_of_simu_index.Psi_k_global_cut_edge

        """
        log_prefix='/home/tplab/hoomd-examples_0/'
        str_index=str(int(index))
        str_seed=str(int(randomseed))
        str_index=str_index+"_"+str_seed
        record[(index-start_index).astype(int),7]=log_prefix+'trajectory_auto'+str_index+'.gsd'
        """

    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'klt43'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexkTPsi36Seed(account,start_index,end_index,kt1,step,linear_compression_ratio,kset,randomseed):
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data 
    | simu_index | HarmonicK | LinearCompressionRatio | kT |
      Psi3Global | Psi6Global | RandomSeed |
    """
    diference_index=end_index-start_index+1
    prefix='/home/'+account+'/Downloads/'
    record=numpy.zeros((diference_index,7))

    for index in numpy.linspace(start_index,end_index,diference_index):
        data_filename=prefix+'index'+str(index.astype(int))+"_"+str(randomseed)
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=kset

        record[(index-start_index).astype(int),2]=linear_compression_ratio
        
        record[(index-start_index).astype(int),3]=kt1+step*(index-start_index)
        
        obj_of_simu_index.get_bond_orientational_order(k_set=3)
        record[(index-start_index).astype(int),4]=obj_of_simu_index.Psi_k_global_cut_edge

        obj_of_simu_index.get_bond_orientational_order(k_set=6)
        record[(index-start_index).astype(int),5]=obj_of_simu_index.Psi_k_global_cut_edge

        record[(index-start_index).astype(int),6]=randomseed

    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'klt36'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexPsi3Psi6R(account,start_index,end_index,k1,step,linear_compression_ratio,randomseed):
    R"""
    This function will save a txt file named 'start index - end index kl', which contains
    n rows of data [| simu_index | HarmonicK | LinearCompressionRatio | 
                      Psi3Global | Psi6Global| RandomSeed]
    """
    #start_index=193
    #end_index=205
    diference_index=end_index-start_index+1
    prefix='/home/'+account+'/Downloads/'
    record=numpy.zeros((diference_index,6))

    for index in numpy.linspace(start_index,end_index,diference_index):
        #index=index0+(i-1.0)
        data_filename=prefix+'index'+str(index.astype(int))+'_'+str(randomseed)
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)
        
        
        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=k1+step*(index-start_index)

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        png_filename=prefix+'index'+str(index)+'Psi3.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=3)#,plot=True,png_filename=png_filename
        record[(index-start_index).astype(int),3]=obj_of_simu_index.Psi_k_global_cut_edge

        png_filename=prefix+'index'+str(index)+'Psi6.png'
        obj_of_simu_index.get_bond_orientational_order(k_set=6)#,plot=True,png_filename=png_filename
        record[(index-start_index).astype(int),4]=obj_of_simu_index.Psi_k_global_cut_edge
        #print(obj_of_simu_index.Psi_k_global)

        record[(index-start_index).astype(int),5]=randomseed
        
    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'kl'
    numpy.savetxt(save_file_name,record)
    #print(record)
    return save_file_name

def saveIndexPressure(start_index,end_index,k1,step,linear_compression_ratio,seed_set):
    R"""
    This function will save a txt file named 'start index - end index Tlp', which contains
    n rows of data [| simu_index| KBT| LinearCompressionRatio| Pressure| Psi6Global| RandomSeed]
    [x]maybe I should check area of voronoi cell and the square of LCR
    """
    diference_index=end_index-start_index+1
    prefix='/home/tplab/Downloads/'
    record=numpy.zeros((diference_index,6))
    
    #error check part
    error_record=numpy.zeros((20,2))
    error_count=0
    
    #simulation setup
    for index in numpy.linspace(start_index,end_index,diference_index):
        #set parameters
        #linear_compression_ratio=linear_compression_ratio1+step*(i-1.0)
        kset=k1+step*(index-start_index)
        str_index=str(index.astype(int))
        log_prefix='/home/tplab/hoomd-examples_0/'
        file_log=log_prefix+'log-output_auto'+str_index+'.log'
        #file_gsd=log_prefix+'trajectory_auto'+str_index+'.gsd'
        prefix='/home/tplab/Downloads/'
        #png_filename=prefix+'index'+str_index+'.png'

        KBT,pressure=get_KBT_pressure(file_log)
        #error check part
        check=abs((KBT-kset)/KBT)
        if check>0.015:
            print("error:temperature fluctuating!")
            error_record[error_count,0]=index
            error_record[error_count,1]=check
            error_count=error_count+1
            
        

        data_filename=prefix+'index'+str_index
        obj_of_simu_index = points_analysis_2D.static_points_analysis_2d(filename=data_filename)

        record[(index-start_index).astype(int),0]=index

        record[(index-start_index).astype(int),1]=kset#KBT

        record[(index-start_index).astype(int),2]=linear_compression_ratio

        record[(index-start_index).astype(int),3]=pressure
        
        obj_of_simu_index.get_bond_orientational_order(k_set=6)#plot=True
        record[(index-start_index).astype(int),4]=obj_of_simu_index.Psi_k_global_cut_edge#20221202 updated

        record[(index-start_index).astype(int),5]=seed_set

        #print(obj_of_simu_index.Psi_k_global)
    save_file_name=prefix+str(start_index)+'-'+str(end_index)+'Tlp'
    numpy.savetxt(save_file_name,record)
    #error check part
    error_file_name=prefix+str(start_index)+'-'+str(end_index)+'Terror'
    numpy.savetxt(error_file_name,error_record)
    #print(record)
    return save_file_name



def get_KBT_pressure(file_log):
    data = numpy.genfromtxt(fname=file_log, skip_header=True)
    average_KBT=numpy.average(data[300:-1,2])
    average_pressure=numpy.average(data[300:-1,3])

    #check temperature
    return average_KBT,average_pressure


def rearrange_data():
    R"""
    Purpose:
        merge multiple txt file [index psi_k_global psi_k_rate] into one,
        and at last get 
        [1     2 3                        4            5         ]
        [index k linear_compression_ratio psi_k_global psi_k_rate]
    """
    lcr_list=numpy.linspace(0.77,0.88,12)
    #control pointers, not data.
    #linear_compression_ratio2=0.80
    fn1='/home/tplab/Downloads/91-102'
    #fn2='/home/tplab/Downloads/103-132'
    #fn3='/home/tplab/Downloads/133-192'
    data1=readdata(fn1)
    #data2=readdata(fn2)
    #data3=readdata(fn3)
    sp1=numpy.shape(data1)
    #sp2=numpy.shape(data2)
    #sp3=numpy.shape(data3)

    data_ex=numpy.zeros((sp1[0],sp1[1]+2))#+sp2[0]+sp3[0]
    data_ex[0:sp1[0],0]=data1[:,0]
    data_ex[0:sp1[0],1]=800
    data_ex[0:sp1[0],2]=lcr_list
    data_ex[0:sp1[0],3:5]=data1[:,1:3]

    print(data_ex)
    #save_filename=fn1+'kl'
    save_filename='/home/tplab/Downloads/91-102kl'
    numpy.savetxt(save_filename,data_ex)
    '''
    #merge two 5-column list
    fn1='/home/tplab/Downloads/103-132kl'
    fn2='/home/tplab/Downloads/133-192kl'
    data1=readdata(fn1)
    data2=readdata(fn2)
    sp1=numpy.shape(data1)
    sp2=numpy.shape(data2)
    data_ex=numpy.zeros((sp1[0]+sp2[0],sp1[1]))
    data_ex[0:sp1[0]]=data1[:]
    data_ex[sp1[0]:(sp1[0]+sp2[0])]=data2[:]
    '''
    '''
    #merge two 3-column list
    data_ex[0:30,1]=k_list
    data_ex[30:60,1]=k_list
    data_ex[0:30,2]=linear_compression_ratio1
    data_ex[30:60,2]=linear_compression_ratio2
    '''
    '''
    #add kl for 3-column list, generate 5-column list
    lcr_list=numpy.linspace(0.77,0.88,12)
    fn1='/home/tplab/Downloads/91-102'
    data1=readdata(fn1)
    sp1=numpy.shape(data1)
    data_ex=numpy.zeros((sp1[0],sp1[1]+2))
    data_ex[0:sp1[0],0]=data1[:,0]
    data_ex[0:sp1[0],1]=800
    data_ex[0:sp1[0],2]=lcr_list
    data_ex[0:sp1[0],3:5]=data1[:,1:3]
    '''

    


def readdata(fn):
    data=numpy.loadtxt(fn)
    return data

def plotHarmonicKAndPsi():
    data=readdata("/home/tplab/Downloads/193-205kl")
    data_psi3=data[:,[1,3]]
    data_psi6=data[:,[1,4]]
    plt.figure()
    p3=plt.scatter(data_psi3[:,0],data_psi3[:,1],c='r')
    p6=plt.scatter(data_psi6[:,0],data_psi6[:,1],c='b')
    plt.legend(handles=[p3,p6],labels=['Psi3','Psi6'],loc='best')
    plt.xlabel('Harmonic_K');
    plt.ylabel('Psi_3 & Psi6');
    plt.show()

def phase_diagram(select=True):
    if select:
        tt='Psi_3_ratio'#tt='Psi_3_ratio'
        prefix='/home/tplab/Downloads/'
        name='103-192kl'
        filename=prefix+name
        data=readdata(filename)

        k_unit=0.5*data[:,1]*numpy.multiply(1,1)
        k_unit=k_unit.astype(int)
        plt.figure()
        plt.scatter(data[:,2]*100,k_unit,c=data[:,4])
        plt.title(tt+'(>0.9)(E_yukawa=233,substrate=400)')
        plt.colorbar()
        plt.xlabel('linear_compression_ratio/%')
        plt.ylabel('depth of trap/kT')
        png_filename=filename+tt
        plt.savefig(png_filename)
        plt.show()
    else:
        tt='Psi_3'#tt='Psi_3_ratio'
        prefix='/home/tplab/Downloads/'
        name='103-192kl'
        filename=prefix+name
        data=readdata(filename)

        k_unit=0.5*data[:,1]*numpy.multiply(1,1)
        k_unit=k_unit.astype(int)
        plt.figure()
        plt.scatter(data[:,2]*100,k_unit,c=data[:,3])
        plt.title(tt+'(E_yukawa=233,substrate=400)')
        plt.colorbar()
        plt.xlabel('linear_compression_ratio/%')
        plt.ylabel('depth of trap/kT')
        png_filename=filename+tt
        plt.savefig(png_filename)
        plt.show()

def phase_diagram_line(select=True):
    prefix='/home/tplab/Downloads/'
    name='91-102kl'
    filename=prefix+name
    data=readdata(filename)

    #k_unit=0.5*data[:,1]*numpy.multiply(1,1)
    #k_unit=k_unit.astype(int)

    if select:
        tt='Psi_3_ratio'#tt='Psi_3_ratio'
        plt.figure()
        plt.plot(data[:,2]*100,data[:,4])
        plt.title(tt+'(>0.9)(E_yukawa=233,)')
        plt.xlabel('linear_compression_ratio/%')
        plt.ylabel(tt)
        png_filename=filename+tt
        plt.savefig(png_filename)
        plt.show()
    else:
        tt='Psi_3'#tt='Psi_3_ratio'
        plt.figure()
        plt.plot(data[:,2]*100,data[:,3])
        plt.title(tt+'(E_yukawa=233)')
        plt.xlabel('linear_compression_ratio/%')
        plt.ylabel(tt)
        png_filename=filename+tt
        plt.savefig(png_filename)
        plt.show()

def save_from_gsd(simu_index=None,seed=None,frame_cut=0,
                    trajectory=False,
                    save_result_txt=False,
                    displacement_field=False,
                    final_cut=False,
                    psik=False,
                    psik_plot=None,
                    neighbor_cloud=False,
                    coordination_number=False,
                    coordination_number3_plot=False,
                    bond_plot=False,bond_plot_gr=False,nb_change=False,show_traps=False,trap_filename=None,trap_lcr=None,
                    gr=False,
                    sk=False,log_sk=False,
                    msd=False,single_particle=False,
                    account='tplab'):
    R"""
    Introduction:
        Read a gsd file and save a series of analyzed results as follow.
        trajectory: 
        displacement_field:
        psik: global psi_k vs time. 
        neighbor_cloud:
        coordination_number:
        coordination_number3_plot:
        final_cut: true to proceed the last frame only.
        
    Format:
        [Psi_3_global,Psi_6_global]

    example:
        import data_analysis_cycle as da
        i=3003
        while i<3062:
            da.save_from_gsd(simu_index=i,seed=9,
                                bond_plot =True,
                                show_traps=True,
                                trap_filename="/home/tplab/hoomd-examples_0/testkagome_part3-11-6",
                                trap_lcr=0.86)
            i=i+1
    example2:
        i=1233
        while i<1243:
            da.save_from_gsd(simu_index=i,final_cut=True,
                        bond_plot =True,
                        gr=True,
                        sk=True)#,seed=9
            i+=1
    example3:
        import data_analysis_cycle as da
        da.save_from_gsd(simu_index=4728,seed=9,
        coordination_number=True,
        bond_plot=True,
        show_traps=True,
        trap_filename='/home/tplab/hoomd-examples_0/testkagome3-11-6',
        trap_lcr=0.853,
        account='remote',)
    """
    import freud
    prefix='/home/'+account+'/Downloads/'#'/home/tplab/Downloads/'
    log_prefix='/home/'+account+'/hoomd-examples_0/'#'/home/tplab/hoomd-examples_0/'
    #load time steps
    if seed is None:
        str_index=str(int(simu_index))
        gsd_data = points_analysis_2D.proceed_gsd_file(simu_index=simu_index)
    else:
        str_index=str(int(simu_index))+'_'+str(seed)
        file_gsd = log_prefix+'trajectory_auto'+str_index+'.gsd'#+'_'+str(seed)
        gsd_data = points_analysis_2D.proceed_gsd_file(filename_gsd_seed=file_gsd,account=account)
        
    file_log=log_prefix+'log-output_auto'+str_index+'.log'#+'_'+str(seed)
    log_data = numpy.genfromtxt(fname=file_log, skip_header=True)
    time_steps = log_data[:,0]

    
    #print(gsd_data.num_of_frames)
    #gsd_data.trajectory.
    if msd:
        #print('this function msd is fault!')
        #pass
        #load particle trajectory [N_frames,N_particles,system_dimension]
        
        """
        iframes = 0
        nframes=gsd_data.num_of_frames
        Np_edge_cut=numpy.shape(gsd_data.edge_cut_positions_list)
        pos_list = numpy.zeros([nframes,Np_edge_cut[1],3])#gsd_data.trajectory[0].particles.N,
        while iframes < nframes:
            af = gsd_data.trajectory.read_frame(iframes)
            pos_list[iframes] = af.particles.position[gsd_data.edge_cut_positions_list]
            iframes = iframes + 1
        #select particles never walking across boundaries.
        cut = points_analysis_2D.select_particle_in_box(pos_list)
        particles_reasonable=cut.compute()
        msds = freud.msd.MSD(af.configuration.box)#the class is fault,,'direct'
        msds.compute(positions=pos_list[:,particles_reasonable,:])#,images=pos_list
        #print(msds.msd)#print(msds.particle_msd)
        
        """
        iframes = 0
        nframes=gsd_data.num_of_frames
        af = gsd_data.trajectory.read_frame(iframes)
        pos_list = numpy.zeros([nframes,af.particles.N,3])#gsd_data.trajectory[0].particles.N,
        while iframes < nframes:
            af = gsd_data.trajectory.read_frame(iframes)
            pos_list[iframes] = af.particles.position
            iframes = iframes + 1
        msds = freud.msd.MSD(af.configuration.box)#the class is fault,,'direct'
        msds.compute(positions=pos_list)#,images=pos_list
        #print(msds.msd)#print(msds.particle_msd)
        # Plot the MSD
        
        plt.figure()
        plt.plot(msds.msd)
        plt.title("Mean Squared Displacement")
        plt.xlabel("$t$")
        plt.ylabel("MSD$(t)$")
        png_filename = prefix +'msd_'+'index'+str_index+'.png'
        plt.savefig(png_filename)#png_filename
        plt.close()
        
        if single_particle:
            plt.figure()
            lenth=numpy.shape(msds.particle_msd[0,:])
            for i in range(lenth[0]):
                plt.semilogy(msds.particle_msd[10:-20,i],base=10)#,msds.particle_msd[:,1],msds.particle_msd[:,2])
            plt.title("Mean Squared Displacement")
            plt.xlabel("$t$")
            plt.ylabel("MSD$(t)$")
            #plt.set_yscale("log",base=2)
            png_filename = prefix +'msd_'+'index'+str_index+'_single.png'
            print(png_filename)
            plt.savefig(png_filename)#png_filename
            plt.close()
            
    
    if trajectory:
        gsd_data.get_trajectory()
    
    for i in range(gsd_data.num_of_frames):
        if final_cut:
            i = gsd_data.num_of_frames-1#i=9#!!! 23
        
        a_frame = points_analysis_2D.static_points_analysis_2d(points=gsd_data.read_a_frame(i))#hide_figure=False
        snap = gsd_data.trajectory.read_frame(i)
        
        if save_result_txt:
            result_filename=prefix+'index'+str_index 
            points=snap.particles.position[:]#temp
            numpy.savetxt(result_filename,points)#temp


        if displacement_field:
            png_filename1 = prefix +'Displacement_Field_xy_'+'index'+str_index+'_'+str(int(i))+'.png'
            gsd_data.get_displacement_field_xy(i,plot=True,png_filename=png_filename1)#
            png_filename2 = prefix +'Displacement_Field_hist_log_'+'index'+str_index+'_'+str(int(i))+'.png'
            gsd_data.get_displacement_field_distribution(i,log_mode=True,png_filename=png_filename2)
            png_filename3 = prefix +'Displacement_Field_hist_'+'index'+str_index+'_'+str(int(i))+'.png'
            gsd_data.get_displacement_field_distribution(i,png_filename=png_filename3)
            
        if psik:
            if not "record_psik" in locals():#check if the variable exists
                #load Psi_k s
                record_psik = numpy.zeros((gsd_data.num_of_frames,3))#[time_steps,psi3,psi6]
                record_psik[:,0] = time_steps#[0:25]*20
            a_frame.get_bond_orientational_order(k_set=3)
            record_psik[i,1] = a_frame.Psi_k_global_cut_edge
            a_frame.get_bond_orientational_order(k_set=6)
            record_psik[i,2] = a_frame.Psi_k_global_cut_edge

        if  not psik_plot is None:
            png_filename_psik = prefix +'bond_orientational_order_'+str(int(psik_plot))+'_'+'index'+str_index+'_'+str(int(i))+'.png'
            a_frame.get_bond_orientational_order(k_set=psik_plot,plot=True,png_filename=png_filename_psik)

        if neighbor_cloud:
            folder_name=prefix+"record_"+str_index#+"/"
            #check if the folder exists
            isExists=os.path.exists(folder_name)
            if isExists:
                pass
            else:
                os.makedirs(folder_name)
            png_filename = folder_name+"/"+'neighbor_cloud_1st_minima_index'+str_index+'_'+str(int(i))+'.png'
            #a_frame.get_neighbor_cloud(png_filename=png_filename)
            a_frame.get_neighbor_cloud_method_1st_minima_bond(png_filename=png_filename)

        if coordination_number:
            R"""
            CN0 % should be 0 for all the particles must be linked by bond.
            CN1 % is likely to be edge?
            CN2 % in body(edge-cutted) shows the mechanical unstability
            CN3 % shows the proportion of honeycomb.
            CN4 % shows the proportion of kagome.
            CN6 % shows the proportion of hexagonal.
            CN5/7 % shows the proportion of disclination.
            """
            #print('index '+str(i))
            #print(snap.particles.position[137])
            a_frame.get_coordination_number_conditional()#cut edge to remove CN012
            ccn = a_frame.count_coordination_ratio#[time_steps,psi3,psi6]
            ccn = numpy.transpose(ccn)
            if not "record_cn" in locals():#check if the variable exists
                #load CN_k s
                record_cn = numpy.zeros((gsd_data.num_of_frames,numpy.shape(ccn)[1]+1))
                record_cn[:,0] = time_steps#range(10)##gsd frame is different from log frame for period set 100 vs 2e3
            #print(numpy.shape(ccn)[1])
            record_cn[i,1:numpy.shape(ccn)[1]+1] = ccn#[0:numpy.shape(ccn)[1]-1]
        
        if coordination_number3_plot:
            pass

        if bond_plot:
            if final_cut:
                #bond_plot+trap_plot
                png_filename1 = prefix +'bond_hist_index'+str_index+'_'+str(int(i))+'.png'
                png_filename2 = prefix +'bond_plot_1st_minima_index'+str_index+'_'+str(int(i))+'.png'
            else:
                folder_name=prefix+"record_"+str_index#+"/"
                #check if the folder exists
                isExists=os.path.exists(folder_name)
                if isExists:
                    pass
                else:
                    os.makedirs(folder_name)
                #bond_plot+trap_plot
                png_filename1 = folder_name+"/" +'bond_hist_index'+str_index+'_'+str(int(i))+'.png'
                png_filename2 = folder_name+"/" +'bond_plot_1st_minima_index'+str_index+'_'+str(int(i))+'.png'
            
            a_frame.get_first_minima_bond_length_distribution(lattice_constant=3)#,png_filename=png_filename1
            a_frame.draw_bonds_conditional_bond(check=[0.4, a_frame.bond_first_minima_left], png_filename=png_filename2,
                                            show_traps=show_traps,LinearCompressionRatio=trap_lcr,trap_filename=trap_filename)
        
        if bond_plot_gr:
            if final_cut:
                #bond_plot+trap_plot
                png_filename1 = prefix +'bond_gr_index'+str_index+'_'+str(int(i))+'.png'
                png_filename2 = prefix +'bond_plot_gr_1st_minima_index'+str_index+'_'+str(int(i))+'.png'
            else:
                folder_name=prefix+"record_"+str_index#+"/"
                #check if the folder exists
                isExists=os.path.exists(folder_name)
                if isExists:
                    pass
                else:
                    os.makedirs(folder_name)
                #bond_plot+trap_plot
                #png_filename1 = folder_name+"/" +'bond_hist_index'+str_index+'_'+str(int(i))+'.png'
                png_filename2 = folder_name+"/" +'bond_plot_1st_minima_index'+str_index+'_'+str(int(i))+'.png'
            
            rdf = freud.density.RDF(bins=150, r_max=15.0,r_min=1.0)#
            rdf.compute(system=snap)
            a_frame.draw_radial_distribution_function_and_first_minima(rdf,lattice_constant=3,png_filename=png_filename1)#
            a_frame.draw_bonds_conditional_bond(check=[0.4, a_frame.bond_first_minima_left], png_filename=png_filename2,
                                            show_traps=show_traps,LinearCompressionRatio=trap_lcr,trap_filename=trap_filename)
            """
            rdf = freud.density.RDF(bins=150, r_max=15.0)#
            rdf.compute(system=snap)
            #print(rdf.bin_centers) print(rdf.bin_counts)
            rdf.plot()
            fig_type = 'gr'
            data_filename=prefix+fig_type+'_index'+str_index+'_'+str(int(i))+'.png'
            plt.savefig(data_filename)
            plt.close()
            """
            """
            #checked right
            import gsd.hoomd
            import freud
            traj = gsd.hoomd.open('/home/tplab/hoomd-examples_0/trajectory_auto5208_9.gsd')
            rdf = freud.density.RDF(bins=50,r_max=10)
            rdf.compute(system=traj[-1])
            r =rdf.bin_centers
            y = rdf.rdf
            import matplotlib.pyplot as plt
            fig,ax = plt.subplots()
            rdf.plot(ax=ax)
            plt.savefig('/home/tplab/Downloads/gr.png')
            """
    
        if sk:
            sk = freud.diffraction.DiffractionPattern()
            sk.compute(system=snap)
            if log_sk:
                fig_type = 'log_sk'        
                data_filename=prefix+fig_type+'_index'+str_index+'_'+str(int(i))+'.png'
                ax = sk.plot(vmin=0.01,vmax=1)
            else:
                fig_type = 'sk'        
                data_filename=prefix+fig_type+'_index'+str_index+'_'+str(int(i))+'.png'
                fig,ax = plt.subplots()
                im = ax.pcolormesh(sk.k_values,sk.k_values,sk.diffraction,cmap='afmhot')#im = 
                #ax.colorbar().remove()
                fig.colorbar(im)
                ax.axis('equal')
            """
            #method1
            ax = sk.plot()
            ax.pcolormesh(sk.k_values,sk.k_values,sk.diffraction,cmap='afmhot')

            #method2
            fig,ax = plt.subplots()
            #print(sk.k_values)
            X, Y = numpy.meshgrid(sk.k_values, sk.k_values)
            im = ax.pcolormesh(X,Y,sk.diffraction,cmap='summer')#'afmhot' im = 
            #https://matplotlib.org/stable/tutorials/colors/colormaps.html
            # ?highlight=afmhot
            #ax.colorbar().remove()
            fig.colorbar(im)
            ax.axis('equal')
            """
            plt.savefig(data_filename)
            plt.close()
            #ax.pcolormesh(X, Y, Z,cmap="plasma",)
            """
            maybe that the loglog map is not suitable for my sk.
            linear colorbar is the right choice
            """
        if final_cut:
            break
    if psik:
        plt.figure()
        plt.plot(record_psik[:,0],record_psik[:,1],label='Psi_3')#psi3
        plt.plot(record_psik[:,0],record_psik[:,2],label='Psi_6')#psi6
        plt.legend()
        plt.title('Psi_3 VS Psi_6 '+'index:'+str_index)
        plt.xlabel('time(steps)')
        plt.ylabel('Psi_k(1)')
        #plt.show()
        png_filename = prefix +'T_VS_Psi_k_'+'index'+str_index+'.png'
        plt.savefig(png_filename)
        plt.close()
    if coordination_number:
        plt.figure()
        if frame_cut == 0:#frame_cut is set to abstract a part of the process to watch in detail
            #plt.plot(record_cn[:,0],record_cn[:,1],label='CN_0')
            #plt.plot(record_cn[:,0],record_cn[:,2],label='CN_1')
            #plt.plot(record_cn[:,0],record_cn[:,3],label='CN_2')
            plt.plot(record_cn[:,0],record_cn[:,4],label='CN_3')
            plt.plot(record_cn[:,0],record_cn[:,5],label='CN_4')
            plt.plot(record_cn[:,0],record_cn[:,6],label='CN_5')
            plt.plot(record_cn[:,0],record_cn[:,7],label='CN_6')
            plt.plot(record_cn[:,0],record_cn[:,8],label='CN_7')
            #plt.plot(record_cn[:,0],record_cn[:,9],label='CN_8')
            #plt.plot(record_cn[:,0],record_cn[:,-1],label='CN_9')
            png_filename = prefix +'T_VS_CN_k'+'index'+str_index+'egcut'+'.png'
        else:
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,1],label='CN_0')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,2],label='CN_1')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,3],label='CN_2')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,4],label='CN_3')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,5],label='CN_4')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,6],label='CN_5')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,7],label='CN_6')
            plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,8],label='CN_7')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,9],label='CN_8')
            #plt.plot(record_cn[0:frame_cut,0],record_cn[0:frame_cut,-1],label='CN_9')
            png_filename = prefix +'T_VS_CN_k_tcut'+'index'+str_index+'egcut'+'.png'
        plt.legend()
        plt.title('CN_k '+'index:'+str_index)
        plt.xlabel('time(steps)')
        plt.ylabel('CN_k(1)')
        #plt.show()
        plt.savefig(png_filename)
        record_filename = prefix +'T_VS_CN_k_cut'+'index'+str_index+'.txt'
        numpy.savetxt(record_filename,record_cn)
        plt.close()


class analyze:
    def __init__(self):
        pass
    
    def get_msd(self,pixel_to_um=3.0/32.0,um_to_sigma=1.0/2.0):
        from points_analysis_2D import dynamic_points_analysis_2d
        import numpy 
        txyz_npy_filename = self.path_to_results+'/'+'txyz_stable'
        traj = numpy.load(txyz_npy_filename)
        traj_um = traj*pixel_to_um# pixel to um
        traj_sigma = traj_um*um_to_sigma# um to sigma
        pa = dynamic_points_analysis_2d(traj_sigma,mode='exp')
        ts_id_dxy, average_1st_bond_length = pa.compute_nearest_neighbor_displacements(unit='sigma')
        import pandas as pd
        ts_id_dxy = pd.read_csv('ts_id_dxy.csv')
        ts_id_dxy['z'] = 0
        pa = dynamic_points_analysis_2d(ts_id_dxy,mode='exp')
        txyz_ids_stable = pa.compute_nearest_neighbor_displacements_stable(pa.txyz_stable)
        pa = dynamic_points_analysis_2d(txyz_ids_stable,mode='exp')
        pa.compute_atmsd_t_chips(0.9)
        time_log = self.path_to_results+'DefaultVideo_5.txt'
        time_log = numpy.loadtxt(time_log)
        pa.plot_lindemann_msd(pa.record_msd,average_1st_bond_length,time_log)
        print('average_1st_bond_length\n',average_1st_bond_length)
        
    
