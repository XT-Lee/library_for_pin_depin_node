def workdlow_simu_to_mysql(index1,lcr):
    #set parameters
   
    #index1=206
    #lcr=0.79
    #print("thread 1 is running\n")


    k1=0.0
    stp=10.0

    #get simulation results
    import symmetry_transformation_auto_honeycomb as sa_h
    end_index=sa_h.workflow(index1=index1,k1=k1,step=stp,k_end=90.0,linear_compression_ratio=lcr)
    'get file index123'

    #get analyzed data
    import data_analysis_cycle as da
    filename_kl=da.saveIndexPsi(start_index=index1,end_index=end_index,k1=k1,step=stp,linear_compression_ratio=lcr)
    #filename_kl=da.saveIndexPsi(start_index=index1,end_index=index1+9,k1=k1,step=stp,linear_compression_ratio=lcr)
    'get file named index1 index2 kl'

    #loadDataToMysql
    R"""
    Note: the format of table_name='depin_from_honeycomb'
        simu_index | HarmonicK | LinearCompressionRatio | Psi3Global | Psi6Global
    """
    import opertateOnMysql as osql
    osql.loadDataToMysql(path_to_file_name=filename_kl,table_name='depin_from_honeycomb')#"/home/tplab/Downloads/193-205kl"




#da.saveIndexPsi()
#da.plotHarmonicKAndPsi()
#sa_h.workflow()
#data = osql.getDataFromMysql('/home/tplab/Downloads/testopsql',table_name='honeycomb')
#print(data)
'''
mysql_conn = pymysql.connect(host='localhost',user='root',password='4KjnuS99n2tZO#ghs',db='hoomd_data')
sql = "insert into honeycomb values (90,800,0.77,6.39212489e-01,1.9444444e-01);"
try:
    with mysql_conn.cursor() as cursor:
        cursor.execute(sql)
    mysql_conn.commit()
except Exception as e:
    mysql_conn.rollback()
'''

'''
df=pd.DataFrame(np.random.randn(10,4))
print(df)#.tail()
df.to_csv('/home/tplab/Downloads/test_csv',)
df.to_sql('/home/tplab/Downloads/test_sql')
'''
#da.rearrange_data()
#da.phase_diagram_line(select=False)

#sa.workflow()
#da.save_index_psi





