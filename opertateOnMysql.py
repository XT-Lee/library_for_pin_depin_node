import numpy as np
from numpy.lib.function_base import select
import pymysql
def introduction():
    R"""
    Function list:
        connectToMysql: it helps program connect to Mysql database.
        loadDataToMysql: it is used to load data to Mysql database.
        getDataFromMysql: it is used to extract data from Mysql database through 'data=getDataFromMysql()'.
    Reference：
        https://zhuanlan.zhihu.com/p/51553625
    """
    pass
def loadDataToMysql(path_to_file_name,table_name):
    R"""
    Examples:
        operateOnMysql.loadDataToMysql(
        path_to_file_name='/home/tplab/Downloads/91-192kl',
        table_name='honeycomb')
    """
    
    mysql_conn = connectToMysql()
    data=np.loadtxt(path_to_file_name)
    shape_data=np.shape(data)

    for row in np.linspace(1,shape_data[0],shape_data[0]):
        data_to_string='('
        #print(data[row.astype(int)-1])
        for column in np.linspace(1,shape_data[1],shape_data[1]):
            data_to_string=data_to_string+str(data[row.astype(int)-1,column.astype(int)-1])
            if column < shape_data[1]:
                data_to_string=data_to_string+','
        
        data_to_string=data_to_string+');'
        #print(data_to_string)
        sql = 'insert into '+table_name+' values'+data_to_string
        try:
            with mysql_conn.cursor() as cursor:
                cursor.execute(sql)
            mysql_conn.commit()
        except Exception as e:
            print(e)
            mysql_conn.rollback()

    mysql_conn.close()

def getDataFromMysql(path_to_save_file=None,table_name='',search_condition='',):
    R"""
    Examples:
        data = operateOnMysql.getDataFromMysql(
                path_to_save_file='/home/tplab/Downloads/91-192kl'
                table_name='honeycomb',
                search_condition='where HarmonicK > 100')
    """
    
    mysql_conn = connectToMysql()

    sql = 'select * from '+table_name+' '+search_condition+';'#where column_name = xx
    try:
        with mysql_conn.cursor() as cursor:
            cursor.execute(sql)
        select_result = cursor.fetchall()
        #np.savetxt(path_to_save_file,select_result)
        return select_result
    except Exception as e:
        print(e)
        
        mysql_conn.rollback()

    mysql_conn.close()

def connectToMysql():
    mysql_conn = pymysql.connect(host='localhost',user='root',password='4KjnuS99n2tZO#ghs',db='hoomd_data')
    return mysql_conn