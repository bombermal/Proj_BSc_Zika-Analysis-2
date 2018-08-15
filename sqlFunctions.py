import mysql.connector
from mysql.connector import Error
import readConfig as rc

def insertDataframe(df, dbName):
    temp = []
    for ii in df.iterrows():
        temp.append(tuple(ii[1]))
        
    columns = returnColumnsList(dbName)
    insertData(dbName, columns, temp)


def returnColumnsList(dbName):
    """
        Retorna a lista com os nomes da tabela desejada
    """ 
    temp = []
    query = "SHOW COLUMNS FROM "+dbName
    try:
        conn = mysql.connector.connect(**rc.read_db_config("zyka"))

        cursor = conn.cursor()
        cursor.execute(query)

        row = cursor.fetchone()

        while row is not None:
            temp.append(row[0])
            row = cursor.fetchone()
    except Error as error:
        print(error)

    finally:
        cursor.close()
        conn.close()
        
    return temp

def insertData(dbName, columnsList, valuesList):
    valueRepetition = ["%s" for _ in columnsList]
    
    query = "INSERT INTO "+dbName+"("+','.join(columnsList)+") VALUES("+','.join(valueRepetition)+")"
    
    try:
        conn = mysql.connector.connect(**rc.read_db_config("zyka"))
 
        cursor = conn.cursor()
        cursor.executemany(query, valuesList)
        conn.commit()
        
        print("Dados inseridos em: "+dbName)
    except Error as error:
        print(error)
 
    finally:
        cursor.close()
        conn.close()