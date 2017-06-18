from database import mysqlConnection
from repConfig import getConfig

def getmiRnaAcc(mirna):
    #cnx = mysqlConnection(user=getConfig("db", "user"),
    #                              password=getConfig("db", "password"),
    #                              database=getConfig("db", "db_name"),
    #                              host=getConfig("db", "host"),
    #                              port=getConfig("db", "port"))
    #cursor = cnx.cursor(buffered=True)    
    connection, cursor = mysqlConnection()
    asMature1 = "SELECT mature1_acc FROM mirna WHERE mature1_id = '%s'" % mirna
    cursor.execute(asMature1)
    acc = cursor.fetchone()
    if (acc == None):
        asMature2 = "SELECT mature2_acc FROM mirna WHERE mature2_id = '%s'" % mirna
        cursor.execute(asMature2)    
        acc = cursor.fetchone()        
    
    return str(acc[0])

def getmiRnaSymbol(con, cur, mirna):
    #format mirna symbol
    mirna = mirna.replace("-5p", "").replace("-3p", "")
    query = """SELECT `symbol` FROM %s WHERE `mir_id`='%s';""" % (getConfig("dispatch", "mirInfo"), mirna.lower())
    cur.execute(query)
    syn = cur.fetchone()
    if not (syn):
        return 1
    else:
        return syn[0]