import sqlite3

def conSQLite(db):
    if (db == 'memory'):
        con = sqlite3.connect(':memory:')
    else:
        try:
            con = sqlite3.connect(db)
        except:
            return 0

    cursor = con.cursor()
    return con, cursor