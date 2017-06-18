#!/usr/bin/env python
from database import mysqlConnection


def cache_single(input_file):
    pool = []
    ini = 1
    tf = open(input_file, mode='r')
    for line in tf.readlines():
        if ini:
            ini = 0
        else:
            items = line.split('\n')[0].split('\t')
            pool.append(items[0]+','+items[1])
    pool = list(set(pool))
    tf.close()
    return pool


def make_cache():
    import pickle
    big_pool = dict()
    big_pool['mir'] = cache_single(raw_input('For miRNA: '))
    big_pool['utr'] = cache_single(raw_input('For 3\' UTR: '))
    big_pool['mis'] = cache_single(raw_input('For missense: '))
    big_pool['alt'] = cache_single(raw_input('For alternative splicing: '))
    output_file = raw_input('Save to: ')
    cf = open(output_file, mode='wb')
    pickle.dump(big_pool, cf)
    cf.close()


def load_cache(cache_file='/mnt/rep/dispatch/ref/cache.dat'):
    import pickle
    cf = open(cache_file, mode='rb')
    return pickle.load(cf)


def copy_mir(chromosome, pos, job_id):
    try:
        cnx, cursor = mysqlConnection()
        overview_sql = """SELECT `mirna`, `accession`, `edit_pos_raw`, `edit_pos_chr`, `event`, `role`, `sequence`, `old_t`, `new_t`, `com_t`, `sig` FROM `smir` WHERE `edit_pos_chr`='%s' AND `edit_pos_raw`=%s;"""%(chromosome, pos)
        cursor.execute(overview_sql)
        res = cursor.fetchone()
        keyword = ''
        if res:
            res = list(res)
            keyword = res[10]
            res.append(job_id)
            insert_ov_sql = """INSERT INTO `mirediting` (`mirna`, `accession`, `edit_pos_raw`, `edit_pos_chr`, `event`, `role`, `sequence`, `old_t`, `new_t`, `com_t`, `sig`, `job`) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');"""%(res[0], res[1], res[2], res[3], res[4], res[5], res[6], res[7], res[8], res[9], res[10], res[11])
            cursor.execute(insert_ov_sql)
        else:
            print 'ggc'

        seed_tar_sql = """SELECT `mirna`, `sig`, `tag`, `way`, `gene_symbol`, `transcript_id`, `dg_duplex`, `dg_binding`, `dg_duplex_seed`, `dg_binding_seed`, `utr_start`, `utr_end`, `utr3` FROM `sms` WHERE `sig`='%s';"""%(keyword)
        cursor.execute(seed_tar_sql)
        res = cursor.fetchall()
        for raw in res:
            raw = list(raw)
            raw.append(job_id)
            insert_st_sql = """INSERT INTO `mirmap_score` (`mirna`, `sig`, `tag`, `way`, `gene_symbol`, `transcript_id`, `dg_duplex`, `dg_binding`, `dg_duplex_seed`, `dg_binding_seed`, `utr_start`, `utr_end`, `utr3`, `job`) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');"""%(raw[0], raw[1], raw[2], raw[3], raw[4], raw[5], raw[6], raw[7], raw[8], raw[9], raw[10], raw[11], raw[12], raw[13])
            cursor.execute(insert_st_sql)
        mature_tar_sql = """SELECT `mirna`, `sig`, `tag`, `way`, `gene_symbol`, `score`, `energy`, `mi_start`, `mi_end`, `utr_start`, `utr_end`, `match_len`, `identity`, `similarity`, `mir_seq`, `lines`, `utr_seq` FROM `smiranda` WHERE `sig`='%s';"""%(keyword)
        cursor.execute(mature_tar_sql)
        res = cursor.fetchall()
        for raw in res:
            raw = list(raw)
            raw.append(job_id)
            insert_mt_sql = """INSERT INTO `miranda_score` (`mirna`, `sig`, `tag`, `way`, `gene_symbol`, `score`, `energy`, `mi_start`, `mi_end`, `utr_start`, `utr_end`, `match_len`, `identity`, `similarity`, `mir_seq`, `lines`, `utr_seq`, `job`) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', %s);"""%(raw[0], raw[1], raw[2], raw[3], raw[4], raw[5], raw[6], raw[7], raw[8], raw[9], raw[10], raw[11], raw[12], raw[13], raw[14], raw[15], raw[16], raw[17])
            cursor.execute(insert_mt_sql)
        cnx.commit()
        cnx.close()
        return 0
    except Exception, e:
        print e
        return 1


def copy_utr(chromosome, pos, job_id):
    try:
        cnx, cursor = mysqlConnection()
        overview_sql = """SELECT `gene`, `edit_pos_chr`, `edit_pos_raw`, `event`, `old_t`, `new_t`, `com_t`, `sig` FROM `sutr3` WHERE `edit_pos_chr`='%s' AND `edit_pos_raw`='%s';"""%(chromosome, pos)
        cursor.execute(overview_sql)
        res = cursor.fetchone()
        keyword = ''
        if res:
            keyword = res[7]
            insert_ov_sql = """INSERT INTO `utr3editing` (`gene`, `edit_pos_chr`, `edit_pos_raw`, `event`, `old_t`, `new_t`, `com_t`, `sig`, `job`) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', %s);"""%(res[0], res[1], res[2], res[3], res[4], res[5], res[6], res[7], job_id)
            cursor.execute(insert_ov_sql)
        else:
            print 'gga'
            
        tar_sql = """SELECT `gene`, `sig`, `tag`, `way`, `mir`, `dg_duplex`, `dg_binding`, `dg_duplex_seed`, `dg_binding_seed`, `utr_start`, `utr_end`, `utr3` FROM `sum` WHERE `sig`='%s';"""%(keyword)
        cursor.execute(tar_sql)
        res = cursor.fetchall()
        for raw in res:
            raw = list(raw)
            raw.append(job_id)        
            targets_sql = """INSERT INTO `utr_map` (`gene`, `sig`, `job`, `tag`, `way`, `mir`, `dg_duplex`, `dg_binding`, `dg_duplex_seed`, `dg_binding_seed`, `utr_start`, `utr_end`, `utr3`) VALUES ('%s', '%s', %s, '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');"""%(raw[0], raw[1], job_id, raw[2], raw[3], raw[4], raw[5], raw[6], raw[7], raw[8], raw[9], raw[10], raw[11])
            cursor.execute(targets_sql)
        cnx.commit()
        cnx.close()
        return 0
    except Exception, e:
        print e
        return 1


def copy_mis(chromosome, pos, job_id):
    try:
        cnx, cursor = mysqlConnection()
        overview_sql = """SELECT `chromosome`, `position`, `gene`, `transcript`, `relpos`, `fr`, `to` FROM `smis` WHERE `chromosome`='%s' AND `position`='%s';""" % (chromosome, pos)
        cursor.execute(overview_sql)
        res = cursor.fetchone()
        if res:
            insert_mis_sql = """INSERT INTO `missense` (`job`, `chromosome`, `position`, `gene`, `transcript`, `relpos`, `fr`, `to`) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');""" % (job_id, res[0], res[1], res[2], res[3], res[4], res[5], res[6])
            cursor.execute(insert_mis_sql)
        cnx.commit()
        cnx.close()
        return 0
    except Exception, e:
        print e
        return 1


def copy_alter_splicing(chromosome, pos, job_id):
    try:
        cnx, cursor = mysqlConnection()
        ask_at_sql = """SELECT `gene`, `pos`, `type`, `raw_score`, `new_score`, `variation`, `order`, `transcript` FROM `sse` WHERE `chromosome`='%s' AND `position`=%s;""" % (chromosome, pos)
        cursor.execute(ask_at_sql)
        res = cursor.fetchone()
        if res:
            insert_at_sql = """INSERT INTO `splicing_event` (`gene`, `pos`, `type`, `raw_score`, `new_score`, `variation`, `order`, `transcript`, `job`) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', %s);""" % (res[0], res[1], res[2], res[3], res[4], res[5], res[6], res[7], job_id)
            cursor.execute(insert_at_sql)
        cnx.commit()
        cnx.close()
        return 0
    except Exception, e:
        print e
        return 1


def copy_overview(chromosome, pos, job_id):
    try:
        cnx, cursor = mysqlConnection()
        overview_sql = """SELECT `Frequency`, `Position`, `RefSeq_feat`, `RefSeq_gid`, `Region`, `RepMask_gid`, `Strand` FROM `searchable` WHERE `Region`='%s' AND `Position`=%s;"""%(chromosome, pos)
        cursor.execute(overview_sql)
        res = cursor.fetchone()
        if res:
            freq = res[0] if res[0] != None else 0
            insert_ov_sql = """INSERT INTO `candidates` (`Frequency`, `Job`, `Position`, `RefSeq_feat`, `RefSeq_gid`, `Region`, `RepMask_gid`, `Strand`) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');"""%(freq, job_id, res[1], res[2], res[3], res[4], res[5], res[6])
            cursor.execute(insert_ov_sql)
            cnx.commit()
        else:
            print 'ggb'
        cnx.close()
        return 0
    except Exception, e:
        print e
        return 1
    

def main(job_file, job_id, header=True):
    cached = load_cache()
    jf = open(job_file, mode='r')
    cf = open(job_file+'.nd', mode='w')
    ef_mir = 0    
    ef_mis = 0
    ef_utr = 0
    ef_alt = 0
    if header:
        candidates = jf.readlines()
        cf.write(candidates[0])
        candidates = candidates[1:]
    else:
        candidates = jf.readlines()

    for line in candidates:
        items = line.split('\t')
        keyword = items[0]+','+items[1]
        #if keyword in cached['mir']:
        #    ret1 = copy_overview(items[0], items[1], job_id)
        #    ret2 = copy_mir(items[0], items[1], job_id)
        #    if ret1 or ret2:
        #        cf.write(line)
        if keyword in cached['utr']:
            ret1 = copy_overview(items[0], items[1], job_id)
            ret2 = copy_utr(items[0], items[1], job_id)
            ef_utr += 1
            if ret1 or ret2:
                cf.write(line)
        elif keyword in cached['mis']:
            ret1 = copy_overview(items[0], items[1], job_id)
            ret2 = copy_mis(items[0], items[1], job_id)
            ef_mis += 1
            if ret1 or ret2:
                cf.write(line)
            '''
        elif keyword in cached['alt']:
            ret1 = copy_overview(items[0], items[1], job_id)
            ret2 = copy_alter_splicing(items[0], items[1], job_id)
            ef_alt += 1
            if ret1 or ret2:
                cf.write(line)
            '''
        else:
            cf.write(line)
    jf.close()
    cf.close()
    return ef_mir, ef_utr, ef_mis, ef_alt

if __name__ == '__main__':
    make_cache()
