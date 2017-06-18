def makeCache(inputFile, outputFile):
    import pickle
    pool = []
    tf = open(inputFile, mode='r')
    cf = open(outputFile, mode='wb')
    ini = 1
    for line in tf.readlines():
        if ini:
            ini = 0
        else:
            items = line.split('\t')
            pool.append(items[0]+','+items[1])
    pickle.dump(list(set(pool)), cf)
    tf.close()
    cf.close()
    
if __name__ == "__main__":
    makeCache(raw_input('Input file: '), raw_input('Output file: '))