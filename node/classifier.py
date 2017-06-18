from svmutil import *

def train(fileName):
    y, x = svm_read_problem(fileName)
    print y, x    
    m = svm_train(y[:], x[:], '-c 4')
    svm_save_model(fileName+'.model', m)

#def predict():

def predictBySVM(modelFile, x, y):
    m = svm_load_model(modelFile)
    p_label, p_acc, p_val = svm_predict(y, x, m)
    return p_label

def formatInputForLibSVM(row):
    x = []
    #for row in matrix:
    #    items = row.split('\t')
    single = {}
    for k, v in enumerate(row):
        single[k+1] = row[k]
    x.append(single)
    y = [0] * len(x)
    return x, y
