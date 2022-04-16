import joblib
import numpy as np
import dealdata

model_path = r'.\Model\dnamodel.pkl'


fasta_path = r'data/test_dna.txt'

out_apth = open(r'data/results.txt','w')

fasta_path = fasta_path ##文件路径
fasta_data = open(fasta_path).readlines()

data_dict = {} ##用来存放fasta数据
label = np.loadtxt(r'data/testdnalabel.txt')
###deal fasta data，key：sequence id, values： sequence
dealdata.dealdata(fasta_path,r'data/test_dna')

data = np.load(r'data/test_dna.npy')



model = joblib.load(model_path)
key = ''
for i in fasta_data:
    i = i.strip()
    if i[0] == '>':
        key = i
        data_dict[key] = ''
    else:
        data_dict[key] = i
pred = model.predict(data)
proba = model.predict_proba(data)
proba_ls = []
for i,j in zip(pred,proba):
    proba_ls.append(int(j[int(i)]*10))
num = 0
for ids,seq in data_dict.items():
    s = ''
    p = ''
    for i in pred[num:len(seq)+num]:
        s+=str(int(i))
    for j in proba_ls[num:len(seq)+num]:
        p+=str(j)
    print(s,file=out_apth)
    print(p,file=out_apth)
    print(ids,file=out_apth)
    print(seq,file=out_apth)
    num+=len(seq)

