import os
import numpy as np

from sklearn.preprocessing import minmax_scale
from Feature import ExtractFea



def dealdata(fasta_path,out_path,nuc_type = 'DNA'):
    fasta_path = fasta_path
    fussion_path = '%s.npy'%out_path
    fasta_path = fasta_path ##文件路径
    fasta_data = open(fasta_path).readlines()

    data_dict = {} ##用来存放fasta数据

    ###处理fasta数据，key：sequence id, values： sequence
    key = ''
    for i in fasta_data:
        i = i.strip()
        if i[0] == '>':
            key = i
            data_dict[key] = ''
        else:
            data_dict[key] = i
    feature = ExtractFea(nuc_type)
    if nuc_type == 'RNA':
        feature = ExtractFea(nuc_type,windows=5)
    for key,values in data_dict.items():

        single_seq = open(r'%s.fasta'%key[1:],'w')
        print('>%s\n%s'%(key,values),file=single_seq)
        single_seq.close()
        feature.command_pssm(r'%s.fasta'%key[1:],values)
        os.remove(r'%s.fasta'%key[1:])
        if nuc_type == 'DNA':
            feature.cal_pcp(values)
        feature.one_hot(values)
    feature.pssmfea.close()
    if nuc_type == 'DNA':
        feature.pcppath.close()
    feature.hotpath.close()
    feature.pssmlaebl.close()
    print('----------------------------------------Feature extraction is completed----------------------------------------')
    pssmdata = np.loadtxt(r'.\Featuredir\pssmfea.txt')
    os.remove(path=r'.\Featuredir\pssmfea.txt')
    ont_data = np.loadtxt(r'.\Featuredir\one-hot.txt')
    os.remove(r'.\Featuredir\one-hot.txt')
    label = np.loadtxt(r'Featuredir/label.txt')
    all_data = np.hstack((pssmdata, ont_data))
    if nuc_type == 'DNA':
        pcp_data = np.loadtxt(r'.\Featuredir\pcp10.txt')
        pcp_data = minmax_scale(pcp_data)
        all_data = np.hstack((pssmdata,pcp_data,ont_data))
    np.save(fussion_path,all_data)


if __name__ == '__main__':
    fasta_path = r'data/test_dna.txt'
    fussion_path = r'.\Featuredir\test_dna.npy'
    dealdata(fasta_path,fussion_path,nuc_type='DNA')