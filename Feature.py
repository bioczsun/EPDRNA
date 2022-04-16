##该文件用来提取特征参数


import os
import re
import math
import numpy as np
import pandas as pd
import json



class ExtractFea():
    def __init__(self,nuc_type='DNA',windows = 7):
        Blastpath = r'E:\date\epdrna\blast\ncbi-blast-2.12.0+\bin'##Blast 文件bin的路径
        self.dbpath = r'E:\date\epdrna\db\swissprot' ##进行迭代数据库的路径和数据库的名字
        self.pssm_outpath = r'.\PSSM'
        self.Blast = os.path.join(Blastpath,'psiblast')
        self.windows = windows
        self.pssmfea = open(r'.\Featuredir\pssmfea.txt', 'w')
        self.pssmlaebl = open(r'Featuredir/label.txt', 'w')
        self.hotpath = open(r'.\Featuredir\one-hot.txt', 'w')
        if nuc_type == 'DNA':
            self.aaindexpath = open('aaindex.json').read()
            self.pcppath = open(r'.\Featuredir\pcp10.txt','w')





    ##用来生成PSSM
    def command_pssm(self,fastaseq,sequence):
        pssmpath = os.path.join(self.pssm_outpath,fastaseq.split('.')[-2]+'.pssm')
        if os.path.exists(pssmpath) == False:
            os.system('%s -query %s -db %s -num_iterations 3 -num_threads 5 -evalue 0.001  -out_ascii_pssm %s'%(self.Blast,
                                                                                                            fastaseq,
                                                                                                            self.dbpath,
                                                                                                            pssmpath))
        self.deal_pssm(pssmpath,sequence)

    def minmax(self,ls_a):
        mix_ls = []
        for n in ls_a:
            mix_ls.append(1 / (1 + math.exp(-float(n))))
        return mix_ls
    ###计算每个氨基酸残基的PSSM特征
    def deal_pssm(self,pssmpath,sequence):
        windows = self.windows
        data = pd.read_table(pssmpath,delimiter='\t', skiprows=2)
        label = []
        pssm_fea = []
        for n in range(0,len(sequence)):
            all_ls = []#用来存放a和a-i和a+i的pssm数据
            sy_ls = []#用来存放a-i和a+i的pssm数据
            if sequence[n].isupper():
                label.append(1)
            if sequence[n].islower():
                label.append(0)
            if n - int(windows / 2) <= 0:
                all_ls.append(self.minmax(data.values[n][0].split()[2:22]))
                for i in range(int((windows / 2) - n)):
                    sy_ls.append(np.zeros(20, dtype=float))
                for i in range(0, n):
                    if data.values[i][0].split()[2:22]:
                        sy_ls.append(self.minmax(data.values[i][0].split()[2:22]))
                for j in range(n + 1, n + int(windows / 2) + 1):
                    sy_ls.append(self.minmax(data.values[j][0].split()[2:22]))
            if 0 < n - (windows / 2) and len(sequence) - n > int(windows / 2):
                all_ls.append(self.minmax(data.values[n][0].split()[2:22]))
                for i in range(n - int(windows / 2), n):
                    sy_ls.append(self.minmax(data.values[i][0].split()[2:22]))
                for i in range(n + 1, n + int(windows / 2) + 1):
                    sy_ls.append(self.minmax(data.values[i][0].split()[2:22]))
            if n + int(windows / 2) >= len(sequence):
                all_ls.append(self.minmax(data.values[n][0].split()[2:22]))
                for i in range(n - int(windows / 2) - 1, n):
                    sy_ls.append(self.minmax(data.values[i][0].split()[2:22]))
                for i in range(n + 1, len(sequence) - 1):
                    if data.values[i][0].split()[2:22]:
                        sy_ls.append(self.minmax(data.values[i][0].split()[2:22]))
                for i in range(int(n + int(windows / 2) - len(sequence) + 1)):
                    sy_ls.append(np.zeros(20, dtype=float))
            all_ls.append(np.array(sy_ls, dtype=float).mean(axis=0))
            a = np.array(all_ls, dtype=float)
            d = a.flatten()
            c = d.reshape(1, 40)
            pssm_fea.append(c)

        for v in np.array(pssm_fea):
            np.savetxt(self.pssmfea, v, delimiter='\t', fmt="%.2f")
        for l in label:
            print(l,file=self.pssmlaebl)

    ###计算每条序列理化性质
    def cal_pcp(self,sequence,typenuc='DNA'):
        aaindex = json.loads(self.aaindexpath)
        win_ls = [3, 5, 7, 9, 11, 13]
        titls = ['TSAJ990101', 'CHOP780204', 'OOBM770101', 'AURR980107', 'GOLD730102',
                 'WOLR790101', 'NAGK730101', 'NAKH920102', 'TANS770110', 'KLEP840101']
        if typenuc == 'DNA':
            for num in range(0,len(sequence)):
                all_window_ls = []
                pcp_ls = []
                for t in titls:
                    psy_ls = []
                    psy_ls.append(aaindex[t][sequence[num].upper()])
                    for windows in win_ls:
                        psy_win = []
                        winseq = ''
                        if num - int(windows / 2) <= 0:
                            foward = 'X' * int((windows / 2) - num) + sequence[0:num + 1]
                            back = sequence[num + 1:num + int(windows / 2) + 1]
                            winseq = '%s%s' % (foward, back)
                        if 0 < num - (windows / 2) and len(sequence) - num - 1 >= int(windows / 2):
                            winseq = '%s%s' % (
                                sequence[num - int(windows / 2):num], sequence[num:num + int(windows / 2) + 1])
                        if num + int(windows / 2) >= len(sequence):
                            back = sequence[num:].strip() + 'M' * int(num + int(windows / 2) - len(sequence) + 1)
                            winseq = '%s%s' % (sequence[num - int(windows / 2):num], back)
                        for s in winseq:
                            psy_win.append(aaindex[t][s.upper()])
                        psy_ls.append(round(float(np.mean(psy_win)), 2))
                    all_window_ls.append(psy_ls)
                np.savetxt(self.pcppath,np.array(all_window_ls).flatten().reshape(1,70), delimiter='\t', fmt="%.2f")

    def one_hot(self,sequence):
        windows_pp = 3
        for num in range(0,len(sequence)):
            hot_ls = []
            seq = ''
            if num - int(windows_pp / 2) <= 0:
                foward = 'X' * int((windows_pp / 2) - num) + sequence[0:num + 1]
                back = sequence[num + 1:num + int(windows_pp / 2) + 1]
                seq = '%s%s' % (foward, back)
            if 0 < num - (windows_pp / 2) and len(sequence) - num - 1 >= int(windows_pp / 2):
                seq = '%s%s' % (sequence[num - int(windows_pp / 2):num], sequence[num:num + int(windows_pp / 2) + 1])
            if num + int(windows_pp / 2) >= len(sequence):
                back = sequence[num:].strip() + 'X' * int(num + int(windows_pp / 2) - len(sequence) + 1)
                seq = '%s%s' % (sequence[num - int(windows_pp / 2):num], back)
            for s in seq:
                AA_dict = {}
                for a in 'ACDEFGHIKLMNPQRSTVWY':
                    AA_dict[a] = 0
                if s == 'X':
                    pass
                else:
                    AA_dict[s.upper()] += 1
                for value in AA_dict.values():
                    hot_ls.append(value)
            np.savetxt(self.hotpath,np.array(hot_ls).reshape(1,60),delimiter='\t', fmt="%d")



if __name__ == '__main__':
    print('该程序为特征提取所用')

