from pymol import cmd
import os

class PDBparase():
    def __init__(self,pdbid,pdbfilepath,out_path=os.path.abspath('.')):
        self.pdbid = pdbid
        self.path = pdbfilepath
        self.outpath = open(os.path.join(out_path,'Bindsite.txt'),'w')
        self.aadic = {
         'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E',
         'PHE':'F', 'GLY':'G', 'HIS':'H', 'LYS':'K',
         'ILE':'I', 'LEU':'L', 'MET':'M', 'ASN':'N',
         'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S',
         'THR':'T', 'VAL':'V', 'TYR':'Y', 'TRP':'W'}#定义一个字典，将三字母氨基酸转换为单字母


    def extractChain(self):
        '''
        extract pdb file Chain list
        :return:
        '''
        cmd.delete('all')
        cmd.load(self.path)
        cmd.remove('solvent')  ##移除溶剂分子和水分子
        # cmd.select('proteinseq','  (byres polymer & name CA)')
        cmd.select('all_','all')
        prot_atom = cmd.get_model('all_').atom
        chain_prot = [ atom.chain for atom in prot_atom]
        chain_lig = sorted(list(set(chain_prot)))
        return chain_lig


    def extractseq(self):
        '''
        提取pdb文件中每条链的蛋白质序列
        :return: 
        '''''
        chain_ls = self.extractChain()
        prot_chain = []
        for c in chain_ls:
            cmd.select('singleseq', '  (byres polymer & name CA) in chain %s'%c)
            singleseq_atom = cmd.get_model('singleseq').atom
            seq = ''
            aanum = []
            for a in singleseq_atom:
                if a.resn in self.aadic:
                    if a.resi not in aanum:
                        seq += self.aadic[a.resn]
                        aanum.append(a.resi)
            if aanum:
                print('>%s\t%s\n%s'%(self.pdbid,c,seq))
                prot_chain.append(c)
            cmd.delete('singleseq')
        return prot_chain
    ##提取相互作用蛋白质链
    def extractinterchain(self):
        chain_ls = self.extractChain()
        prot_chain = []
        for c in chain_ls:
            cmd.select('singleseq', '  (byres polymer & name CA) in chain %s' % c)
            singleseq_atom = cmd.get_model('singleseq').atom
            seq = ''
            aanum = []
            for a in singleseq_atom:
                if a.resn in self.aadic:
                    if a.resi not in aanum:
                        seq += self.aadic[a.resn]
                        aanum.append(a.resi)
            if aanum:
                prot_chain.append(c)
            cmd.delete('singleseq')
        return prot_chain
    def extractligand(self):
        '''
        提取pdb文件中每条链上的配体
        :return:
        '''
        cmd.remove('solvent')
        chain_ls = self.extractChain()
        ligand_ls = []
        for c in chain_ls:
            cmd.select('dna', ' resn DA+DT+DC+DG in chain %s'%c)
            dna_atom = cmd.get_model('dna').atom
            if dna_atom:
                ligand_ls.append((c, 'DNA'))
            cmd.select('rna', ' resn G+C+U+A in chain %s' % c)
            rna_atom = cmd.get_model('rna').atom
            rna_atom = cmd.get_model('rna').atom
            if rna_atom:
                ligand_ls.append((c, 'RNA'))
            cmd.select('singleseq', '  (byres polymer & name CA) in chain %s' % c)
            cmd.select('ligand', ' (!singleseq & !dna & !rna) in chain %s' % c)
            ligand_atom = cmd.get_model('ligand').atom
            for a in ligand_atom:
                if a.chain == c:
                    ligand_ls.append((c, a.resn))
        return sorted(set(ligand_ls))
    def calbind(self,cutoff=3.5):
        inter_chain = self.extractinterchain()
        ligand_ls = self.extractligand()
        for l in ligand_ls:
            chain = l[0]
            ligand = l[1]
            if l[1] == 'DNA':
                ligand = 'DA+DT+DC+DG'
            if l[1] == 'RNA':
                ligand = 'G+C+U+A'
            cmd.select('lig',' resn %s in chain %s'%(ligand,chain))
            for c in inter_chain:
                cmd.select('interseq', '  (byres polymer & name CA) in chain %s'%c)
                cmd.select('intersite', 'byres interseq around %s ' %cutoff)
                bsnum_ls = []
                site_ls = []
                for a in cmd.get_model('intersite').atom:
                    if a.resn in self.aadic:
                        aazong = self.aadic[a.resn] + a.resi
                        if a.resi not in bsnum_ls:
                            bsnum_ls.append(a.resi)
                        if aazong not in site_ls:
                            site_ls.append(aazong)
                seq_atom = cmd.get_model('interseq').atom
                seq = ''
                aanum = []
                for a in seq_atom:
                    if a.resn in self.aadic:
                        if a.resi not in aanum:
                            if a.resi in bsnum_ls:
                                seq += self.aadic[a.resn].lower()
                            else:
                                seq += self.aadic[a.resn]
                            aanum.append(a.resi)
                if ligand == 'DA+DT+DC+DG':
                    ligand = 'DNA'
                if ligand == 'G+C+U+A':
                    ligand = 'RNA'
                if site_ls:
                    print('%s\t%s\t%s\t%s\t%s\t%s' % (self.pdbid, chain, ligand, c, ' '.join(site_ls), seq),file=self.outpath)
                    print('%s\t%s\t%s\t%s\t%s\t%s'%(self.pdbid,chain,ligand,c,' '.join(site_ls),seq))
        cmd.delete('all')
    def calpip(self,cutoff=3.5):
        inter_chain = self.extractinterchain()
        for c in inter_chain:
            cmd.select('rcpseq', '  (byres polymer & name CA) in chain %s' % c)##定义受体链
            for d in inter_chain:
                if d != c:
                    cmd.select('interpip', '  (byres polymer & name CA) in chain %s' % d)##定义相互作用链
                    cmd.select('protsite',' byres interpip around %s in rcpseq'%cutoff)
                    bsnum_ls = []
                    site_ls = []
                    for a in cmd.get_model('protsite').atom:
                        if a.resn in self.aadic:
                            aazong = self.aadic[a.resn] + a.resi
                            if a.resi not in bsnum_ls:
                                bsnum_ls.append(a.resi)
                            if aazong not in site_ls:
                                site_ls.append(aazong)
                    seq_atom = cmd.get_model('protsite').atom
                    seq = ''
                    aanum = []
                    for a in seq_atom:
                        if a.resn in self.aadic:
                            if a.resi not in aanum:
                                if a.resi in bsnum_ls:
                                    seq += self.aadic[a.resn].lower()
                                else:
                                    seq += self.aadic[a.resn]
                                aanum.append(a.resi)
                    if site_ls:
                        print('%s\t%s\t%s\t%s\t%s\t%s' % (self.pdbid, d,'PIP', c, ' '.join(site_ls), seq),
                              file=self.outpath)
                        print('%s\t%s\t%s\t%s\t%s\t%s' % (self.pdbid, d,'PIP', c, ' '.join(site_ls), seq))
pdb = PDBparase('4qoz','5mxd.pdb')

# print(pdb.extractligand())

pdb.calbind()
