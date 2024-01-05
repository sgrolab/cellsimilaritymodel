# Binding motif initial values 

import os, pickle, numpy as np 
os.chdir('/groups/sgro/sgrolab/mark/comp_proj/motifs')
import motiffunc as mf
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/motifs')


Tcc = 1000
nCells = 1000
rng = np.random.default_rng(seed=1000)

varTccs = np.logspace(0,3,11)
divStates = np.zeros([len(varTccs),5,nCells])

for i in range(len(varTccs)):

    motherCell = mf.Cell(Tcc,varTccs[i])
    motherCell.parameterize('single',[.1])
    motherCell.run(nCells)
    
    divStates[i] = motherCell.getMotherStates()
    
dsis = np.zeros([len(varTccs),nCells,5])
drnd = np.zeros([len(varTccs),nCells,5])

for j in range(len(varTccs)):
    for i in range(nCells):
        cell1 = rng.binomial(divStates[j,:,i].astype('int'),0.5)
        cell2 = rng.binomial(divStates[j,:,rng.integers(0,nCells)].astype('int'),0.5)
        
        dsis[j,i] = divStates[j,:,i] - 2*cell1
        drnd[j,i] = cell1 - cell2

os.chdir('/groups/sgro/sgrolab/mark/comp_proj/varTcc')
with open('varTcc_diffs.pickle','wb') as f:
    pickle.dump([dsis,drnd],f,pickle.HIGHEST_PROTOCOL)
