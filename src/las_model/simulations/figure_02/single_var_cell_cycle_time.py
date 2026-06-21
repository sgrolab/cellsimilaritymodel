# Binding motif initial values 
import pickle
import sys
import numpy as np 
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR


nCells = 1000
rng = np.random.default_rng(seed=1000)

Tcc = 1000
varTccs = np.logspace(0,3,11)
PprodA = 10**-1

divStates = np.zeros([len(varTccs),5,nCells])

for i in range(len(varTccs)):

    motherCell = mf.Cell(Tcc,varTccs[i])
    motherCell.parameterize('single',[PprodA])
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

with open(PROJECT_DIR / 'varTcc/varTcc_diffs.pickle','wb') as f:
    pickle.dump([dsis,drnd],f,pickle.HIGHEST_PROTOCOL)
