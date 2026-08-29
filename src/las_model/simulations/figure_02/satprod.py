# Binding motif initial values 
import pickle
import numpy as np 
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR

Tcc = 1000
nCells = 1000
rng = np.random.default_rng(seed=1000)

prodAs = np.logspace(-3,2,6)

motherCells = []
divStates = np.zeros([len(prodAs),5,nCells])

for i in range(len(prodAs)):

    motherCell = mf.Cell(Tcc,0)
    motherCell.parameterize('prodsat',[prodAs[i],0.01])
    motherCell.run(nCells)
    
    motherCells.append(motherCell)
    
    divStates[i] = motherCell.getMotherStates()
    
dsis = np.zeros([len(prodAs),nCells,5])
drnd = np.zeros([len(prodAs),nCells,5])

for j in range(len(prodAs)):
    for i in range(nCells):
        cell1 = rng.binomial(divStates[j,:,i].astype('int'),0.5)
        cell2 = rng.binomial(divStates[j,:,rng.integers(0,nCells)].astype('int'),0.5)
        
        dsis[j,i] = divStates[j,:,i] - 2*cell1
        drnd[j,i] = cell1 - cell2

with open(PROJECT_DIR / 'analyticalData/motifs_prodsat.pickle','wb') as f:
    pickle.dump([prodAs,motherCells,dsis,drnd],f,pickle.HIGHEST_PROTOCOL)
