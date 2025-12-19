# Saturated Production: Sweep Tcc 
import pickle
import numpy as np 
from utils import motiffunc as mf
from utils.config import PROJECT_DIR

nCells = 1000
rng = np.random.default_rng(seed=1000)

Tccs = np.logspace(2,4,5)
kcatA = 10**-2
PprodA = 10**-2

motherCells = []
divStates = np.zeros([len(Tccs),5,nCells])

for i in range(len(Tccs)):

    motherCell = mf.Cell(Tccs[i],0)
    motherCell.parameterize('prodsat',[PprodA,kcatA])
    motherCell.run(nCells)
    
    motherCells.append(motherCell)
    
    divStates[i] = motherCell.getMotherStates()
    
dsis = np.zeros([len(Tccs),nCells,5])
drnd = np.zeros([len(Tccs),nCells,5])

for j in range(len(Tccs)):
    for i in range(nCells):
        cell1 = rng.binomial(divStates[j,:,i].astype('int'),0.5)
        cell2 = rng.binomial(divStates[j,:,rng.integers(0,nCells)].astype('int'),0.5)
        
        dsis[j,i] = divStates[j,:,i] - 2*cell1
        drnd[j,i] = cell1 - cell2

with open(PROJECT_DIR / 'analyticalData/motifs_prodsat_Tccsweep.pickle','wb') as f:
    pickle.dump([Tccs,motherCells,dsis,drnd],f,pickle.HIGHEST_PROTOCOL)
