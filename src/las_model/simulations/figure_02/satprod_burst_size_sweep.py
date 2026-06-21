# Burst Size Sweep 
import pickle
import numpy as np 
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR

Tcc = 1000
nCells = 1000
rng = np.random.default_rng(seed=1000)

burstSizes = np.linspace(1,20,20)
prodA_std = 10**0
prodAs = prodA_std / burstSizes

Aeqs = np.zeros([len(burstSizes)])
Beqs = np.zeros_like(Aeqs)
normvarAs = np.zeros_like(Aeqs)
normvarBs = np.zeros_like(Aeqs)

for i in range(len(burstSizes)):
    burstSize = burstSizes[i]
    prodA = prodAs[i]
    
    motherCell = mf.Cell(Tcc,0)
    motherCell.parameterize('prodsat_burst',[prodA,0.01,burstSize])
    motherCell.run(nCells)
    
    Aeqs[i] = np.mean(motherCell.A/motherCell.V)
    Beqs[i] = np.mean(motherCell.B/motherCell.V)
    
    divStates = motherCell.getMotherStates()

    dsis = np.zeros([nCells,5])
    drnd = np.zeros([nCells,5])

    for k in range(nCells):
        cell1 = rng.binomial(divStates[:,k].astype('int'),0.5)
        cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
        
        dsis[k] = divStates[:,k] - 2*cell1
        drnd[k] = cell1 - cell2
    
    normvarAs[i] = 1-np.var(dsis[:,0],axis=0)/np.var(drnd[:,0],axis=0)
    normvarBs[i] = 1-np.var(dsis[:,1],axis=0)/np.var(drnd[:,1],axis=0)

with open(PROJECT_DIR / 'burstSize/burstSize_prodA-0.pickle','wb') as f:
    pickle.dump([burstSizes,prodAs,Aeqs,Beqs,normvarAs,normvarBs],f,pickle.HIGHEST_PROTOCOL)