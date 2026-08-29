# Saturated Production: 2D Sweep PprodA, kcatA

import os, sys, pickle, numpy as np 
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR

nCells = 1000
rng = np.random.default_rng(seed=1000)

Tccs = np.logspace(2,4,5)
prodAs = np.logspace(-2,0,3)
kcatAs = np.logspace(-4,0,9)

prodAindex = int(sys.argv[1])
prodA = prodAs[prodAindex]
kcatAindex = int(sys.argv[2])
kcatA = kcatAs[kcatAindex]

divStates = np.zeros([len(Tccs),5,nCells])
dsis = np.zeros([len(Tccs),nCells,5])
drnd = np.zeros([len(Tccs),nCells,5])

for i in range(len(Tccs)):

    motherCell = mf.Cell(Tccs[i],0)
    motherCell.parameterize('prodsat',[prodA,kcatA])
    motherCell.run(nCells)
    
    divStates[i] = motherCell.getMotherStates()
    
    for j in range(nCells):
        cell1 = rng.binomial(divStates[i,:,j].astype('int'),0.5)
        cell2 = rng.binomial(divStates[i,:,rng.integers(0,nCells)].astype('int'),0.5)
        
        dsis[i,j] = divStates[i,:,j] - 2*cell1
        drnd[i,j] = cell1 - cell2


with open(PROJECT_DIR / 'prodsat_sweep/sweep1/motifs_prodsat_2dsweep_prodA_' + str(prodAindex) + '_kcatA_' + str(kcatAindex) + '.pickle','wb') as f:
    pickle.dump([prodA,kcatA,Tccs,divStates,dsis,drnd],f,pickle.HIGHEST_PROTOCOL)

