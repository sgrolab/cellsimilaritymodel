# Saturated Production sweep PprodA
import pickle
import sys
import numpy as np 
from utils import motiffunc as mf
from utils.config import PROJECT_DIR

nCells = 1000
rng = np.random.default_rng(seed=1000)

PprodAs = np.logspace(-3,2,6)
PprodAindex = int(sys.argv[1])
PprodA = PprodAs[PprodAindex]
Tcc = 1000
kcatA = 10**-1

motherAs = np.zeros(nCells)
motherBs = np.zeros_like(motherAs)
drnd = np.zeros([nCells,6])
dsis = np.zeros_like(drnd)

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('prodsat',[PprodA,kcatA])
motherCell.equilibrate(20)
motherCell.run(nCells)

divStates = motherCell.motherStates[2::]
motherAs = divStates[0]
motherBs = divStates[1]

for k in range(nCells):
    cell1 = rng.binomial(divStates[:,k].astype('int'),0.5)
    cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    dsis[k] = divStates[:,k] - 2*cell1
    drnd[k] = cell1 - cell2

with open(PROJECT_DIR / f'satprod/n1000/satprod_prodAsweep_{PprodAindex}.pickle','wb') as f:
    pickle.dump([PprodA,Tcc,kcatA,motherAs,motherBs,drnd,dsis],f,pickle.HIGHEST_PROTOCOL)
