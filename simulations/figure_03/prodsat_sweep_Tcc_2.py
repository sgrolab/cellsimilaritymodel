# Satured Production: Sweep Tcc (high PprodA, kcatA)
import sys 
import pickle
import numpy as np 
from utils import motiffunc as mf
from utils.config import PROJECT_DIR

Tccs = [500,1000,2000,5000,10000]
Tccindex = int(sys.argv[1])
Tcc = Tccs[Tccindex]

nCells = 1000
rng = np.random.default_rng(seed=1000)

PprodA = 10**-1
kcatA = 10**-1

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('prodsat',[PprodA,kcatA])
motherCell.equilibrate(20)
motherCell.run(nCells)

divStates = motherCell.getMotherStates()

dsis = np.zeros([nCells,6])
drnd = np.zeros([nCells,6])

for i in range(nCells):
    cell1 = rng.binomial(divStates[:,i].astype('int'),0.5)
    cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    dsis[i] = divStates[:,i] - 2*cell1
    drnd[i] = cell1 - cell2

normvarA = 1-np.var(dsis[:,0],axis=0)/np.var(drnd[:,0],axis=0)
normvarB = 1-np.var(dsis[:,1],axis=0)/np.var(drnd[:,1],axis=0)

with open(PROJECT_DIR / 'prodsat_sweep/prodsat_Tccsweep/prodsat_Tccsweep_Tcc_%.2i.pickle' % (Tccindex),'wb') as f:
    pickle.dump([Tcc,normvarA,normvarB],f,pickle.HIGHEST_PROTOCOL)
