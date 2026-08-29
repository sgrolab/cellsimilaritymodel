# Binding motif initial values 
import pickle
import sys 
import numpy as np 
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR

nCells = 1000
rng = np.random.default_rng(seed=1000)

Tcc = 1000
PprodA = 10**-1
kcats = np.logspace(-3,3,13)
kcatAindex = int(sys.argv[1])
kcatA = kcats[kcatAindex]

motherCells = []
divStates = np.zeros([len(kcats),6,nCells])

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('prodsat',[PprodA,kcatA])
motherCell.equilibrate(20)
motherCell.run(nCells)

motherCells.append(motherCell)

divStates = motherCell.motherStates[2::]

dsis = np.zeros([nCells,6])
drnd = np.zeros([nCells,6])

for i in range(nCells):
    cell1 = rng.binomial(divStates[:,i].astype('int'),0.5)
    cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    dsis[i] = divStates[:,i] - 2*cell1
    drnd[i] = cell1 - cell2

normvarA = 1-np.var(dsis[:,0],axis=0)/np.var(drnd[:,0],axis=0)
normvarB = 1-np.var(dsis[:,1],axis=0)/np.var(drnd[:,1],axis=0)

with open(PROJECT_DIR / 'prodsat_sweep/prodsat_kcatsweep/prodsat_kcatsweep_kcatA_%.2i.pickle' % (kcatAindex),'wb') as f:
    pickle.dump([kcatA,normvarA,normvarB],f,pickle.HIGHEST_PROTOCOL)
