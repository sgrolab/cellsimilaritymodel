# Phosphorylation Monocycle: Sweep PprodA and PprodB
import sys 
import pickle
import numpy as np 
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR

nCells = 2000
rng = np.random.default_rng(seed=1000)

Tcc = 1000

PprodAs = np.logspace(-3,2,31)
PprodAindex = int(sys.argv[1])
PprodA = PprodAs[PprodAindex]

PprodBs = np.logspace(-3,2,31)
PprodBindex = int(sys.argv[2])
PprodB = PprodBs[PprodBindex]

PprodC = 10**-1
ka = 10**-1
kc = 10**-1

params = [PprodA,PprodB,PprodC,ka,kc]

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('phos_cycle',params)
motherCell.run(nCells)

Aeq = np.mean(motherCell.A/motherCell.V)
Beq = np.mean(motherCell.B/motherCell.V)
Ceq = np.mean(motherCell.C/motherCell.V)
Deq = np.mean(motherCell.D/motherCell.V)
Eqs = [Aeq,Beq,Ceq,Deq]

varA = np.var(motherCell.A/motherCell.V)
varB = np.var(motherCell.B/motherCell.V)
varC = np.var(motherCell.C/motherCell.V)
varD = np.var(motherCell.D/motherCell.V)
var = [varA,varB,varC,varD]

divStates = motherCell.getMotherStates()

dsis = np.zeros([nCells,6])
drnd = np.zeros([nCells,6])

for k in range(nCells):
    cell1 = rng.binomial(divStates[:,k].astype('int'),0.5)
    cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    dsis[k] = divStates[:,k] - 2*cell1
    drnd[k] = cell1 - cell2

normvarA = 1-np.var(dsis[:,0],axis=0)/np.var(drnd[:,0],axis=0)
normvarB = 1-np.var(dsis[:,1],axis=0)/np.var(drnd[:,1],axis=0)
normvarC = 1-np.var(dsis[:,2],axis=0)/np.var(drnd[:,2],axis=0)
normvarD = 1-np.var(dsis[:,3],axis=0)/np.var(drnd[:,3],axis=0)
normvars = [normvarA,normvarB,normvarC,normvarD]

with open(PROJECT_DIR / 'phos_cycle/phos_cycle2/phoscycle2_PprodA_%.2i_PprodB_%.2i.pickle' % (PprodAindex,PprodBindex),'wb') as f:
    pickle.dump([params,Eqs,var,normvars],f,pickle.HIGHEST_PROTOCOL)


