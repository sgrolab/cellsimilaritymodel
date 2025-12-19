# Unsaturated production sweep enzyme and substrate production rates
import sys 
import pickle
import numpy as np 
from utils import motiffunc as mf
from utils.config import PROJECT_DIR

nCells = 1000
rng = np.random.default_rng(seed=1000)

Tcc = 1000
PprodAs = np.logspace(-2,1,16)
PprodAindex = int(sys.argv[1])
PprodA = PprodAs[PprodAindex]
PprodBs = np.logspace(-2,4,31)
PprodBindex = int(sys.argv[2])
PprodB = PprodBs[PprodBindex]
kcatA = 10**-1
Km = 10**3

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('produnsat',[PprodB,PprodA,kcatA,Km])
motherCell.equilibrate()
motherCell.run(nCells)

Aeq = np.mean(motherCell.A/motherCell.V)
Beq = np.mean(motherCell.B/motherCell.V)
Ceq = np.mean(motherCell.C/motherCell.V)
means = [Aeq,Beq,Ceq]

Avar = np.var(motherCell.A/motherCell.V)
Bvar = np.var(motherCell.B/motherCell.V)
Cvar = np.var(motherCell.C/motherCell.V)
variances = [Avar,Bvar,Cvar]

divStates = motherCell.getMotherStates()

dsis = np.zeros([nCells,6])
drnd = np.zeros_like(dsis)

for k in range(nCells):
    cell1 = rng.binomial(divStates[:,k].astype('int'),0.5)
    cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    dsis[k] = divStates[:,k] - 2*cell1
    drnd[k] = cell1 - cell2

vardrnd = np.var(drnd,axis=0)
vardsis = np.var(dsis,axis=0)
normvar = 1-vardsis/vardrnd

with open(PROJECT_DIR / 'production/prodAprodBsweep1/prodAprodBsweep1_PprodA_%.2i_PprodB_%.2i.pickle' % (PprodAindex,PprodBindex),'wb') as f:
    pickle.dump([means,variances,divStates,dsis,drnd,vardrnd,vardsis,normvar],f,pickle.HIGHEST_PROTOCOL)
