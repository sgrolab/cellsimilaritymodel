# Unsaturated production sweep substrate production rate 
import sys 
import pickle
import numpy as np 
from utils import motiffunc as mf
from utils.config import PROJECT_DIR

nCells = 1000
rng = np.random.default_rng(seed=1000)

Tcc = 1000
PprodA = 10**-1
PprodBs = np.logspace(-2,4,31)
PprodBindex = int(sys.argv[1])
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

Avar = np.var(motherCell.A/motherCell.V)
Bvar = np.var(motherCell.B/motherCell.V)
Cvar = np.var(motherCell.C/motherCell.V)

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

with open(PROJECT_DIR / 'production/prodBsweep5/prodBsweep5_PprodB_%.2i.pickle' % PprodBindex,'wb') as f:
    pickle.dump([Aeq,Beq,Ceq,Avar,Bvar,Cvar,vardrnd,vardsis,normvar],f,pickle.HIGHEST_PROTOCOL)
