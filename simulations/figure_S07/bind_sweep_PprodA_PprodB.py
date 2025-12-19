# Binding Motif: Sweep PprodA and PproB'
import sys 
import pickle
import numpy as np 
from utils import motiffunc as mf
from utils.config import PROJECT_DIR

nCells = 2000
rng = np.random.default_rng(seed=1000)

Tcc = 1000

prodAs = np.logspace(-3,2,31)
prodAindex = int(sys.argv[1])
prodA = prodAs[prodAindex]

prodBs = np.logspace(-3,2,31)
prodBindex = int(sys.argv[2])
prodB = prodBs[prodBindex]

k = 10**-3

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('bind2',[prodA,prodB,k])
motherCell.run(nCells)

Aeq = np.mean(motherCell.A/motherCell.V)
Beq = np.mean(motherCell.B/motherCell.V)
Ceq = np.mean(motherCell.C/motherCell.V)

varA = np.var(motherCell.A/motherCell.V)
varB = np.var(motherCell.B/motherCell.V)
varC = np.var(motherCell.C/motherCell.V)

divStates = motherCell.getMotherStates()

dsis = np.zeros([nCells,5])
drnd = np.zeros([nCells,5])

for k in range(nCells):
    cell1 = rng.binomial(divStates[:,k].astype('int'),0.5)
    cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    dsis[k] = divStates[:,k] - 2*cell1
    drnd[k] = cell1 - cell2

normvarA = 1-np.var(dsis[:,0],axis=0)/np.var(drnd[:,0],axis=0)
normvarB = 1-np.var(dsis[:,1],axis=0)/np.var(drnd[:,1],axis=0)
normvarC = 1-np.var(dsis[:,2],axis=0)/np.var(drnd[:,2],axis=0)

with open(PROJECT_DIR / 'prodRateSat/bind3/bind3_prodA_%.2i_prodB_%.2i.pickle' % (prodAindex,prodBindex),'wb') as f:
    pickle.dump([prodA,prodB,Aeq,Beq,Ceq,varA,varB,varC,normvarA,normvarB,normvarC],f,pickle.HIGHEST_PROTOCOL)


