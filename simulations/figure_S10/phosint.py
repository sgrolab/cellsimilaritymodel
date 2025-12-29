# Phosphorylation with Intermediate: Sweep PprodA and PprodB 
import sys 
import pickle
import numpy as np 
from utils import motiffunc as mf
from utils.config import PROJECT_DIR

Tcc = 1000
nCells = 2000
rng = np.random.default_rng(seed=1000)

prodAs = np.logspace(-3,2,31)
prodBs = np.logspace(-3,2,31)

prodAindex = int(sys.argv[1])
prodA = prodAs[prodAindex]
prodBindex = int(sys.argv[2])
prodB = prodBs[prodBindex]

ka = 10**-2
kb = 10**-4
kt = 10**-2

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('phos_int',[prodA,prodB,ka,kb,kt])
motherCell.equilibrate()
motherCell.run(nCells)

Aeq = np.mean(motherCell.A/motherCell.V)
Beq = np.mean(motherCell.B/motherCell.V)
Ceq = np.mean(motherCell.C/motherCell.V)
Deq = np.mean(motherCell.D/motherCell.V)
Eeq = np.mean(motherCell.E/motherCell.V) 

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
normvarE = 1-np.var(dsis[:,4],axis=0)/np.var(drnd[:,4],axis=0)

with open(PROJECT_DIR / 'prodRateSat/phos_int3/phos_int3_prodA_%.2i_prodB_%.2i.pickle' % (prodAindex,prodBindex),'wb') as f:
    pickle.dump([prodA,prodB,Aeq,Beq,Ceq,Deq,Eeq,normvarA,normvarB,normvarC,normvarD,normvarE],f,pickle.HIGHEST_PROTOCOL)


