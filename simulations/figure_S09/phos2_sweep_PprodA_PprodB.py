# Phosphorylation Motif 2: Sweep PprodA and PprodB
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

kk = 0.01
krevk = 0.01
k1 = 0.0001
krev1 = 0.01
k2 = 0.0001
krev2 = 0.01
kt = 0.01
kp = 0.01

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('phos2',[prodA,prodB,kk,krevk,k1,krev1,kt,k2,krev2,kp])
motherCell.run(nCells)

Aeq = np.mean(motherCell.A/motherCell.V)
Beq = np.mean(motherCell.B/motherCell.V)
Ceq = np.mean(motherCell.C/motherCell.V)
Deq = np.mean(motherCell.D/motherCell.V)
Eeq = np.mean(motherCell.E/motherCell.V)
Feq = np.mean(motherCell.F/motherCell.V)

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
normvarF = 1-np.var(dsis[:,5],axis=0)/np.var(drnd[:,5],axis=0)

with open(PROJECT_DIR / 'prodRateSat/phos2_0/phos2_0_prodA_%.2i_prodB_%.2i.pickle' % (prodAindex,prodBindex),'wb') as f:
    pickle.dump([prodA,prodB,Aeq,Beq,Ceq,Deq,Eeq,Feq,normvarA,
                 normvarB,normvarC,normvarD,normvarE,normvarF],
                f,pickle.HIGHEST_PROTOCOL)
