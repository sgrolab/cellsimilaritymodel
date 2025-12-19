# 3 Step Cascade Time Simulation 
import pickle
import numpy as np 
from utils import motiffunc as mf
from utils.config import PROJECT_DIR

nCells = 1000
nCycles = 10
rng = np.random.default_rng(seed=1000)

Tcc = 1000
prodA = 10**-1
ka = 10**-1
kb = 10**-1

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('cascade',[prodA,ka,kb])
motherCell.equilibrate(20)
motherCell.run(nCells)

divStates = motherCell.motherStates[2::]

concentrations = np.zeros([3,nCells,6,int(nCycles*Tcc/10+1)])

for i in range(nCells):
    
    sis1state = rng.binomial(divStates[:,i].astype('int'),0.5)
    sis2state = divStates[:,i] - sis1state
    rnd1state = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    sis1 = mf.Cell(Tcc,0)
    sis1.inherit(motherCell,sis1state)
    sis1.run(nCycles)
    concentrations[0,i] = sis1.molecules
    
    sis2 = mf.Cell(Tcc,0)
    sis2.inherit(motherCell,sis2state)
    sis2.run(nCycles)
    concentrations[1,i] = sis2.molecules
    
    rnd1 = mf.Cell(Tcc,0)
    rnd1.inherit(motherCell,rnd1state)
    rnd1.run(nCycles)
    concentrations[2,i] = rnd1.molecules

with open(PROJECT_DIR / 'cascade/motifs_cascade.pickle','wb') as f:
    pickle.dump([motherCell,concentrations],f,pickle.HIGHEST_PROTOCOL)
