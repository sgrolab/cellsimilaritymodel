# Two Component System Time Simulation  
import pickle
import numpy as np 
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR

nCells = 1000
nCycles = 10
rng = np.random.default_rng(seed=1000)

Tcc = 1000
prodA = 10**-1
prodB = 10**1
ka = 10**-2
kb = 10**-4
kt = 10**-2
k4 = 10**-2

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('tcs',[prodA,prodB,ka,kb,kt,k4])
motherCell.equilibrate()
motherCell.run(nCells)

divStates = motherCell.getMotherStates()

concentrations = np.zeros([3,nCells,6,int(nCycles*Tcc/10+1)])

for i in range(nCells):
    
    sis1state = rng.binomial(divStates[:,i].astype('int'),0.5)
    sis2state = divStates[:,i] - sis1state
    rnd1state = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    sis1 = mf.Cell(Tcc,0)
    sis1.inherit(motherCell,sis1state)
    sis1.run(nCycles)
    concentrations[0,i] = sis1.getMolecules()
    
    sis2 = mf.Cell(Tcc,0)
    sis2.inherit(motherCell,sis2state)
    sis2.run(nCycles)
    concentrations[1,i] = sis2.getMolecules()
    
    rnd1 = mf.Cell(Tcc,0)
    rnd1.inherit(motherCell,rnd1state)
    rnd1.run(nCycles)
    concentrations[2,i] = rnd1.getMolecules()

vardsis = np.var(concentrations[0] - concentrations[1],axis=0)
vardrnd = np.var(concentrations[0] - concentrations[2],axis=0)
normvar = 1-vardsis/vardrnd

with open(PROJECT_DIR / 'motifs/motifs_tcs3.pickle','wb') as f:
    pickle.dump(normvar,f,pickle.HIGHEST_PROTOCOL)
