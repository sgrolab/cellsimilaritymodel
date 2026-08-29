# Saturated Production Time Run
import pickle
import sys
import numpy as np 

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR

Tcc = 1000
nCells = 1000
nCycles = 10
rng = np.random.default_rng(seed=1000)

PprodA = 10**-1
kcatA = 10**-1

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('prodsat',[PprodA,kcatA])
motherCell.equilibrate(10)
motherCell.run(nCells)

divStates = motherCell.getMotherStates()

concentrations = np.zeros([3,nCells,6,int(nCycles*Tcc/10+1)])

for i in range(nCells):
    print(f'Simulating cell {i}...')
    
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

# vardsis = np.var(concentrations[0] - concentrations[1],axis=0)
# vardrnd = np.var(concentrations[0] - concentrations[2],axis=0)
# normvar = 1-vardsis/vardrnd

with open(PROJECT_DIR / 'correlationcoef/prodsat_time.pickle','wb') as f:
    pickle.dump(concentrations,f,pickle.HIGHEST_PROTOCOL)
