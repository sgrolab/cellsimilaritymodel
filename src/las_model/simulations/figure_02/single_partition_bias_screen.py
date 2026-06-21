# Assymetric partitioning simulation 
import pickle
import numpy as np 
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR

Tcc = 1000
nCells = 1000
rng = np.random.default_rng(seed=1000)

biases = np.linspace(0.55,1,10)
prodA = 10**-1

Aeqs = np.zeros(len(biases))
normvarAs = np.zeros_like(Aeqs)

for j in range(len(biases)):

    motherCell = mf.Cell(Tcc,0)
    motherCell.parameterize('single',[prodA])
    motherCell.run(nCells,'asymmetric',biases[j])
    
    Aeqs[j] = np.mean(motherCell.A/motherCell.V)
    
    divStates = motherCell.getMotherStates()
    
    dsis = np.zeros([nCells,6])
    drnd = np.zeros([nCells,6])
    
    for k in range(nCells):
        if rng.integers(2)==0:
            cell1 = divStates[:,k].astype('int') * biases[j]
        else:
            cell1 = divStates[:,k].astype('int') * (1-biases[j])
        
        if rng.integers(2) == 0:
            cell2 = divStates[:,rng.integers(0,nCells)].astype('int') * biases[j]
        else:
            cell2 = divStates[:,rng.integers(0,nCells)].astype('int') * (1-biases[j])
        
        dsis[k] = divStates[:,k] - 2*cell1
        drnd[k] = cell1 - cell2
    
    normvarAs[j] = 1-np.var(dsis[:,0],axis=0)/np.var(drnd[:,0],axis=0)


with open(PROJECT_DIR / 'asymmetric/asymmetric_bias_screen.pickle','wb') as f:
    pickle.dump([biases,Aeqs,normvarAs],f,pickle.HIGHEST_PROTOCOL)
