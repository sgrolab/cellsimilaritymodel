# Assymetric partitioning simulation 
import pickle
import numpy as np 
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR

# Pin random seed 
rng = np.random.default_rng(seed=1000)

# Initialize cell parameters 
Tcc = 1000
circuit = 'single'
prodA = 10**-1
biases = np.linspace(0.55,1,10)

# Set simulation parameters
nCells = 1000
nCycles_equilibrate = 10

vardsis = []
vardrnd = []
normvar = []

# === Iterate over biases and simulate mother cells ==========================================
for bias in biases:

    print(f"Simulating bias: {bias}")

    motherCell = mf.Cell(Tcc,0)
    motherCell.parameterize('single',[prodA])
    motherCell.equilibrate(nCycles_equilibrate)

    # Run simulation 
    motherCell.run(nCells,'asymmetric',bias)

    # Save mother cell state 
    divStates = motherCell.getMotherStates()

    dsis = np.zeros([6,nCells])
    drnd = np.zeros_like(dsis)

    for k in range(nCells):

        # Choose bias direction 
        cell1_bias = rng.integers(2)
        cell2_bias = rng.integers(2)

        # Get mother cell state 
        motherCellState = divStates[:,k].astype('int')
        # print(f'mother cell state: {motherCellState}')
        
        # Partition sister cells 
        if cell1_bias==0:
            cell1 = (motherCellState * bias).astype('int')
        else:
            cell1 = (motherCellState * (1-bias)).astype('int')

        # Get random cell state 
        rndCellIndx = rng.integers(0,nCells)
        rndCellState = divStates[:,rndCellIndx].astype('int')

        # Partition random cells 
        if cell2_bias == 0:
            cell2 = (rndCellState * bias).astype('int')
        else:
            cell2 = (rndCellState * (1-bias)).astype('int')

        # Compute pairwise differences
        dsis[:,k] = motherCellState - 2*cell1
        drnd[:,k] = cell1 - cell2

    # Compute variance of pairwise differences
    vardsis.append(np.var(dsis, axis=1))
    vardrnd.append(np.var(drnd, axis=1))
    normvar.append(1-vardsis[-1][0]/vardrnd[-1][0])

with open(PROJECT_DIR / 'asymmetric/asymmetric_bias_screen_2.pickle','wb') as f:
    pickle.dump([biases,vardsis,vardrnd,normvar],f,pickle.HIGHEST_PROTOCOL)
