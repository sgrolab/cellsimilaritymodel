# Saturated Production: Sweep Tcc (low)
from datetime import datetime 
import numpy as np
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR
from las_model.utils.analyze import calculate_division_differences
from las_model.utils.output import save_experiment 

# Experiment metadata
metadata = {
    'experiment_name': 'prodsat_sweep_Tcc_low',
    'experiment_directory': 'satprod',
    'created': datetime.now().isoformat(),
    'seed': 1000,
    'nCells': 1000,
    'nCells_equilibrium': 10,
    'Tccs': list(np.logspace(2,4,5)),
    'varTcc': 0,
    'circuit': 'prodsat',
    'PprodA': 10**-2,
    'kcatA': 10**-2,
}

# Pin random seed 
rng = np.random.default_rng(seed=metadata['seed'])

# Accumulate results 
results = {
    'dsis': [],
    'drnd': [],
    'vardsis': [],
    'vardrnd': [],
    'normvar': [],
}

# === Iterate over Tcc values and simulate cells ==========
for Tcc in metadata['Tccs']:

    print(f"Simulating for Tcc = {Tcc}")

    motherCell = mf.Cell(Tcc,metadata['varTcc'],rng)
    motherCell.parameterize(metadata['circuit'],[metadata['PprodA'],metadata['kcatA']])
    motherCell.equilibrate(metadata['nCells_equilibrium'])

    motherCell.run(metadata['nCells'])

    # Get mother states and calculate division differences 
    divStates = motherCell.getMotherStates()
    dsis, drnd, vardsis, vardrnd, normvar = calculate_division_differences(divStates,rng)

    results['dsis'].append(dsis)
    results['drnd'].append(drnd)
    results['vardsis'].append(vardsis)
    results['vardrnd'].append(vardrnd)
    results['normvar'].append(normvar)

# Stack results 
results = {k: np.stack(v,axis=0) for k, v in results.items()}

# Save results 
exp_dir = save_experiment(
    experiment_name=metadata['experiment_name'],
    data = [metadata['Tccs'],results],
    metadata=metadata,
    base_dir=PROJECT_DIR / metadata['experiment_directory']
)
print(f"Experiment saved to f{exp_dir}")

# nCells = 1000
# rng = np.random.default_rng(seed=1000)

# Tccs = np.logspace(2,4,5)
# kcatA = 10**-2
# PprodA = 10**-2

# motherCells = []
# divStates = np.zeros([len(Tccs),5,nCells])

# for i in range(len(Tccs)):

#     motherCell = mf.Cell(Tccs[i],0)
#     motherCell.parameterize('prodsat',[PprodA,kcatA])
#     motherCell.run(nCells)
    
#     motherCells.append(motherCell)
    
#     divStates[i] = motherCell.getMotherStates()
    
# dsis = np.zeros([len(Tccs),nCells,5])
# drnd = np.zeros([len(Tccs),nCells,5])

# for j in range(len(Tccs)):
#     for i in range(nCells):
#         cell1 = rng.binomial(divStates[j,:,i].astype('int'),0.5)
#         cell2 = rng.binomial(divStates[j,:,rng.integers(0,nCells)].astype('int'),0.5)
        
#         dsis[j,i] = divStates[j,:,i] - 2*cell1
#         drnd[j,i] = cell1 - cell2

# with open(PROJECT_DIR / 'analyticalData/motifs_prodsat_Tccsweep.pickle','wb') as f:
#     pickle.dump([Tccs,motherCells,dsis,drnd],f,pickle.HIGHEST_PROTOCOL)
