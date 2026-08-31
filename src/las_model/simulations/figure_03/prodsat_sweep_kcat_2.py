# Binding motif initial values 
from datetime import datetime 
import numpy as np
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR
from las_model.utils.analyze import calculate_division_differences
from las_model.utils.output import save_experiment 

# Experiment metadata
metadata = {
    'experiment_name': 'prodsat_sweep_kcat_high',
    'experiment_directory': 'satprod',
    'created': datetime.now().isoformat(),
    'seed': 1000,
    'nCells': 1000,
    'nCells_equilibrium': 10,
    'Tcc': 1000,
    'varTcc': 0,
    'circuit': 'prodsat',
    'PprodA': 10**-1,
    'kcats': list(np.logspace(-3,0,7)),
}

# Pin random seed 
rng = np.random.default_rng(metadata['seed'])

# Accumulate results 
results = {
    'dsis': [],
    'drnd': [],
    'vardsis': [],
    'vardrnd': [],
    'normvar': [],
}

# === Iterate over kcatA values and simulate cells ==========
for kcatA in metadata['kcats']:

    print(f"Simulating for kcatA = {kcatA}")

    motherCell = mf.Cell(metadata['Tcc'],metadata['varTcc'],rng)
    motherCell.parameterize(metadata['circuit'],[metadata['PprodA'],kcatA])
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
    data = [metadata['kcats'],results],
    metadata=metadata,
    base_dir=PROJECT_DIR / metadata['experiment_directory']
)
print(f"Experiment saved to f{exp_dir}")

# nCells = 1000
# rng = np.random.default_rng(seed=1000)

# Tcc = 1000
# PprodA = 10**-1
# kcats = np.logspace(-3,3,13)
# kcatAindex = int(sys.argv[1])
# kcatA = kcats[kcatAindex]

# motherCells = []
# divStates = np.zeros([len(kcats),6,nCells])

# motherCell = mf.Cell(Tcc,0)
# motherCell.parameterize('prodsat',[PprodA,kcatA])
# motherCell.equilibrate(20)
# motherCell.run(nCells)

# motherCells.append(motherCell)

# divStates = motherCell.motherStates[2::]

# dsis = np.zeros([nCells,6])
# drnd = np.zeros([nCells,6])

# for i in range(nCells):
#     cell1 = rng.binomial(divStates[:,i].astype('int'),0.5)
#     cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
#     dsis[i] = divStates[:,i] - 2*cell1
#     drnd[i] = cell1 - cell2

# normvarA = 1-np.var(dsis[:,0],axis=0)/np.var(drnd[:,0],axis=0)
# normvarB = 1-np.var(dsis[:,1],axis=0)/np.var(drnd[:,1],axis=0)

# with open(PROJECT_DIR / 'prodsat_sweep/prodsat_kcatsweep/prodsat_kcatsweep_kcatA_%.2i.pickle' % (kcatAindex),'wb') as f:
#     pickle.dump([kcatA,normvarA,normvarB],f,pickle.HIGHEST_PROTOCOL)
