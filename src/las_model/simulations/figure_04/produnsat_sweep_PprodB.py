# Unsaturated production sweep substrate production rate 
from datetime import datetime 
import numpy as np
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR
from las_model.utils.analyze import calculate_division_differences
from las_model.utils.output import save_experiment 

# Experiment metadata
metadata = {
    'experiment_name': 'produnsat_sweep_PprodB',
    'experiment_directory': 'production',
    'created': datetime.now().isoformat(),
    'seed': 1000,
    'nCells': 1000,
    'nCells_equilibrium': 10,
    'Tcc': 1000,
    'varTcc': 0,
    'circuit': 'produnsat',
    'PprodA': 10**-1,
    'kcatA': 10**-1,
    # 'PprodBs': list(np.logspace(-2,4,31)), commented out for small test 
    'PprodBs': list(np.logspace(-2,-1,2)),
    'Km': 10**3
}

# Pin random seed 
rng = np.random.default_rng(seed=metadata['seed'])

# Accumulate results
results = {
    'means': [],
    'variances': [],
    'dsis': [],
    'drnd': [],
    'vardrnd': [],
    'vardsis': [],
    'normvar': []
}

for PprodB in metadata['PprodBs']:

    print(f"Running simulation for PprodB={PprodB}")

    # Initialize and run mother cell 
    motherCell = mf.Cell(metadata['Tcc'],metadata['varTcc'],rng)
    motherCell.parameterize(metadata['circuit'],[PprodB,metadata['PprodA'],metadata['kcatA'],metadata['Km']])
    motherCell.equilibrate(metadata['nCells_equilibrium'])
    motherCell.run(metadata['nCells'])

    # Get molecule amounts 
    molecules = motherCell.getMolecules()
    means = [np.mean(molecules[0]), np.mean(molecules[1]), np.mean(molecules[2])]
    variances = [np.var(molecules[0]), np.var(molecules[1]), np.var(molecules[2])]

    # Get mother states and calculate division differences 
    divStates = motherCell.getMotherStates()
    dsis, drnd, vardsis, vardrnd, normvar = calculate_division_differences(divStates,rng)

    # Store results
    results['means'].append(means)
    results['variances'].append(variances)
    results['dsis'].append(dsis)
    results['drnd'].append(drnd)
    results['vardrnd'].append(vardrnd)
    results['vardsis'].append(vardsis)
    results['normvar'].append(normvar)

# Stack results 
results = {k: np.stack(v,axis=0) for k, v in results.items()}

# Save results 
exp_dir = save_experiment(
    experiment_name=metadata['experiment_name'],
    data = [metadata['PprodBs'],results],
    metadata=metadata,
    base_dir=PROJECT_DIR / metadata['experiment_directory']
)
print(f"Experiment saved to f{exp_dir}")


# nCells = 1000
# rng = np.random.default_rng(seed=1000)

# Tcc = 1000
# PprodA = 10**-1
# PprodBs = np.logspace(-2,4,31)
# PprodBindex = int(sys.argv[1])
# PprodB = PprodBs[PprodBindex]
# kcatA = 10**-1
# Km = 10**3

# motherCell = mf.Cell(Tcc,0)
# motherCell.parameterize('produnsat',[PprodB,PprodA,kcatA,Km])
# motherCell.equilibrate()
# motherCell.run(nCells)

# Aeq = np.mean(motherCell.A/motherCell.V)
# Beq = np.mean(motherCell.B/motherCell.V)
# Ceq = np.mean(motherCell.C/motherCell.V)

# Avar = np.var(motherCell.A/motherCell.V)
# Bvar = np.var(motherCell.B/motherCell.V)
# Cvar = np.var(motherCell.C/motherCell.V)

# divStates = motherCell.getMotherStates()

# dsis = np.zeros([nCells,6])
# drnd = np.zeros_like(dsis)

# for k in range(nCells):
#     cell1 = rng.binomial(divStates[:,k].astype('int'),0.5)
#     cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
#     dsis[k] = divStates[:,k] - 2*cell1
#     drnd[k] = cell1 - cell2

# vardrnd = np.var(drnd,axis=0)
# vardsis = np.var(dsis,axis=0)
# normvar = 1-vardsis/vardrnd

# with open(PROJECT_DIR / 'production/prodBsweep5/prodBsweep5_PprodB_%.2i.pickle' % PprodBindex,'wb') as f:
#     pickle.dump([Aeq,Beq,Ceq,Avar,Bvar,Cvar,vardrnd,vardsis,normvar],f,pickle.HIGHEST_PROTOCOL)
