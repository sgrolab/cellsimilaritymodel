# TODO: run with full PprodA range 

# Saturated Production sweep PprodA
import numpy as np 
from datetime import datetime
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR
from las_model.utils.output import save_experiment
from las_model.utils.analyze import calculate_division_differences

# Experiment metadata
metadata = {
    'experiment_name': 'satprod_PprodAsweep',
    'experiment_directory': 'satprod',
    'created': datetime.now().isoformat(),
    'seed': 1000,
    'nCells': 1000,
    'nCells_equilibrium': 10,
    'Tcc': 1000,
    'varTcc': 0,
    'circuit': 'prodsat',
    'PprodAs': list(np.logspace(-3,2,6)),
    'kcatA': 10**-1,
}

# Pin random seed 
rng = np.random.default_rng(seed=metadata['seed'])

# Accumulate results 
results = {
    'mother_As': [],
    'mother_Bs': [],
    'dsis': [],
    'drnd': [],
    'vardsis': [],
    'vardrnd': [],
    'normvar': [],
}

# === Iterate over PprodA values and simulate mother cells =========
for PprodAindex, PprodA in enumerate(metadata['PprodAs']):

    print(f"Simulating for PprodAindex={PprodAindex}, PprodA={PprodA}")

    motherCell = mf.Cell(metadata['Tcc'],metadata['varTcc'],rng)
    motherCell.parameterize(metadata['circuit'],[PprodA,metadata['kcatA']])
    motherCell.equilibrate(metadata['nCells_equilibrium'])

    # Run simulation 
    motherCell.run(metadata['nCells'])

    divStates = motherCell.getMotherStates()
    dsis, drnd, vardsis, vardrnd, normvar = calculate_division_differences(divStates,rng)

    results['mother_As'].append(divStates[0])
    results['mother_Bs'].append(divStates[1])
    results['dsis'].append(dsis)
    results['drnd'].append(drnd)
    results['vardsis'].append(vardsis)
    results['vardrnd'].append(vardrnd)
    results['normvar'].append(normvar)

# Stack results 
results = {k: np.stack(v,axis=0) for k, v in results.items()}

exp_dir = save_experiment(
    experiment_name=metadata['experiment_name'],
    data=[metadata['PprodAs'],results],
    metadata=metadata,
    base_dir=PROJECT_DIR / metadata['experiment_directory'],
)
print(f"Experiment saved to {exp_dir}")
