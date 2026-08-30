# Saturated Production sweep PprodA
import numpy as np 
from datetime import datetime
from las_model.utils import motiffunc as mf
from las_model.utils.config import PROJECT_DIR
from las_model.utils.output import save_experiment

# Experiment metadata
metadata = {
    'experiment_name': 'satprod_PprodAsweep',
    'experiment_directory': 'satprod/static_PprodA_sweep',
    'created': datetime.now().isoformat(),
    'seed': 1000,
    'nCells': 100,
    'nCells_equilibrium': 10,
    'Tcc': 1000,
    'circuit': 'prodsat',
    'PprodAs': list(np.logspace(-3,2,6)),
    'kcatA': 10**-3,
}

# Pin random seed 
rng = np.random.default_rng(seed=metadata['seed'])

# === Iterate over PprodA values and simulate mother cells =========
for PprodAindex, PprodA in enumerate(metadata['PprodAs']):

    print(f"Simulating for PprodAindex={PprodAindex}, PprodA={PprodA}")

    motherCell = mf.Cell(metadata['Tcc'],0,rng)
    motherCell.parameterize(metadata['circuit'],[PprodA,metadata['kcatA']])
    motherCell.equilibrate(metadata['nCells_equilibrium'])

    # Run simulation 
    motherCell.run(metadata['nCells'])

    divStates = motherCell.getMotherStates()

    motherAs = divStates[0]
    motherBs = divStates[1]

    dsis = np.zeros([metadata['nCells'],6])
    drnd = np.zeros_like(dsis)

    for k in range(metadata['nCells']):
        cell1 = rng.binomial(divStates[:,k].astype('int'),0.5)
        cell2 = rng.binomial(divStates[:,rng.integers(0,metadata['nCells'])].astype('int'),0.5)
        
        dsis[k] = divStates[:,k] - 2*cell1
        drnd[k] = cell1 - cell2

    exp_dir = save_experiment(
        experiment_name=f"{metadata['experiment_name']}_PprodAindex_{PprodAindex}",
        data=[PprodA,metadata['Tcc'],metadata['kcatA'],motherAs,motherBs,drnd,dsis],
        metadata=metadata,
        base_dir=PROJECT_DIR / metadata['experiment_directory'],
    )
    print(f"Experiment saved to {exp_dir}")
