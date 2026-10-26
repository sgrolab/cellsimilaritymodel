# TODO: run this simulation 

# 3 Step Cascade Time Simulation 
from datetime import datetime 
import numpy as np
from las_model.simulations.figure_06 import parameters as p
from las_model.utils.cell import Cell
from las_model.utils.config import PROJECT_DIR
from las_model.utils.analyze import calculate_offspring_similarity_time
from las_model.utils.output import save_experiment

# Experiment metadata
metadata = {
    'experiment_name': 'cascade_time',
    'experiment_directory': 'cascade',
    'created': datetime.now().isoformat(),
    'seed': 1000,
    'nCells': 10000,
    'nCells_equilibrium': 10,
    'nCycles': p.nGenerations,      # as long as the spatial colony grows, from parameters.py
    # the same circuit as the spatial simulations, from parameters.py
    'Tcc': p.Tcc,
    'varTcc': p.varTcc,
    'circuit': p.circuit,
    'PprodA': p.PprodA,
    'kcatA': p.kcatA,
    'kcatB': p.kcatB,
}

# Pin random seed 
rng = np.random.default_rng(seed=metadata['seed'])

# Accumulate results
results = {
    'vardrnd': [],
    'vardsis': [],
    'normvar': []
}

# Intialize Mother cell and run simulation 
motherCell = Cell(metadata['Tcc'],metadata['varTcc'],rng)
motherCell.parameterize(metadata['circuit'],[metadata['PprodA'],metadata['kcatA'],metadata['kcatB']])
motherCell.equilibrate(metadata['nCells_equilibrium'])
motherCell.run(metadata['nCells'])

# Get mother states and calculate division differences 
divStates = motherCell.getMotherStates()
_, _, vardsis, vardrnd, normvar = calculate_offspring_similarity_time(motherCell,metadata,rng)

# Store results
results['vardrnd'].append(vardrnd)
results['vardsis'].append(vardsis)
results['normvar'].append(normvar)

# Stack results 
results = {k: np.stack(v,axis=0) for k, v in results.items()}

# Save results 
exp_dir = save_experiment(
    experiment_name=metadata['experiment_name'],
    data = results,
    metadata=metadata,
    base_dir=PROJECT_DIR / metadata['experiment_directory']
)
print(f"Experiment saved to f{exp_dir}")
