# 3 Step Cascade Spatial Simulation
from datetime import datetime
import numpy as np
from las_model.simulations.figure_06 import parameters as p
from las_model.utils import gridfunc as gf
from las_model.utils.config import PROJECT_DIR
from las_model.utils.output import save_experiment

# Experiment metadata, from parameters.py
metadata = {
    'experiment_name': 'cascade_spatial',
    'experiment_directory': 'cascade',
    'created': datetime.now().isoformat(),
    'seed': p.seed,
    'maxCells': p.maxCells,
    'gridSize': p.gridSize,
    'Tcc': p.Tcc,
    'varTcc': p.varTcc,
    'circuit': p.circuit,
    'PprodA': p.PprodA,
    'kcatA': p.kcatA,
    'kcatB': p.kcatB,
}

# Pin random seed 
rng = np.random.default_rng(seed=metadata['seed'])

# Initiate and seed grid
grid = gf.Grid(metadata['gridSize'],metadata['gridSize'],metadata['maxCells'],rng)
grid.seed(metadata['circuit'],[metadata['PprodA'],metadata['kcatA'],metadata['kcatB']],metadata['Tcc'],metadata['varTcc'])

# Run grid
grid.run()

# Save results
exp_dir = save_experiment(
    experiment_name=metadata['experiment_name'],
    data=grid,
    metadata=metadata,
    base_dir=PROJECT_DIR / metadata['experiment_directory']
)
print(f"Experiment saved to {exp_dir}")
