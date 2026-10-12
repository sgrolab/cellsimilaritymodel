# Saturated Production: Sweep kcatA (high)
from datetime import datetime
import numpy as np
from las_model.utils.cell import Cell
from las_model.utils.config import PROJECT_DIR
from las_model.utils.analyze import calculate_division_differences
from las_model.utils.output import save_experiment
from las_model.utils.parallel import run_pool

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
    'kcats': list(np.logspace(-3,3,13)),
}

def _simulate_single_kcat(task_args):
    """Worker task simulating one kcatA value."""
    idx, kcatA, meta, seed = task_args
    task_rng = np.random.default_rng(seed)

    print(f"Simulating for kcatA = {kcatA}")

    motherCell = Cell(meta['Tcc'], meta['varTcc'], task_rng)
    motherCell.parameterize(meta['circuit'], [meta['PprodA'], kcatA])
    motherCell.equilibrate(meta['nCells_equilibrium'])
    motherCell.run(meta['nCells'])

    # Get mother states and calculate division differences
    divStates = motherCell.getMotherStates()
    dsis, drnd, vardsis, vardrnd, normvar = calculate_division_differences(divStates, task_rng)

    return {
        'dsis': dsis,
        'drnd': drnd,
        'vardsis': vardsis,
        'vardrnd': vardrnd,
        'normvar': normvar,
    }

if __name__ == '__main__':
    # Generate statistically independent, deterministic seeds for each parameter
    master_rng = np.random.default_rng(metadata['seed'])
    seeds = master_rng.integers(0, 2**63 - 1, size=len(metadata['kcats']))

    tasks = [
        (idx, kcatA, metadata, seeds[idx])
        for idx, kcatA in enumerate(metadata['kcats'])
    ]

    sweep_results = run_pool(_simulate_single_kcat, tasks, desc=metadata['experiment_name'],
                             sort_key=lambda task: task[1])   # start the costly points first

    # Stack results
    results = {
        k: np.stack([res[k] for res in sweep_results], axis=0)
        for k in sweep_results[0].keys()
    }

    # Save results
    exp_dir = save_experiment(
        experiment_name=metadata['experiment_name'],
        data=[metadata['kcats'], results],
        metadata=metadata,
        base_dir=PROJECT_DIR / metadata['experiment_directory']
    )
    print(f"Experiment saved to {exp_dir}")
