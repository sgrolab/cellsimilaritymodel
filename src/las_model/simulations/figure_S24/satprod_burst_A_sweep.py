# Saturated Production: effect of burst size on LAS duration
from datetime import datetime
import numpy as np
from las_model.utils.cell import Cell
from las_model.utils.config import PROJECT_DIR
from las_model.utils.analyze import calculate_offspring_similarity_time
from las_model.utils.output import save_experiment

# Experiment metadata shared by both runs
metadata = {
    'experiment_directory': 'satprod',
    'created': datetime.now().isoformat(),
    'seed': 1000,
    'nCells': 10000,
    'nCells_equilibrium': 10,
    'nCycles': 10,
    'Tcc': 1000,
    'varTcc': 0,
    'circuit': 'prodsat_burst',
}

# Four runs, each saved as its own experiment, at kcatA*Tcc = 100, 0.1, 10 and 1. At 100 the
# LAS of B is near 1 for every burst size; the smaller kcatA*Tcc, the smaller the LAS but the
# longer it keeps growing with burst size. PprodA does not change the LAS, so the runs after
# the first set it to 0.1/(kcatA*Tcc) to keep B at ~300 molecules. New runs go at the end so
# the earlier ones keep their random draws.
runs = [
    {'experiment_name': 'satprod_burst_time_sweep', 'PprodA': 10**-2, 'kcatA': 10**-1,
     'burstSizes': [1, 2, 5, 10, 20]},
    {'experiment_name': 'satprod_burst_time_sweep_lowkappa', 'PprodA': 10**0, 'kcatA': 10**-4,
     'burstSizes': [1, 2, 5, 10, 20, 50]},
    {'experiment_name': 'satprod_burst_time_sweep_kappa10', 'PprodA': 10**-2, 'kcatA': 10**-2,
     'burstSizes': [1, 2, 5, 10, 20, 50]},
    {'experiment_name': 'satprod_burst_time_sweep_kappa1', 'PprodA': 10**-1, 'kcatA': 10**-3,
     'burstSizes': [1, 2, 5, 10, 20, 50]},
]

# Pin random seed
rng = np.random.default_rng(seed=metadata['seed'])

for run in runs:
    run_metadata = {**metadata, **run}

    # Accumulate results
    results = {
        'vardsis': [],
        'vardrnd': [],
        'normvar': [],
    }

    # === Iterate over burst sizes and simulate mother cells ==========
    for burstSize in run_metadata['burstSizes']:

        print(f"Simulating {run_metadata['experiment_name']} for burst size = {burstSize}")

        # Bursts of burstSize molecules at PprodA / burstSize, so the mean production is the same for every burst size
        motherCell = Cell(run_metadata['Tcc'], run_metadata['varTcc'], rng)
        motherCell.parameterize(run_metadata['circuit'], [run_metadata['PprodA'] / burstSize, run_metadata['kcatA'], burstSize])
        motherCell.equilibrate(run_metadata['nCells_equilibrium'])
        motherCell.run(run_metadata['nCells'])

        # Offspring similarity over time
        dsis, drnd, vardsis, vardrnd, normvar = calculate_offspring_similarity_time(motherCell, run_metadata, rng)

        results['vardsis'].append(vardsis)
        results['vardrnd'].append(vardrnd)
        results['normvar'].append(normvar)

    # Stack results
    results = {k: np.stack(v, axis=0) for k, v in results.items()}

    # Save results
    exp_dir = save_experiment(
        experiment_name=run_metadata['experiment_name'],
        data=[run_metadata['burstSizes'], results],
        metadata=run_metadata,
        base_dir=PROJECT_DIR / run_metadata['experiment_directory']
    )
    print(f"Experiment saved to {exp_dir}")
