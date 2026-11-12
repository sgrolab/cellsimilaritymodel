# Cascade: scramble inherited B only
from datetime import datetime
import numpy as np
from las_model.utils.cell import Cell
from las_model.utils.config import PROJECT_DIR
from las_model.utils.analyze import simulate_offspring_scramble_time, calculate_offspring_differences
from las_model.utils.output import save_experiment

# Experiment metadata shared by all runs
metadata = {
    'experiment_directory': 'cascade_scramble',
    'created': datetime.now().isoformat(),
    'seed': 1000,
    'nCells': 10000,
    'nCells_equilibrium': 10,
    'nCycles': 10,
    'Tcc': 1000,
    'varTcc': 0,
    'circuit': 'cascade',
    'PprodA': 10**-3,
    'kcatA': 10**-2,
    'kcatB': 10**-2,
    'scrambled_molecules': [1],   # B only
}

# Three runs, each saved as its own experiment, which differ in where the scrambled molecules
# come from (see simulate_offspring_scramble_time). Every run starts from the same seed, so
# all three scramble the same mothers.
runs = [
    {'experiment_name': 'cascade_scramble_B_only',
     'scramble_source': 'different_mothers'},
    {'experiment_name': 'cascade_scramble_B_only_same_mother',
     'scramble_source': 'same_mother'},
    {'experiment_name': 'cascade_scramble_B_only_mean_mother',
     'scramble_source': 'mean_mother'},
]

for run in runs:
    run_metadata = {**metadata, **run}
    print(f"Simulating {run_metadata['experiment_name']}")

    # Pin random seed
    rng = np.random.default_rng(seed=run_metadata['seed'])

    # Initialize and run Mother Cell
    motherCell = Cell(run_metadata['Tcc'],run_metadata['varTcc'],rng)
    params = [run_metadata['PprodA'],run_metadata['kcatA'],run_metadata['kcatB']]
    motherCell.parameterize(run_metadata['circuit'],params)
    motherCell.equilibrate(run_metadata['nCells_equilibrium'])
    motherCell.run(run_metadata['nCells'])

    # Simulate inherited and scrambled offspring, then pairwise differences over time
    sis1, sis2, rnd1, sis1scr, sis2scr, rnd1scr = simulate_offspring_scramble_time(
        motherCell, run_metadata, rng)

    results = {}
    for pairing, cells in [('inherited', (sis1, sis2, rnd1)),
                           ('scrambled', (sis1scr, sis2scr, rnd1scr))]:
        _, _, vardsis, vardrnd, normvar = calculate_offspring_differences(*cells)
        results[pairing] = {'vardsis': vardsis, 'vardrnd': vardrnd, 'normvar': normvar}

    # Save results
    exp_dir = save_experiment(
        experiment_name=run_metadata['experiment_name'],
        data = results,
        metadata=run_metadata,
        base_dir=PROJECT_DIR / run_metadata['experiment_directory']
    )
    print(f"Experiment saved to {exp_dir}")
