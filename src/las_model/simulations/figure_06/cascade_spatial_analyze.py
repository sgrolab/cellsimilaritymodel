# 3 Step Cascade Spatial Simulation: relatedness, cousin maps, concentration maps, Moran's I
import pickle
from datetime import datetime 
import numpy as np
from las_model.simulations.figure_06 import parameters as p
from las_model.utils.config import PROJECT_DIR
from las_model.utils.output import save_experiment

# Experiment metadata, from parameters.py. Generations become times generation*Tcc.
centre = p.gridSize//2
metadata = {
    'experiment_name': 'cascade_spatial_analyze',
    'experiment_directory': 'cascade',
    'created': datetime.now().isoformat(),
    'source_experiment': 'cascade_spatial',
    'moranI_experiment': 'cascade_spatial_moranI',
    'maxCells': p.maxCells,
    'nGenerations': p.nGenerations,
    'Tcc': p.Tcc,
    'maxRadius': p.maxRadius,
    'relatedness_t': p.nGenerations*p.Tcc,
    'crop': [centre-p.cropHalfWidth, centre+p.cropHalfWidth],
    'cousinCell': p.cousinCell,
    'cousin_generations': p.cousin_generations,
    'cousin_ts': [g*p.Tcc for g in p.cousin_generations],
    'cousinCells_all': p.cousinCells_all,
    'cousin_generations_all': p.cousin_generations_all,
    'cousin_ts_all': [g*p.Tcc for g in p.cousin_generations_all],
    'cousinMax': p.nGenerations-2,      # cousin maps hold cousinMax minus the cousin number
    'molecules': p.molecules,
    'panelB_generations': p.panelB_generations,
    'mol_ts': [g*p.Tcc for g in p.panelB_generations],
    'moranI_shape': p.moranI_shape,
    'moranI_radii': p.neighborhoodsizes,
}

base_dir = PROJECT_DIR / metadata['experiment_directory']
crop = slice(metadata['crop'][0],metadata['crop'][1])

# Load grid from the spatial simulation
with open(base_dir / metadata['source_experiment'] / f"{metadata['source_experiment']}.pickle",'rb') as f:
    grid = pickle.load(f)

if grid.maxCells != p.maxCells:
    raise ValueError(f"{metadata['source_experiment']} grew to {grid.maxCells} cells but parameters.py "
                     f"has maxCells = {p.maxCells}: rerun cascade_spatial.py")
late = [t for t in metadata['cousin_ts'] + metadata['cousin_ts_all'] + metadata['mol_ts'] if t >= grid.timepoints[-1]]
if late:
    raise ValueError(f"times {late} are past the end of {metadata['source_experiment']} "
                     f"(t = {grid.timepoints[-1]:.0f}): check the generations in parameters.py")

# Collect Moran's I across neighborhood sizes (cheap; fail early if any run is missing)
morIs = []
for r in metadata['moranI_radii']:
    name = f"{metadata['moranI_experiment']}_{metadata['moranI_shape']}_r{r}"
    with open(base_dir / name / f'{name}.pickle','rb') as f:
        morIs.append(pickle.load(f))
morIs = np.stack(morIs,axis=0)

# Compute relatedness curve 
relatedness = grid.calcCollectiveLocalRelatedness(metadata['maxRadius'],metadata['relatedness_t'])

# Pull cousin maps: cousinMax minus the cousin number, NaN where there is no cell
# (grid.cousinMap marks those -2)
def cousinMap(cellNum,t):
    img = metadata['cousinMax']-grid.cousinMap(cellNum,t)[crop,crop]
    img[np.where(img==metadata['cousinMax']+2)] = 'NaN'
    return img

cousinmaps = []
for t in metadata['cousin_ts']:
    cousinmaps.append(cousinMap(metadata['cousinCell'],t))

cousinmaps_all = []
for cellNum in metadata['cousinCells_all']:
    cousinmaps_all.append([cousinMap(cellNum,t) for t in metadata['cousin_ts_all']])

# Pull molecular concentration maps 
molConcMaps = []
for molecule in metadata['molecules']:
    molConcMaps.append([grid.getFrame(t,molecule)[crop,crop] for t in metadata['mol_ts']])

# Save results 
results = {
    'relatedness': relatedness,
    'cousinmaps': cousinmaps,
    'cousinmaps_all': cousinmaps_all,
    'molConcMaps': molConcMaps,
    'moranIs': morIs,
}
exp_dir = save_experiment(
    experiment_name=metadata['experiment_name'],
    data=results,
    metadata=metadata,
    base_dir=base_dir
)
print(f"Experiment saved to {exp_dir}")
