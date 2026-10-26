# Parameters of the spatial cascade simulations, the single source of truth for
# figure_06/cascade_spatial.py, cascade_spatial_moranI.py, cascade_spatial_analyze.py and
# figure_S15/cascade_spatial_seeds.py. Each script copies what it uses into its saved
# metadata, and plot_figures.py reads it from there, so a figure always matches its data.
# Not a simulation itself: run_simulations.py skips it.
import math

# Colony size: the grid stops growing at maxCells cells, i.e. after nGenerations generations
maxCells = 2**12
nGenerations = int(math.log2(maxCells))
assert maxCells == 2**nGenerations, 'maxCells must be a power of two'

# The colony grows as a rough disc of radius sqrt(maxCells/pi) around the grid centre, and
# its farthest cell sits about 1.4 radii out (1.33 and 1.39 in two 1024-cell runs). The grid
# leaves 2 radii on each side so the colony never reaches the edge, and the concentration and
# cousin maps are cropped to 1.4 radii (25 cells for 1024 cells, as before).
colonyRadius = math.sqrt(maxCells/math.pi)
gridSize = 2*math.ceil(2*colonyRadius) + 1
cropHalfWidth = round(1.4*colonyRadius)

# Cascade circuit
seed = 1000
Tcc = 1000
varTcc = 10
circuit = 'cascade'
PprodA = 10**-1
kcatA = 10**-2
kcatB = 10**-2

# Figure S15: grids rerun with these seeds
seeds = list(range(1000, 1008))

# Relatedness curve (figures 06 and S15), at the end of the last generation
maxRadius = 9

# Moran's I over time (figures 06 and S16)
neighborhoodsizes = list(range(1, 10))
shapes = ['discdist', 'discstep', 'donut', 'gausdist']
moranI_timestep = 100
molecules = ['A', 'B', 'C']

# Analysis for figures 06 and S14. Generations count from the first division, at time
# generation*Tcc.
panelB_generations = [6, 8, 10, 12]     # figure 06 panel B: the concentration maps
cousinCell = 100
cousin_generations = list(range(nGenerations-8, nGenerations+1, 2))
cousinCells_all = [100, 200, 300, 400, 500]     # figure S14 rows
cousin_generations_all = list(range(nGenerations-10, nGenerations+1, 2))     # figure S14 columns
moranI_shape = 'discdist'     # figure 06 panel D
