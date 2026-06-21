# Spatial Simulation 
import pickle 
import utils.gridfunc as gf
from las_model.utils.config import PROJECT_DIR

# set parameters 
maxCells = 2**10
PprodA = 10**-1
kcatA = 10**-2
kcatB = 10**-2
Tcc = 1000
varTcc = 10

# intiate grid 
grid = gf.Grid(101,101,maxCells)

# seed grid 
params = [PprodA,kcatA,kcatB]
grid.seed('cascade',params,Tcc,varTcc)

# run grid 
grid.run()

with open(PROJECT_DIR / 'gridcells/mac_cascade/cascade_10gen2.pickle','wb') as f:
    pickle.dump(grid,f,pickle.HIGHEST_PROTOCOL)
