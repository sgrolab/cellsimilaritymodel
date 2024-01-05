# Cells on a grid

import os, numpy as np 
from matplotlib import pyplot as plt 
os.chdir('//dm11.hhmi.org/sgrolab/mark/comp_proj/gridcells')
import gridfunc as gf
import importlib
importlib.reload(gf)
import pickle 

maxCells = 10


# intiate grid 
grid = gf.Grid(71,71,maxCells)

# seed grid 
params = [10**-2,10**-2]
grid.seed('prodsat',params)

# run grid 
grid.run()

#%%

# grid.makeVideo(10,10)
grid.cousinVideo(0,10,10)

# relatedness = grid.calcCollectiveLocalRelatedness(8, 10)



# grid.Mvideo(10,10)
# grid.Avideo(10,10)

#%%

for i in range(len(grid.Cells)):
    plt.plot(grid.Cells[i].t,grid.Cells[i].C/grid.Cells[i].V)

#%%

grid.Avideo(10,10)

#%% 

with open('grid_1000.pickle','wb') as f:
    pickle.dump(grid,f,pickle.HIGHEST_PROTOCOL)

