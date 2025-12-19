# Calculate Moran Is 
import sys
import pickle
import numpy as np
from utils.config import PROJECT_DIR

filename = 'cascade_10gen2.pickle'
fileprefix = filename.split('.')[0]

with open(PROJECT_DIR / 'gridcells/mac_cascade/' + filename,'rb') as f:
    grid = pickle.load(f)

neighborhoodsize = int(sys.argv[1])
shape = str(sys.argv[2])

timepoints = range(0,int(grid.timepoints[-1]),100)

morIs = np.zeros([5,len(timepoints)])

for i in range(len(timepoints)):
    morIs[0,i] = grid.calcMoranI(neighborhoodsize,timepoints[i],'A',shape)
    morIs[1,i] = grid.calcMoranI(neighborhoodsize,timepoints[i],'B',shape)
    morIs[2,i] = grid.calcMoranI(neighborhoodsize,timepoints[i],'C',shape)

with open(PROJECT_DIR / 'gridcells/mac_cascade/cascade_10gen2_moranIs/' + fileprefix+ '_morIs_' + shape + '_r' + str(neighborhoodsize) + '.pickle','wb') as f:
    pickle.dump(morIs,f,pickle.HIGHEST_PROTOCOL)
