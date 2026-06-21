# Analyze Spatial simluation data 
import os
import pickle
import numpy as np
from las_model.utils.config import PROJECT_DIR

# import grid data 
with open(PROJECT_DIR / 'gridcells/mac_cascade/cascade_10gen2.pickle','rb') as f:
    grid = pickle.load(f)

# compute relatedness curve 
maxRadius = 9
relatedness = grid.calcCollectiveLocalRelatedness(maxRadius, 10000)
with open(PROJECT_DIR / 'gridcells/mac_cascade/cascade_10gen2_relatedness.pickle','wb') as f:
    pickle.dump(relatedness,f,pickle.HIGHEST_PROTOCOL)

# pull cousin maps 
crop = [25,75]
ts = [2000,4000,6000,8000,10000]
imgs = []
for i in range(len(ts)):
    imgs.append(8-grid.cousinMap(100,ts[i])[crop[0]:crop[1],crop[0]:crop[1]])
    imgs[i][np.where(imgs[i]==10)] = 'NaN'
with open(PROJECT_DIR / 'gridcells/mac_cascade/cascade_10gen2_cousinmaps.pickle','wb') as f:
    pickle.dump(imgs,f,pickle.HIGHEST_PROTOCOL)

cousinNums = [100,200,300,400,500]
ts = [0,2000,4000,6000,8000,10000]
imgs = []
for j in range(len(cousinNums)):
    onecousinimgs = []
    for i in range(len(ts)):
        onecousinimgs.append(8-grid.cousinMap(cousinNums[j],ts[i])[crop[0]:crop[1],crop[0]:crop[1]])
        onecousinimgs[i][np.where(onecousinimgs[i]==10)] = 'NaN'
    imgs.append(onecousinimgs)
with open(PROJECT_DIR / 'gridcells/mac_cascade/cascade_10gen2_cousinmaps_all.pickle','wb') as f:
    pickle.dump(imgs,f,pickle.HIGHEST_PROTOCOL)

# get molecular concentration maps 
molAImgs = []
molBImgs = []
molCImgs = []
molTimes = [4000,6000,8000,10000]
for i in range(len(molTimes)):
    molAImgs.append(grid.getFrame(molTimes[i],'A')[crop[0]:crop[1],crop[0]:crop[1]])
    molBImgs.append(grid.getFrame(molTimes[i],'B')[crop[0]:crop[1],crop[0]:crop[1]])
    molCImgs.append(grid.getFrame(molTimes[i],'C')[crop[0]:crop[1],crop[0]:crop[1]])
with open(PROJECT_DIR / 'gridcells/mac_cascade/cascade_10gen2_molConcMaps.pickle','wb') as f:
    pickle.dump([molAImgs,molBImgs,molCImgs],f,pickle.HIGHEST_PROTOCOL)

# get Moran Is
morIs_discdist = np.zeros([9,5,101])
filenames = []
for file in os.listdir(PROJECT_DIR / 'gridcells/mac_cascade/cascade_10gen2_moranIs'):
    if 'cascade_10gen2_morIs_discdist_r' in file:
        filenames.append(file)

for i in range(1,len(filenames)):
    with open(PROJECT_DIR / 'gridcells/mac_cascade/cascade_10gen2_moranIs' / filenames[i],'rb') as f:
        morIs_discdist[i-1] = pickle.load(f)
        
with open(PROJECT_DIR / 'gridcells/mac_cascade/cascade_10gen_moranIs2.pickle','wb') as f:
    pickle.dump(morIs_discdist,f,pickle.HIGHEST_PROTOCOL)
