# Binding motif initial values 

import os, pickle, time, numpy as np 
# os.chdir('/groups/sgro/sgrolab/mark/comp_proj/motifs')
os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/motifs')
import motiffunc as mf

Tcc = 1000
nCells = 10
rng = np.random.default_rng(seed=1000)

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('prodsat',[0.01,1])

t0 = time.time()
motherCell.run(nCells)
t1 = time.time()

print(t1-t0)

#%% plot 

from matplotlib import pyplot as plt 


plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.plot(motherCell.t/motherCell.Tcc,motherCell.A/motherCell.V)
plt.subplot(1,2,2)
plt.plot(motherCell.t/motherCell.Tcc,motherCell.B/motherCell.V)


#%% export 

os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/motifs')
with open('samplerun.pickle','wb') as f:
    pickle.dump([mothers,divStates],f,pickle.HIGHEST_PROTOCOL)