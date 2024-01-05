# Binding motif initial values 

import os, sys, pickle, numpy as np 
os.chdir('/groups/sgro/sgrolab/mark/comp_proj/motifs')
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/motifs')
import motiffunc as mf

Tcc = 1000
nCells = 1000
rng = np.random.default_rng(seed=1000)

prodAs = np.logspace(-3,3,7)
kBs = np.logspace(-3,3,7)

prodAindex = int(sys.argv[1])
prodA = prodAs[prodAindex]
kBindex = int(sys.argv[2])
kB = kBs[kBindex]

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('produnsat',[prodA,0.01,kB,1000])
motherCell.run(nCells)

divStates = motherCell.getMotherStates()

dsis = np.zeros([nCells,5])
drnd = np.zeros([nCells,5])

for i in range(nCells):
    cell1 = rng.binomial(divStates[:,i].astype('int'),0.5)
    cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    dsis[i] = divStates[:,i] - 2*cell1
    drnd[i] = cell1 - cell2

os.chdir('/groups/sgro/sgrolab/mark/comp_proj/subsaturation')
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/subsaturation')
with open('motifs_produnsat_prodA_' + str(prodAindex) + '_kB_' + str(kBindex) + '.pickle','wb') as f:
    pickle.dump([prodA,kB,motherCell,dsis,drnd],f,pickle.HIGHEST_PROTOCOL)


# #%%

# prodC = np.zeros(len(motherCell.t))
# for i in range(len(prodC)):
#     prodC[i] = motherCell.k1/2 * (motherCell.k2+motherCell.A[i]+motherCell.B[i]-np.sqrt((motherCell.k2+motherCell.A[i]+motherCell.B[i])**2-4*motherCell.A[i]*motherCell.B[i]))

# plt.plot(motherCell.t,prodC)

# #%%

# plt.figure(figsize=(16,16))
# for i in range(len(motherCells)):
#     plt.subplot(len(motherCells),4,i*4+1)
#     plt.plot(motherCells[i].t/motherCells[i].Tcc,motherCells[i].V)
#     plt.subplot(len(motherCells),4,i*4+2)
#     plt.plot(motherCells[i].t/motherCells[i].Tcc,motherCells[i].A/motherCells[i].V)
#     plt.subplot(len(motherCells),4,i*4+3)
#     plt.plot(motherCells[i].t/motherCells[i].Tcc,motherCells[i].B/motherCells[i].V)
#     plt.subplot(len(motherCells),4,i*4+4)
#     plt.plot(motherCells[i].t/motherCells[i].Tcc,motherCells[i].C/motherCells[i].V)
