# Binding motif initial values 

import os, pickle, numpy as np 
os.chdir('/groups/sgro/sgrolab/mark/comp_proj/motifs')
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/motifs')
import motiffunc as mf


Tcc = 1000
nCells = 100
rng = np.random.default_rng(seed=1000)

prodAs = np.logspace(-3,2,11)
prodBs = np.logspace(-2,0,3)

Aeqs = np.zeros([len(prodBs),len(prodAs)])
Beqs = np.zeros_like(Aeqs)
Ceqs = np.zeros_like(Aeqs)
normvarAs = np.zeros_like(Aeqs)
normvarBs = np.zeros_like(Aeqs)
normvarCs = np.zeros_like(Aeqs)

for j in range(len(prodBs)):
    for i in range(len(prodAs)):
        prodA = prodAs[i]
        prodB = prodBs[j]
        
        motherCell = mf.Cell(Tcc,0)
        motherCell.parameterize('bind',[prodA,prodB,0.001])
        motherCell.run(nCells)
        
        Aeqs[j,i] = np.mean(motherCell.A/motherCell.V)
        Beqs[j,i] = np.mean(motherCell.B/motherCell.V)
        Ceqs[j,i] = np.mean(motherCell.C/motherCell.V)
        
        divStates = motherCell.getMotherStates()

        dsis = np.zeros([nCells,5])
        drnd = np.zeros([nCells,5])

        for k in range(nCells):
            cell1 = rng.binomial(divStates[:,k].astype('int'),0.5)
            cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
            
            dsis[k] = divStates[:,k] - 2*cell1
            drnd[k] = cell1 - cell2
        
        normvarAs[j,i] = 1-np.var(dsis[:,0],axis=0)/np.var(drnd[:,0],axis=0)
        normvarBs[j,i] = 1-np.var(dsis[:,1],axis=0)/np.var(drnd[:,1],axis=0)
        normvarCs[j,i] = 1-np.var(dsis[:,2],axis=0)/np.var(drnd[:,2],axis=0)


os.chdir('/groups/sgro/sgrolab/mark/comp_proj/prodRateSat')
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/subsaturation')
with open('bind.pickle','wb') as f:
    pickle.dump([prodAs,prodBs,Aeqs,Beqs,Ceqs,normvarAs,normvarBs,normvarCs],f,pickle.HIGHEST_PROTOCOL)