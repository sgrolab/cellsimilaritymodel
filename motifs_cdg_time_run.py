# Binding motif initial values 

import os, pickle, numpy as np 
os.chdir('/groups/sgro/sgrolab/mark/comp_proj/motifs')
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/motifs')
import motiffunc as mf

Tcc = 1000
nCells = 1000
nCycles = 10
rng = np.random.default_rng(seed=1000)

prodA = 0.01
prodB = 0.01

kcatA = 0.1
kcatB = 0.1
KM = (kcatA * prodA - kcatB * prodB /2) * Tcc**2
kprod = 0.01

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('cdg',[prodA,prodB,kcatA,kcatB,KM,kprod])
motherCell.run(nCells)

divStates = motherCell.getMotherStates()

concentrations = np.zeros([3,nCells,5,int(nCycles*Tcc/10+1)])

for i in range(nCells):
    
    sis1state = rng.binomial(divStates[:,i].astype('int'),0.5)
    sis2state = divStates[:,i] - sis1state
    rnd1state = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    sis1 = mf.Cell(Tcc,0)
    sis1.inherit(motherCell,sis1state)
    sis1.run(nCycles)
    concentrations[0,i] = sis1.getMolecules()
    
    sis2 = mf.Cell(Tcc,0)
    sis2.inherit(motherCell,sis2state)
    sis2.run(nCycles)
    concentrations[1,i] = sis2.getMolecules()
    
    rnd1 = mf.Cell(Tcc,0)
    rnd1.inherit(motherCell,rnd1state)
    rnd1.run(nCycles)
    concentrations[2,i] = rnd1.getMolecules()

os.chdir('/groups/sgro/sgrolab/mark/comp_proj/cdg')
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/timeData')
with open('motifs_cdg.pickle','wb') as f:
    pickle.dump([motherCell,concentrations],f,pickle.HIGHEST_PROTOCOL)
