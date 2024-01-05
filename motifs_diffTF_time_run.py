# Binding motif initial values 

import os, pickle, numpy as np 
os.chdir('/groups/sgro/sgrolab/mark/comp_proj/motifs')
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/motifs')
import motiffunc as mf
# from matplotlib import pyplot as plt 

Tcc = 1000
nCells = 1000
nCycles = 10
rng = np.random.default_rng(seed=1000)

prodA = 0.1
prodB = 0.01
kbind = 0.001
kprod = 0.01

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('diffTF',[prodA,prodB,kbind,kprod])
motherCell.run(nCells)

# plt.figure(figsize=(16,4))
# plt.subplot(1,4,1)
# plt.plot(motherCell.t/motherCell.Tcc,motherCell.A/motherCell.V)
# plt.subplot(1,4,2)
# plt.plot(motherCell.t/motherCell.Tcc,motherCell.B/motherCell.V)
# plt.subplot(1,4,3)
# plt.plot(motherCell.t/motherCell.Tcc,motherCell.C/motherCell.V)
# plt.subplot(1,4,4)
# plt.plot(motherCell.t/motherCell.Tcc,motherCell.D/motherCell.V)

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

os.chdir('/groups/sgro/sgrolab/mark/comp_proj/diffTF')
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/timeData')
with open('motifs_diffTF.pickle','wb') as f:
    pickle.dump([motherCell,concentrations],f,pickle.HIGHEST_PROTOCOL)
