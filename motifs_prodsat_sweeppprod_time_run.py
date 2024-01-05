# Binding motif initial values 

import os, pickle, numpy as np 
os.chdir('/groups/sgro/sgrolab/mark/comp_proj/motifs')
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/motifs')
import motiffunc as mf

Tcc = 1000
nCells = 1000
nCycles = 10
rng = np.random.default_rng(seed=1000)

partitions = ['binomial','perfect','correlated']
concentrations = np.zeros([3,3,nCells,5,int(nCycles*Tcc/10+1)])

for j in range(len(partitions)):
    motherCell = mf.Cell(Tcc,0)
    motherCell.parameterize('prodsat',[0.01,0.01])
    motherCell.run(nCells,partitions[j])
    divStates = motherCell.getMotherStates()
    
    for i in range(nCells):
        if i % 10 == 0:
            print('partition method: %i, cell: %i' % (j,i))
        
        if partitions[j] == 'binomial':
            sis1state = rng.binomial(divStates[:,i].astype('int'),0.5)
            sis2state = divStates[:,i] - sis1state
            rnd1state = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
        elif partitions[j] == 'perfect':
            sis1state = divStates[:,i] // 2
            sis2state = divStates[:,i] - sis1state
            rnd1state = divStates[:,rng.integers(0,nCells).astype('int')] // 2
        elif partitions[j] == 'correlated':
            coef = rng.binomial(10000,0.5)/10000
            coef2 = rng.binomial(10000,0.5)/10000
            sis1state = (divStates[:,i] * coef).astype('int')
            sis2state = divStates[:,i] - sis1state
            rnd1state = (divStates[:,rng.integers(0,nCells).astype('int')] * coef2).astype('int')
        
        sis1 = mf.Cell(Tcc,0)
        sis1.inherit(motherCell,sis1state)
        sis1.run(nCycles,partitions[j])
        concentrations[j,0,i] = sis1.getMolecules()
        
        sis2 = mf.Cell(Tcc,0)
        sis2.inherit(motherCell,sis2state)
        sis2.run(nCycles,partitions[j])
        concentrations[j,1,i] = sis2.getMolecules()
        
        rnd1 = mf.Cell(Tcc,0)
        rnd1.inherit(motherCell,rnd1state)
        rnd1.run(nCycles,partitions[j])
        concentrations[j,2,i] = rnd1.getMolecules()

# export grid 
# os.chdir('/groups/sgro/sgrolab/mark/comp_proj/durationSweep')
os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/durationSweep')
with open('partitionsweep.pickle','wb') as f:
    pickle.dump(concentrations,f,pickle.HIGHEST_PROTOCOL)


# from matplotlib import pyplot as plt 

# plt.figure(figsize=(12,7))
# plt.subplot(3,3,1)
# plt.plot(motherCell_bin.t/motherCell_bin.Tcc,motherCell_bin.V)
# plt.subplot(3,3,2)
# plt.plot(motherCell_bin.t/motherCell_bin.Tcc,motherCell_bin.A/motherCell_bin.V)
# plt.subplot(3,3,3)
# plt.plot(motherCell_bin.t/motherCell_bin.Tcc,motherCell_bin.B/motherCell_bin.V)

# plt.subplot(3,3,4)
# plt.plot(motherCell_per.t/motherCell_per.Tcc,motherCell_per.V)
# plt.subplot(3,3,5)
# plt.plot(motherCell_per.t/motherCell_per.Tcc,motherCell_per.A/motherCell_per.V)
# plt.subplot(3,3,6)
# plt.plot(motherCell_per.t/motherCell_per.Tcc,motherCell_per.B/motherCell_per.V)

# plt.subplot(3,3,7)
# plt.plot(motherCell_cor.t/motherCell_cor.Tcc,motherCell_cor.V)
# plt.subplot(3,3,8)
# plt.plot(motherCell_cor.t/motherCell_cor.Tcc,motherCell_cor.A/motherCell_cor.V)
# plt.subplot(3,3,9)
# plt.plot(motherCell_cor.t/motherCell_cor.Tcc,motherCell_cor.B/motherCell_cor.V)


# plt.tight_layout()
