# Binding motif initial values 

import os, sys, pickle, numpy as np 
os.chdir('/groups/sgro/sgrolab/mark/comp_proj/motifs')
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/motifs')
import motiffunc as mf

nCells = 1000
nCycles = 10
rng = np.random.default_rng(seed=1000)

Tccs = [500,1000,2000,5000,10000]
Tccindex = int(sys.argv[1])
# Tccindex = 3
Tcc = Tccs[Tccindex]

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('prodsat',[0.01,0.01])
motherCell.run(nCells)
divStates = motherCell.getMotherStates()

concentrations = np.zeros([3,nCells,5,int(nCycles*Tcc/10+1)])

for i in range(nCells):
    # if i % 10 == 0:
    #     print('%2.f %% done!' % (i/nCells*100))
    
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

# export grid 
os.chdir('/groups/sgro/sgrolab/mark/comp_proj/durationSweep')
# os.chdir('//prfs.hhmi.org/sgrolab/mark/comp_proj/durationSweep')
with open('Tccsweep_' + str(Tccindex) + '.pickle','wb') as f:
    pickle.dump(concentrations,f,pickle.HIGHEST_PROTOCOL)

