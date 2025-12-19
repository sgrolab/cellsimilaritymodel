# Fixed Reactant 5: Sweep kcatA and PprodB 
import sys
import pickle
import numpy as np 
from utils.config import PROJECT_DIR

# fixed parameters
nCells = 1000
Tcc = 1000

PprodA = 10**-1
Km = 10**3

# swept parameter 
kcatAs = np.logspace(-2,2,5)
kcatAindex = int(sys.argv[1])
kcatA = kcatAs[kcatAindex]
PprodBs = np.logspace(-3,3,31)
PprodBindex = int(sys.argv[2])
PprodB = PprodBs[PprodBindex]

# storage arrays 
molecules = np.zeros([3,int(1e9)])
volume = np.zeros(int(1e9))
times = np.zeros_like(volume)

enzyme_i = PprodA * Tcc
reactant_i = PprodB * Tcc
product_i = 3/2 * kcatA/2*(Km+enzyme_i+reactant_i-np.sqrt((Km+enzyme_i+reactant_i)**2-4*enzyme_i*reactant_i))

molecules[:,0] = [enzyme_i,reactant_i,product_i]
volume[0] = 1

# initialize counters
gen = 0
n = 1

# initialize random generator 
rng = np.random.default_rng(seed=1000)

while gen < nCells:
    
    while volume[n-1] < 2:
        
        # get enzyme amount from previous timestep 
        A = molecules[0,n-1]
        
        # set amount of substrate based on volume 
        B = volume[n-1] * PprodB * Tcc
        
        # calculate reaction probabilities 
        prodA = PprodA
        prodB = kcatA/2 * (Km+A+B-np.sqrt((Km+A+B)**2-4*A*B))
        
        Rtot = prodA + prodB
        
        # generate random numbers
        r1 = rng.uniform()
        r2 = rng.uniform() * Rtot
        
        # calculate time step
        tau = 1/Rtot*np.log(1/r1)
        
        # pick reaction 
        if r2 < prodA: # produce enzyme 
            molecules[0,n] = molecules[0,n-1] + 1
            molecules[2,n] = molecules[2,n-1]
        else:           # produce product 
            molecules[0,n] = molecules[0,n-1]
            molecules[2,n] = molecules[2,n-1] + 1
    
        # update volume 
        volume[n] = volume[n-1] + tau/Tcc
        
        # update time
        times[n] = times[n-1] + tau
        
        # update amount of substrate 
        molecules[1,n] = volume[n] * PprodB * Tcc
        
        # update counter 
        n+=1
    
    # divide 
    volume[n-1] = 1
    molecules[0,n-1] = rng.binomial(molecules[0,n-1],0.5)
    molecules[1,n-1] = volume[n-1] * PprodB*Tcc
    molecules[2,n-1] = rng.binomial(molecules[2,n-1],0.5)
    
    # update generation counter 
    gen += 1

# trim zeros 
molecules = molecules[:,0:n]
volume = volume[0:n]
times = times[0:n]

# calculate molecule concentration statistics 
means = np.mean(molecules/volume,axis=1)
variances = np.var(molecules/volume,axis=1)

# get mother states 
motherTimes = np.where(volume==1)[0][1::]-1
motherMolecules = molecules[:,motherTimes]
motherMolecules[1] = np.round(motherMolecules[1])

dsis = np.zeros([nCells,3])
drnd = np.zeros_like(dsis)

for i in range(nCells):
    cell1 = rng.binomial(motherMolecules[:,i].astype('int'),0.5)
    cell2 = rng.binomial(motherMolecules[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    dsis[i] = motherMolecules[:,i] - 2*cell1
    drnd[i] = cell1 - cell2

vardsis = np.var(dsis,axis=0)
vardrnd = np.var(drnd,axis=0)
normvar = 1-vardsis/vardrnd

with open(PROJECT_DIR / 'fixed_reactant/fixedReactant5/fixedreactant5_kcatA_%.2i_PprodB_%.2i.pickle' % (kcatAindex,PprodBindex),'wb') as f:
    pickle.dump([molecules,volume,times,vardsis,vardrnd,normvar],f,pickle.HIGHEST_PROTOCOL)