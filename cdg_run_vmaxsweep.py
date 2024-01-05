# kcat sweep 

import pickle, os, sys, numpy as np

def perfectDivide(V,A,B,C,D,E):
    return [V//2,A//2,B//2,C//2,D//2,E//2]

def binomialDivide(V,A,B,C,D,E):
    return [1,rng.binomial(A,0.5),rng.binomial(B,0.5),rng.binomial(C,0.5),rng.binomial(D,0.5),rng.binomial(E,0.5)]

def corBinomialDivide(V,A,B,C,D,E):
    divProd = rng.binomial(A, .5) / A
    return [1,round(A*divProd),round(B*divProd),round(C*divProd),round(D*divProd),round(E*divProd)]

def stochCycle(y0):
    
    t_array = np.zeros(int(1e8))
    V_array = np.zeros_like(t_array)
    A_array = np.zeros_like(t_array)
    B_array = np.zeros_like(t_array)
    C_array = np.zeros_like(t_array)
    D_array = np.zeros_like(t_array)
    E_array = np.zeros_like(t_array)
    
    V_array[0] = y0[0]
    A_array[0] = y0[1]
    B_array[0] = y0[2]
    C_array[0] = y0[3]
    D_array[0] = y0[4]
    E_array[0] = y0[5]
    
    n = 1
    while t_array[n-1] < cycleTime:
        
        # update arrays 
        V_array[n] = V_array[n-1]
        A_array[n] = A_array[n-1]
        B_array[n] = B_array[n-1]
        C_array[n] = C_array[n-1]
        D_array[n] = D_array[n-1]
        E_array[n] = E_array[n-1]
        
        # get values 
        DGC = A_array[n]
        PDE = B_array[n]
        cdg = C_array[n]
        
        
        # calculate probabilities
        prodA = prodDGC
        prodB = prodPDE
        prod_cdg = k_dgc * DGC
        deg_cdg = k_pde*PDE*cdg/(kM_pde+cdg)
        prodP = k_gene * cdg
        
        Rtot = prodA + prodB + prod_cdg + deg_cdg + prodP
        
        # generate random numbers
        r1 = rng.uniform()
        r2 = rng.uniform() * Rtot
        
        # calculate time step
        tau = 1/Rtot*np.log(1/r1)
        
        # pick reaction 
        if r2 < prodA:
            A_array[n] = A_array[n] + 1
        elif r2 < prodA + prodB:
            B_array[n] = B_array[n] + 1
        elif r2 < prodA + prodB + prod_cdg:
            C_array[n] = C_array[n] + 1 
        elif r2 < prodA + prodB + prod_cdg + deg_cdg:
            C_array[n] = C_array[n] - 1
        else:
            D_array[n] = D_array[n] + 1
            
        # calculate cell growth 
        V_array[n] = V_array[n] + tau*growthRate
        
        # update time 
        t_array[n] = t_array[n-1] + tau
        
        # update counter
        n = n+1
        
    return [t_array[0:n],V_array[0:n],A_array[0:n],B_array[0:n],C_array[0:n],D_array[0:n],E_array[0:n]]

def stochModel(y0,nCycles,divProb='binomial'):
    
    sol_t = []
    sol_V = []
    sol_A = []
    sol_B = []
    sol_C = []
    sol_D = []
    sol_E = []
    
    t0 = 0
    
    # steady-state M and B 
    for i in range(0,nCycles):
        nextGenArrs = stochCycle(y0)
        if divProb == 'binomial':
            y0 = binomialDivide(nextGenArrs[1][-1],nextGenArrs[2][-1],nextGenArrs[3][-1],nextGenArrs[4][-1],nextGenArrs[5][-1],nextGenArrs[6][-1])
        elif divProb == 'corBinomial':
            y0 = corBinomialDivide(nextGenArrs[1][-1],nextGenArrs[2][-1],nextGenArrs[3][-1],nextGenArrs[4][-1],nextGenArrs[5][-1],nextGenArrs[6][-1])
        else:
            y0 = perfectDivide(nextGenArrs[1][-1],nextGenArrs[2][-1],nextGenArrs[3][-1],nextGenArrs[4][-1],nextGenArrs[5][-1],nextGenArrs[6][-1])
        
        sol_t = np.concatenate((sol_t,nextGenArrs[0]+t0))
        t0 = sol_t[-1]
        sol_V = np.concatenate((sol_V,nextGenArrs[1]))
        sol_A = np.concatenate((sol_A,nextGenArrs[2]))
        sol_B = np.concatenate((sol_B,nextGenArrs[3]))
        sol_C = np.concatenate((sol_C,nextGenArrs[4]))
        sol_D = np.concatenate((sol_D,nextGenArrs[5]))
        sol_E = np.concatenate((sol_E,nextGenArrs[6]))
    
    return [sol_t,sol_V,sol_A,sol_B,sol_C,sol_D,sol_E]

def createSisCells(ys,divProb='binomial'):
    if divProb=='binomial':
        cell1_A = rng.binomial(ys[1],0.5)
        cell1_B = rng.binomial(ys[2],0.5)
        cell1_C = rng.binomial(ys[3],0.5)
        cell1_D = rng.binomial(ys[4],0.5)
        cell1_E = rng.binomial(ys[5],0.5)
    elif divProb == 'corBinomial':
        divPrct = rng.binomial(ys[0],0.5) / ys[0]
        cell1_A = round(ys[0] * divPrct)
        cell1_B = round(ys[1] * divPrct)
        cell1_C = round(ys[2] * divPrct)
        cell1_D = round(ys[3] * divPrct)
        cell1_E = round(ys[4] * divPrct)
    else:
        cell1_A = ys[0] // 2
        cell1_B = ys[1] // 2
        cell1_C = ys[2] // 2
        cell1_D = ys[3] // 2
        cell1_E = ys[4] // 2
    
    cell2_A = ys[1] - cell1_A
    cell2_B = ys[2] - cell1_B 
    cell2_C = ys[3] - cell1_C
    cell2_D = ys[4] - cell1_D
    cell2_E = ys[5] - cell1_E
    
    return [1,cell1_A,cell1_B,cell1_C,cell1_D,cell1_E],[1,cell2_A,cell2_B,cell2_C,cell2_D,cell2_E]

def reducedIndices(t,nCycles):
    times = np.linspace(0,cycleTime*nCycles,cycleTime*nCycles)
    red_indices = np.zeros(len(times))
    for i in range(len(red_indices)):
        red_indices[i] = np.argmin(abs(t - times[i]))
    return red_indices.astype('int')

# generate initial conditions for sister cells
# model parameters
cycleTime = 1000

# get value for kcatM
kcatindex = int(sys.argv[1])

prodDGC = 0.01
prodPDE = 0.01

k_dgc = np.logspace(np.log10(0.05),0,20)[kcatindex]
k_pde = 2*k_dgc-0.1

k_gene = 0.002

nCells= 1000
nCycles = 10
rng = np.random.default_rng(seed=1000)

growthRate = 1/cycleTime
Aeq = prodDGC * cycleTime
Beq = prodPDE * cycleTime
Ceq = (k_dgc*prodDGC*cycleTime-1/2*k_pde*prodPDE*cycleTime)*cycleTime
Deq = 3/2*Ceq*k_gene*cycleTime
Eeq = 0

kM_pde = Ceq

y0i = [1,Aeq,Beq,Ceq,Deq,Eeq]

# generate seed cells
seedcell = stochModel(y0i,nCells)

divindices = np.where(seedcell[1]>2)
predivstates = np.asarray((seedcell[1][divindices],seedcell[2][divindices],seedcell[3][divindices],seedcell[4][divindices],seedcell[5][divindices],seedcell[6][divindices]))

diffs_abs = np.zeros([nCells,2,5,cycleTime*nCycles])
diffs_conc = np.zeros_like(diffs_abs)

for i in range(nCells):
    # generate daughter cells 
    sis1_t0,sis2_t0 = createSisCells(predivstates[:,i])
    rnd1_t0,rnd2_t0 = createSisCells(predivstates[:,rng.integers(0,nCells)])
    
    # run daughter cells for 10 cycles 
    sis1 = stochModel(sis1_t0,nCycles)
    sis2 = stochModel(sis2_t0,nCycles)
    rnd1 = stochModel(rnd1_t0,nCycles)
    
    # get reduced indices
    indices_sis1 = reducedIndices(sis1[0],nCycles)
    indices_sis2 = reducedIndices(sis2[0],nCycles)
    indices_rnd1 = reducedIndices(rnd1[0],nCycles)
    
    # calculate deltas betweeen sister cells
    for j in range(2,7):
        diffs_abs[i,0,j-2] = sis1[j][indices_sis1] - sis2[j][indices_sis2]
        diffs_abs[i,1,j-2] = sis1[j][indices_sis1] - rnd1[j][indices_rnd1]
        diffs_conc[i,0,j-2] = sis1[j][indices_sis1]/sis1[1][indices_sis1] - sis2[j][indices_sis2]/sis2[1][indices_sis2]
        diffs_conc[i,1,j-2] = sis1[j][indices_sis1]/sis1[1][indices_sis1] - rnd1[j][indices_rnd1]/rnd1[1][indices_rnd1]
    

params = [cycleTime,prodDGC,prodPDE,k_dgc,k_pde,kM_pde,k_gene]

os.chdir('/groups/sgro/sgrolab/mark/comp_proj/cdg')
with open('cdg_vmaxsweep_' + str(kcatindex) + '.pickle','wb') as f:
    pickle.dump([params,predivstates,diffs_abs,diffs_conc],f,pickle.HIGHEST_PROTOCOL)
