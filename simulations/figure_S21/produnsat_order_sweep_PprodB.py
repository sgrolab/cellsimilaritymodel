# Unsaturated Production: Order Analysis with varying PprodB 
import sys 
import pickle
import numpy as np 
from utils import motiffunc as mf
from utils.config import PROJECT_DIR

def calcProdRate(A,B,kcat,Km):
    return kcat/2*(A+B+Km-np.sqrt((A+B+Km)**2-4*A*B))

def calcOrder(B,kcat,Km,A):
    
    B1 = B+1 
    rate0 = kcat/2*(A+B+Km-np.sqrt((A+B+Km)**2-4*A*B))
    rate1 = kcat/2*(A+B1+Km-np.sqrt((A+B1+Km)**2-4*A*B1))
       
    logRate0 = np.log10(rate0)
    logRate1 = np.log10(rate1)
    
    logB = np.log10(B)
    logB1 = np.log10(B1)
    
    return (logRate1-logRate0)/(logB1-logB)

nCells = 1000
Tcc = 1000
rng = np.random.default_rng(seed=1000)

PprodA = 10**-1
kcatA = 10**-1
PprodBs = np.logspace(-2,4,31)
PprodBindex = int(sys.argv[1])
PprodB = PprodBs[PprodBindex]
Km = 10**3

motherCell = mf.Cell(Tcc,0)
motherCell.parameterize('produnsat',[PprodB,PprodA,kcatA,Km])
motherCell.equilibrate()
motherCell.run(nCells)

As = motherCell.A/motherCell.V
Bs = motherCell.B/motherCell.V
Cs = motherCell.C/motherCell.V

mols = [As,Bs,Cs]

divStates = motherCell.getMotherStates()
dsis = np.zeros([nCells,6])
drnd = np.zeros([nCells,6])

for k in range(nCells):
    cell1 = rng.binomial(divStates[:,k].astype('int'),0.5)
    cell2 = rng.binomial(divStates[:,rng.integers(0,nCells)].astype('int'),0.5)
    
    dsis[k] = divStates[:,k] - 2*cell1
    drnd[k] = cell1 - cell2

vardrnd = np.var(drnd,axis=0)
vardsis = np.var(dsis,axis=0)
normvar = 1-vardsis/vardrnd
vards = [vardrnd,vardsis,normvar]

# compute reaction order at each time step 
order = np.zeros(len(motherCell.A))
for j in range(len(order)):
    
    # get amount of A
    A = motherCell.B[j]
    B = motherCell.A[j]
    
    # get reaction curve 
    order[j] = calcOrder(B,kcatA,Km,A)
    

with open(PROJECT_DIR / 'orderAnalysis/calcOrder4/order4_PprodB_%.2i.pickle' % (PprodBindex),'wb') as f:
    pickle.dump([mols,vards,order],f,pickle.HIGHEST_PROTOCOL)