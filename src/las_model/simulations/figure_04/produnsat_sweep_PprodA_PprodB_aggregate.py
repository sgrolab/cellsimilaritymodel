# Aggregate results from PprodA and PprodB sweep simulations 
import os 
import pickle
import numpy as np 
from las_model.utils.config import PROJECT_DIR

PprodAs = np.logspace(-2,1,16)
PprodBs = np.logspace(-2,4,31)

Ameans = np.zeros([len(PprodAs),len(PprodBs)])
Avars = np.zeros_like(Ameans)
Bmeans = np.zeros_like(Ameans)
Bvars = np.zeros_like(Ameans)
Cmeans = np.zeros_like(Ameans)
Cvars = np.zeros_like(Ameans)
motherMols = np.zeros([len(Ameans),3,1000])
dsiss = np.zeros([len(Ameans),1000,6])
drnds = np.zeros_like(dsiss)
vardsiss = np.zeros([len(PprodAs),len(PprodBs),3])
vardrnds = np.zeros_like(vardsiss)
normvarAs = np.zeros_like(Ameans)
normvarBs = np.zeros_like(Ameans)
normvarCs = np.zeros_like(Ameans)

for file in os.listdir(PROJECT_DIR / 'production/prodAprodBsweep1'):
    print('processing file %s' % file)
    prodAindex = int(file.split('_')[2])
    prodBindex = int(file.split('_')[4].split('.')[0])
    
    with open(PROJECT_DIR / 'production/prodAprodBsweep1' / file,'rb') as f:
        means,variances,divStates,dsis,drnd,vardrnd,vardsis,normvar = pickle.load(f)

    Ameans[prodAindex,prodBindex] = means[0]
    Avars[prodAindex,prodBindex] = variances[0]
    Bmeans[prodAindex,prodBindex] = means[1]
    Bvars[prodAindex,prodBindex] = variances[1]
    Cmeans[prodAindex,prodBindex] = means[2]
    Cvars[prodAindex,prodBindex] = variances[2]
    vardsiss[prodAindex,prodBindex] = vardsis[0:3]
    vardrnds[prodAindex,prodBindex] = vardrnd[0:3]
    normvarAs[prodAindex,prodBindex] = normvar[0]
    normvarBs[prodAindex,prodBindex] = normvar[1]
    normvarCs[prodAindex,prodBindex] = normvar[2]

with open(PROJECT_DIR / 'production/production_prodAprodBsweep1.pickle','wb') as f:
    pickle.dump([Ameans,Avars,Bmeans,Bvars,Cmeans,Cvars,normvarAs,normvarBs,normvarCs],f,pickle.HIGHEST_PROTOCOL)
