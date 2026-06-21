# Consolidate data for production motif, PprodB sweep 
import os 
import pickle
import numpy as np 
from las_model.utils.config import PROJECT_DIR

PprodBs = np.logspace(-2,4,31)

Ameans = np.zeros(31)
Avars = np.zeros_like(Ameans)
Bmeans = np.zeros_like(Ameans)
Bvars = np.zeros_like(Ameans)
Cmeans = np.zeros_like(Ameans)
Cvars = np.zeros_like(Ameans)
vardsiss = np.zeros([len(Ameans),3])
vardrnds = np.zeros_like(vardsiss)
normvarAs = np.zeros_like(Ameans)
normvarBs = np.zeros_like(Ameans)
normvarCs = np.zeros_like(Ameans)

for file in os.listdir(PROJECT_DIR / 'production/prodBsweep5'):
    print('processing file %s' % file)
    prodBindex = int(file.split('_')[2].split('.')[0])
    
    filepath = os.path.join(PROJECT_DIR / 'production/prodBsweep5', file)

    with open(filepath,'rb') as f:
        Aeq,Beq,Ceq,Avar,Bvar,Cvar,vardrnd,vardsis,normvar = pickle.load(f)

    Ameans[prodBindex] = Aeq
    Avars[prodBindex] = Avar
    Bmeans[prodBindex] = Beq
    Bvars[prodBindex] = Bvar
    Cmeans[prodBindex] = Ceq
    Cvars[prodBindex] = Cvar
    vardsiss[prodBindex] = vardsis[0:3]
    vardrnds[prodBindex] = vardrnd[0:3]
    normvarAs[prodBindex] = normvar[0]
    normvarBs[prodBindex] = normvar[1]
    normvarCs[prodBindex] = normvar[2]

with open(PROJECT_DIR / 'production/production_prodBsweep5.pickle','wb') as f:
    pickle.dump([Ameans,Avars,Bmeans,Bvars,Cmeans,Cvars,normvarAs,normvarBs,normvarCs],f,pickle.HIGHEST_PROTOCOL)