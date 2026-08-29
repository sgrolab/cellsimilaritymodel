import os
import pickle
import numpy as np
from las_model.utils.config import PROJECT_DIR

prodAs = np.zeros([31,31])
prodBs = np.zeros([31,31])
Aeqs = np.zeros_like(prodBs)
Beqs = np.zeros_like(prodBs)
Ceqs = np.zeros_like(prodBs)
Deqs = np.zeros_like(prodBs)
varAs = np.zeros_like(prodBs)
varBs = np.zeros_like(prodBs)
varCs = np.zeros_like(prodBs)
varDs = np.zeros_like(prodBs)
normvarAs = np.zeros_like(prodBs)
normvarBs = np.zeros_like(prodBs)
normvarCs = np.zeros_like(prodBs)
normvarDs = np.zeros_like(prodBs)

for file in os.listdir(PROJECT_DIR / 'prodRateSat/phos_cycle2'):
    print('processing file %s' % file)
    prodAindex = int(file.split('_')[2])
    prodBindex = int(file.split('_')[4].split('.')[0])
    
    with open(PROJECT_DIR / 'prodRateSat/phos_cycle2/' / file,'rb') as f:
       params,Eqs,var,normvars = pickle.load(f)

    prodAs[prodAindex,prodBindex] = params[0]
    prodBs[prodAindex,prodBindex] = params[1]
    Aeqs[prodAindex,prodBindex] = Eqs[0]
    Beqs[prodAindex,prodBindex] = Eqs[1]
    Ceqs[prodAindex,prodBindex] = Eqs[2]
    Deqs[prodAindex,prodBindex] = Eqs[3]
    varAs[prodAindex,prodBindex] = var[0]
    varBs[prodAindex,prodBindex] = var[1]
    varCs[prodAindex,prodBindex] = var[2]
    varDs[prodAindex,prodBindex] = var[3]
    normvarAs[prodAindex,prodBindex] = normvars[0]
    normvarBs[prodAindex,prodBindex] = normvars[1]
    normvarCs[prodAindex,prodBindex] = normvars[2]
    normvarDs[prodAindex,prodBindex] = normvars[3]

with open(PROJECT_DIR / 'phos_cycle/phoscycle2.pickle','wb') as f:
    pickle.dump([prodAs,prodBs,Aeqs,Beqs,Ceqs,Deqs,varAs,varBs,varCs,varDs,normvarAs,normvarBs,normvarCs,normvarDs],f,pickle.HIGHEST_PROTOCOL)