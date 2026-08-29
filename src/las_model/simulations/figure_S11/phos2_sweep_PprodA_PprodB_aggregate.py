import os
import pickle
import numpy as np
from las_model.utils.config import PROJECT_DIR

prodAs = np.zeros([31,31])
prodBs = np.zeros_like(prodAs)
Aeqs = np.zeros_like(prodAs)
Beqs = np.zeros_like(prodAs)
Ceqs = np.zeros_like(prodAs)
Deqs = np.zeros_like(prodAs)
Eeqs = np.zeros_like(prodAs)
Feqs = np.zeros_like(prodAs)
varAs = np.zeros_like(prodAs)
varBs = np.zeros_like(prodAs)
varCs = np.zeros_like(prodAs)
varDs = np.zeros_like(prodAs)
varEs = np.zeros_like(prodAs)
varFs = np.zeros_like(prodAs)

for file in os.listdir(PROJECT_DIR / 'prodRateSat/phos2_0'):
    print('processing file %s' % file)
    prodAindex = int(file.split('_')[3])
    prodBindex = int(file.split('_')[5].split('.')[0])

    with open(PROJECT_DIR / 'prodRateSat/phos2_0' / file,'rb') as f:
        prodA,prodB,Aeq,Beq,Ceq,Deq,Eeq,Feq,normvarA,normvarB,normvarC,normvarD,normvarE,normvarF = pickle.load(f)

    prodAs[prodAindex,prodBindex] = prodA
    prodBs[prodAindex,prodBindex] = prodB
    Aeqs[prodAindex,prodBindex] = Aeq
    Beqs[prodAindex,prodBindex] = Beq
    Ceqs[prodAindex,prodBindex] = Ceq
    Deqs[prodAindex,prodBindex] = Deq
    Eeqs[prodAindex,prodBindex] = Eeq
    Feqs[prodAindex,prodBindex] = Feq
    varAs[prodAindex,prodBindex] = normvarA
    varBs[prodAindex,prodBindex] = normvarB
    varCs[prodAindex,prodBindex] = normvarC
    varDs[prodAindex,prodBindex] = normvarD
    varEs[prodAindex,prodBindex] = normvarE
    varFs[prodAindex,prodBindex] = normvarF
    
with open(PROJECT_DIR / 'prodRateSat/phos2_0.pickle','wb') as f:
    pickle.dump([prodAs,prodBs,Aeqs,Beqs,Ceqs,Eeqs,Deqs,Feqs,varAs,varBs,varCs,varDs,varEs,varFs],f,pickle.HIGHEST_PROTOCOL)