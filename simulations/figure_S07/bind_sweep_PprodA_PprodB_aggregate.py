import os
import pickle
import numpy as np
from utils.config import PROJECT_DIR

prodAs = np.zeros([31,31])
prodBs = np.zeros_like(prodAs)
Aeqs = np.zeros_like(prodAs)
Beqs = np.zeros_like(prodAs)
Ceqs = np.zeros_like(prodAs)
varAs = np.zeros_like(prodAs)
varBs = np.zeros_like(prodAs)
varCs = np.zeros_like(prodAs)

for file in os.listdir(PROJECT_DIR / 'prodRateSat/bind3'):
    print('processing file %s' % file)
    prodAindex = int(file.split('_')[2])
    prodBindex = int(file.split('_')[4].split('.')[0])

    with open(PROJECT_DIR / 'prodRateSat/bind3' / file,'rb') as f:
        prodA,prodB,Aeq,Beq,Ceq,normvarA,normvarB,normvarC = pickle.load(f)

    prodAs[prodAindex,prodBindex] = prodA
    prodBs[prodAindex,prodBindex] = prodB
    Aeqs[prodAindex,prodBindex] = Aeq
    Beqs[prodAindex,prodBindex] = Beq
    Ceqs[prodAindex,prodBindex] = Ceq
    varAs[prodAindex,prodBindex] = normvarA
    varBs[prodAindex,prodBindex] = normvarB
    varCs[prodAindex,prodBindex] = normvarC
    

with open(PROJECT_DIR / 'prodRateSat/bind3.pickle','wb') as f:
    pickle.dump([prodAs,prodBs,Aeqs,Beqs,Ceqs,varAs,varBs,varCs],f,pickle.HIGHEST_PROTOCOL)