import os
import pickle
import numpy as np
from las_model.utils.config import PROJECT_DIR

prodAs = np.zeros([31,31])
prodBs = np.zeros_like(prodAs)
Aeqs = np.zeros_like(prodAs)
Beqs = np.zeros_like(prodAs)
Ceqs = np.zeros_like(prodAs)
varAs = np.zeros_like(prodAs)
varBs = np.zeros_like(prodAs)
varCs = np.zeros_like(prodAs)
normvarAs = np.zeros_like(prodAs)
normvarBs = np.zeros_like(prodAs)
normvarCs = np.zeros_like(prodAs)

for file in os.listdir(PROJECT_DIR / 'binding_rev/revbind4'):
    print('processing file %s' % file)
    prodAindex = int(file.split('_')[2])
    prodBindex = int(file.split('_')[4].split('.')[0])

    with open(PROJECT_DIR / 'binding_rev/revbind4' / file,'rb') as f:
        prodA,prodB,Aeq,Beq,Ceq,varA,varB,varC,normvarA,normvarB,normvarC = pickle.load(f)

    prodAs[prodAindex,prodBindex] = prodA
    prodBs[prodAindex,prodBindex] = prodB
    Aeqs[prodAindex,prodBindex] = Aeq
    Beqs[prodAindex,prodBindex] = Beq
    Ceqs[prodAindex,prodBindex] = Ceq
    varAs[prodAindex,prodBindex] = varA
    varBs[prodAindex,prodBindex] = varB
    varCs[prodAindex,prodBindex] = varC
    normvarAs[prodAindex,prodBindex] = normvarA
    normvarBs[prodAindex,prodBindex] = normvarB
    normvarCs[prodAindex,prodBindex] = normvarC
    
with open(PROJECT_DIR / 'binding_rev/revbind4.pickle','wb') as f:
    pickle.dump([prodAs,prodBs,Aeqs,Beqs,Ceqs,varAs,varBs,varCs,normvarAs,normvarBs,normvarCs],f,pickle.HIGHEST_PROTOCOL)
