import os
import pickle
import numpy as np
from utils.config import PROJECT_DIR

PprodBs = np.logspace(-2,4,31)

substrateMean = np.zeros(len(os.listdir()))
substrateVar = np.zeros_like(substrateMean)
substrateLow = np.zeros_like(substrateMean)
substrateHigh = np.zeros_like(substrateMean)
substrateDrnd = np.zeros_like(substrateMean)
substrateDsis = np.zeros_like(substrateMean)
substrateSim = np.zeros_like(substrateMean)
orderMean = np.zeros_like(substrateMean)
orderStd = np.zeros_like(substrateMean)
enzymeMean = np.zeros_like(substrateMean)
enzymeVar = np.zeros_like(substrateMean)
enzymeDrnd = np.zeros_like(substrateMean)
enzymeDsis = np.zeros_like(substrateMean)
enzymeSim = np.zeros_like(substrateMean)
productMean = np.zeros_like(substrateMean)
productVar = np.zeros_like(substrateMean)
productDrnd = np.zeros_like(substrateMean)
productDsis = np.zeros_like(substrateMean)
productSim = np.zeros_like(substrateMean)

for file in os.listdir(PROJECT_DIR / 'orderAnalysis/calcOrder4'):
    
    prodBindex = int(file.split('_')[2].split('.')[0])
    
    with open(PROJECT_DIR / 'orderAnalysis/calcOrder4' / file,'rb') as f:
        mols,vards,order = pickle.load(f)

    substrateMean[prodBindex] = np.mean(mols[0])
    substrateVar[prodBindex] = np.var(mols[0])
    substrateLow[prodBindex] = substrateMean[prodBindex] - np.quantile(mols[0],0.1)
    substrateHigh[prodBindex] = np.quantile(mols[0],0.9) - substrateMean[prodBindex]
    substrateDrnd[prodBindex] = vards[0][0]
    substrateDsis[prodBindex] = vards[1][0]
    substrateSim[prodBindex] = vards[2][0]
    
    enzymeMean[prodBindex] = np.mean(mols[1])
    enzymeVar[prodBindex] = np.var(mols[1])
    enzymeDrnd[prodBindex] = vards[0][1]
    enzymeDsis[prodBindex] = vards[1][1]
    enzymeSim[prodBindex] = vards[2][1]
    
    productMean[prodBindex] = np.mean(mols[2])
    productVar[prodBindex] = np.var(mols[2])
    productDrnd[prodBindex] = vards[0][2]
    productDsis[prodBindex] = vards[1][2]
    productSim[prodBindex] = vards[2][2]
    
    orderMean[prodBindex] = np.nanmean(order)
    orderStd[prodBindex] = np.sqrt(np.nanvar(order))

with open(PROJECT_DIR / 'orderAnalysis/calcOrder4.pickle','wb') as f:
    pickle.dump([substrateMean,substrateVar,substrateLow,substrateHigh,substrateDrnd,substrateDsis,substrateSim,enzymeMean,enzymeVar,enzymeDrnd,enzymeDsis,enzymeSim,productMean,productVar,productDrnd,productDsis,productSim,orderMean,orderStd],f,pickle.HIGHEST_PROTOCOL)
