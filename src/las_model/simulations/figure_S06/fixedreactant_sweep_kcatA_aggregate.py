import os, pickle, numpy as np
from las_model.utils.config import PROJECT_DIR

nCells = 1000
Tcc = 1000
kcatA = 10**-1
PprodAs = np.logspace(-2,2,5)
PprodBs = np.logspace(-3,3,31)
Km = 10**3

means = np.zeros([len(PprodAs),len(PprodBs),3])
variances = np.zeros_like(means)
vardSis = np.zeros_like(means)
vardRnd = np.zeros_like(means)
normvars = np.zeros_like(means)

for file in os.listdir(PROJECT_DIR / 'fixed_reactant/fixedReactant5'):
    print('processing file %s' % file)
    kcatAindex = int(file.split('_')[2])
    prodBindex = int(file.split('_')[4].split('.')[0])

    with open(PROJECT_DIR / 'fixed_reactant/fixedReactant5/' + file,'rb') as f:
        molecules,volume,times,vardsis,vardrnd,normvar = pickle.load(f)

    means[kcatAindex,prodBindex] = np.mean(molecules/volume,axis=1)
    variances[kcatAindex,prodBindex] = np.var(molecules/volume,axis=1)
    vardSis[kcatAindex,prodBindex] = vardsis
    vardRnd[kcatAindex,prodBindex] = vardrnd
    normvars[kcatAindex,prodBindex] = normvar

with open(PROJECT_DIR / 'fixed_reactant/fixedReactant5.pickle','wb') as f:
    pickle.dump([means,variances,vardSis,vardRnd,normvars],f,pickle.HIGHEST_PROTOCOL)
