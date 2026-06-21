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

for file in os.listdir(PROJECT_DIR / 'fixed_reactant/fixedReactant4'):
    print('processing file %s' % file)
    prodAindex = int(file.split('_')[2])
    prodBindex = int(file.split('_')[4].split('.')[0])

    with open(PROJECT_DIR / 'fixed_reactant/fixedReactant4/' + file,'rb') as f:
        molecules,volume,times,vardsis,vardrnd,normvar = pickle.load(f)

    means[prodAindex,prodBindex] = np.mean(molecules/volume,axis=1)
    variances[prodAindex,prodBindex] = np.var(molecules/volume,axis=1)
    vardSis[prodAindex,prodBindex] = vardsis
    vardRnd[prodAindex,prodBindex] = vardrnd
    normvars[prodAindex,prodBindex] = normvar

os.chdir('/groups/sgro/sgrolab/mark/comp_proj/fixed_reactant')
with open(PROJECT_DIR / 'fixed_reactant/fixedReactant4.pickle','wb') as f:
    pickle.dump([means,variances,vardSis,vardRnd,normvars],f,pickle.HIGHEST_PROTOCOL)




