# Cascade Sweep Time: Aggregate Data  
import pickle
import numpy as np 
from utils.config import PROJECT_DIR

with open(PROJECT_DIR / 'cascade/motifs_cascade.pickle','rb') as f:
    motherCell,concentrations = pickle.load(f)

dsis = concentrations[0] - concentrations[1]
drnd = concentrations[0] - concentrations[2]

vardsis = np.var(dsis,axis=0)
vardrnd = np.var(drnd,axis=0)

normvar = 1-vardsis/vardrnd

with open(PROJECT_DIR / 'cascade/cascade_time_normdvar.pickle','wb') as f:
    pickle.dump(normvar,f,pickle.HIGHEST_PROTOCOL)

