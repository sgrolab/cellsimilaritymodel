import numpy as np 
from las_model.utils import motiffunc as mf 

def calculate_division_differences(divStates, rng):
    """
    Calculate division differences for a set of cell division states.

    Parameters:
    -----------
    divStates : np.ndarray
        Array of division states with shape (nVars, nCells).
    rng : np.random.Generator
        Random number generator for stochastic simulations.

    Returns:
    --------
    dsis : np.ndarray
        Array of division state differences for individual cells.
    drnd : np.ndarray
        Array of division state differences for random cell pairs.
    vardsis : np.ndarray
        Variance of dsis across cells.
    vardrnd : np.ndarray
        Variance of drnd across cells.
    normvar : np.ndarray
        Normalized variance (1 - vardsis / vardrnd).
    """

    nCells = divStates.shape[1]
    divStatesInt = divStates.astype('int')

    # One binomial draw per cell
    cell1 = rng.binomial(divStatesInt, 0.5)

    # One random partner index per cell, then gather + draw for all partners at once
    partnerIdx = rng.integers(0, nCells, size=nCells)
    cell2 = rng.binomial(divStatesInt[:, partnerIdx], 0.5)

    dsis = divStates - 2 * cell1
    drnd = cell1 - cell2

    vardsis = np.var(dsis, axis=1)
    vardrnd = np.var(drnd, axis=1)
    normvar = 1 - (vardsis / vardrnd)

    return dsis, drnd, vardsis, vardrnd, normvar

def calculate_offspring_similarity_time(motherCell,metadata,rng):
    # Get division states and create offspring cells 
    divStates = (motherCell.getMotherStates()).astype('int')

    sis1states = rng.binomial(divStates, 0.5)
    sis2states = divStates - sis1states 

    partnerIdx = rng.integers(0, metadata['nCells'], size=metadata['nCells'])
    rnd1states = rng.binomial(divStates[:, partnerIdx], 0.5)

    # preallocate molecules lists
    molcules = {
        'sis1': [],
        'sis2': [],
        'rnd1': []
    }

    # Divide cells and run offspring
    for i in range(metadata['nCells']):
        print(f"Simulating cell {i+1}/{metadata['nCells']}")
        
        sis1 = mf.Cell(metadata['Tcc'],metadata['varTcc'],rng)
        sis1.inherit(motherCell,sis1states[:,i])
        sis1.run(metadata['nCycles'])
        molecules_sis1 = sis1.getMolecules()
        
        sis2 = mf.Cell(metadata['Tcc'],metadata['varTcc'],rng)
        sis2.inherit(motherCell,sis2states[:,i])
        sis2.run(metadata['nCycles'])
        molecules_sis2 = sis2.getMolecules()
        
        rnd1 = mf.Cell(metadata['Tcc'],metadata['varTcc'],rng)
        rnd1.inherit(motherCell,rnd1states[:,i])
        rnd1.run(metadata['nCycles'])
        molecules_rnd1 = rnd1.getMolecules()

        molcules['sis1'].append(molecules_sis1)
        molcules['sis2'].append(molecules_sis2)
        molcules['rnd1'].append(molecules_rnd1)

    # Stack molecule lists and compute pairwise differences for this kcat 
    sis1stack = np.stack(molcules['sis1'],axis=1)
    sis2stack = np.stack(molcules['sis2'],axis=1)
    rnd1stack = np.stack(molcules['rnd1'],axis=1)

    dsis = sis1stack - sis2stack
    drnd = sis1stack - rnd1stack

    vardsis = np.var(dsis,axis=1)
    vardrnd = np.var(drnd,axis=1)
    normvar = 1-vardsis/vardrnd

    return dsis, drnd, vardsis, vardrnd, normvar