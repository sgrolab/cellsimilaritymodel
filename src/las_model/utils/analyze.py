import numpy as np 

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