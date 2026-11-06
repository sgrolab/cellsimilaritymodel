"""
Analytical lineage-associated similarity (LAS) of a linear production chain.

Species 0 is made at mean rate P, in bursts of burstSize molecules; species m >= 1 is made at
rate rates[m-1] per molecule of species m-1; species m is degraded at rate degradation[m]. The
cell cycle has a fixed length Tcc, and at division each molecule goes to either daughter with
probability 1/2. The moment equations of such a linear network close, so the means and
covariances at division follow exactly from one cycle of them plus the division map, and the
LAS of every species at the start of each generation follows from those.
"""
import numpy as np
from scipy.linalg import expm


def chain_matrix(rates, degradation=None):
    """The production rates on the subdiagonal, minus the degradation rates on the diagonal."""
    n = len(rates) + 1
    M = np.zeros((n,n))
    M[np.arange(1,n), np.arange(n-1)] = rates
    if degradation is not None:
        M -= np.diag(degradation)
    return M


def division_moments(P, rates, Tcc, burstSize=1, degradation=None):
    """
    Steady-state mean vector and covariance matrix of the molecule counts at division.

    Within a cycle
        d mu/dt    = M mu + P e_0
        d Sigma/dt = M Sigma + Sigma M^T + diag(q),
        q_0 = P*burstSize + g_0 mu_0,   q_m = k_m mu_{m-1} + g_m mu_m,
    which is linear in z = (mu, vec Sigma, 1), so a cycle is z -> expm(A Tcc) z. A daughter keeps
    mu/2 and Sigma/4 + diag(mu)/4, and the moments at birth are the fixed point of that map.
    Returns (mu_T, Sigma_T).
    """
    M = chain_matrix(rates, degradation)
    n = len(M)
    k = np.concatenate([[0], rates])
    g = np.zeros(n) if degradation is None else np.asarray(degradation, dtype=float)
    mu = np.arange(n)                                    # indices of z
    S = n + np.arange(n*n).reshape(n,n)
    c = n + n*n
    A = np.zeros((c+1,c+1))
    A[np.ix_(mu,mu)] = M
    A[mu[0],c] = P
    for i in range(n):
        for j in range(n):
            A[S[i,j],S[:,j]] += M[i,:]
            A[S[i,j],S[i,:]] += M[j,:]
        A[S[i,i],mu[i]] += g[i]
        if i > 0:
            A[S[i,i],mu[i-1]] += k[i]
    A[S[0,0],c] += P*burstSize

    R = np.zeros((c+1,c+1))                              # division
    R[mu,mu] = 1/2
    R[S.ravel(),S.ravel()] = 1/4
    R[S[mu,mu],mu] += 1/4
    R[c,c] = 1

    Phi = expm(A*Tcc)
    G = R @ Phi
    zb = np.append(np.linalg.solve(np.eye(c) - G[:c,:c], G[:c,c]), 1)
    zT = Phi @ zb
    Sigma = zT[S]
    return zT[mu], (Sigma + Sigma.T)/2


def las(mu, Sigma, M, Tcc, gens):
    """
    LAS of every species at the start of each generation g in gens (g need not be an integer),
    given the moments at division mu, Sigma and the chain matrix M. Sisters share
    (Sigma - diag mu)/4 at birth; the shared surplus feeds downstream through expm(M g Tcc), and
    each division quarters it. Returns an array (len(gens), len(mu)).
    """
    C0 = (Sigma - np.diag(mu))/4        # covariance between sisters at birth
    v = (np.diag(Sigma) + mu)/4         # variance of a newborn
    values = []
    for g in gens:
        E = expm(M*g*Tcc)
        values.append(np.diag(E @ C0 @ E.T)/v/4**g)
    return np.array(values)


def las_chain(mu, Sigma, rates, Tcc, gens):
    """las() for a chain without degradation, given its production rates."""
    return las(mu, Sigma, chain_matrix(rates), Tcc, gens)


def peak_and_duration(gens, values):
    """
    The peak of a LAS curve, and the generation after which it stays below half that peak
    (linearly interpolated between the points of gens).
    """
    peak = np.max(values)
    above = np.nonzero(values >= peak/2)[0][-1]
    if above == len(values) - 1:
        return peak, gens[-1]
    g0, g1, v0, v1 = gens[above], gens[above+1], values[above], values[above+1]
    return peak, g0 + (v0 - peak/2)*(g1 - g0)/(v0 - v1)
