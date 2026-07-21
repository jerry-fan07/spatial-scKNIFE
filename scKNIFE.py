import time

s = time.perf_counter()

import numpy as np
from numpy.random import default_rng
#import matplotlib.pyplot as plt
#import seaborn as sns

import scipy.sparse as sp
import scipy.sparse.linalg as spla
#from scKNIFE_graph.src.scKNIFE_graph.config import PRO_DIR
from matrix_data import get_matrices

print(f"{time.perf_counter() - s : .3f}", "import times")

# R_g maps genes to nodes in unified graph

def NMF(X: sp.csr_array, r: int, iters: int, 
        epsilon: float=1e-9, seed: tuple[int, int]=(10,11)):
    # generate matrices of 1's with dimension N*r, r*p
    #W = np.ones((len(X), r), dtype=float)
    #H = np.ones((r, len(X[0])), dtype=float)
    N, M = X.shape
    W = default_rng(seed[0]).random((N, r), dtype=float)
    H = default_rng(seed[1]).random((r, M), dtype=float)

    coo = X.tocoo()
    rows, cols, x_data = coo.row, coo.col, coo.data

    def build_Q():
        loop = time.perf_counter()

        wh = np.zeros(len(rows), dtype=float)
        for i in range(W.shape[1]):
            wh += W[rows, i] * H[i, cols]

        print(time.perf_counter() - loop, "loop")

        
        Q = sp.csr_matrix((x_data / (wh + epsilon), (rows, cols)), shape=(N,M))
        return Q

    for t in range(iters):
        Q = build_Q()
        # H sum broadcasted over rest of the rows
        W *= Q.dot(H.T) / (H.sum(axis=1)[np.newaxis, :] + epsilon)

        Q = build_Q()
        # must be sparse.dot(dense)
        H *= (Q.T.dot(W)).T / (W.sum(axis=0)[:, np.newaxis] + epsilon)

    return W, H

# W, H dense arrays
def scKNIFE(X: sp.csr_array, R_g: np.ndarray, 
            D: np.ndarray, A: np.ndarray,
            r: int, lam: float, beta_W: float, beta_H: float, 
            NMF_iters: int=1, scKNIFE_iters: int=1):
    
    L = D - A
    interp = True

    start = time.perf_counter()

    if interp:
        #RgT_X = R_g.T @ X
        #RgT_Rg = R_g.T @ R_g

        #DONE: implement the base NMF algorithm to 
        # determine the values of labeled nodes in W and H
        # 07/06/2026: 
        # read ESL book pages, derived NMF algorithm in P-14.23, 
        # proof from monotone convergence
        # implemented NMF

        W_f, H = NMF(X, r, NMF_iters)
        
        # 07/07/2026:
        # (L_uu + \epsilon I) W_u = - L_uf W_f
        # minimizes the Laplacian quadratic form from
        # block matrix [W_u W_f]^T
        #DONE: implement the interpolation

        i_u = L.shape[0] - len(W_f)
        # L_uu + I is symmetric positive definite; give splu a float CSC
        # matrix and the symmetric fill-reducing ordering (MMD_AT_PLUS_A)
        # instead of the default COLAMD -> ~18x fewer nonzeros in the LU.
        L = L.tocsc().astype(float)
        lhs = L[:i_u, :i_u] + 1e-6 * sp.identity(i_u, format='csc')

        invA = spla.splu(lhs, permc_spec='MMD_AT_PLUS_A')
        W_u = invA.solve(-L[:i_u, i_u:] @ W_f)
        W = np.vstack((W_u, W_f))

    print(f"{time.perf_counter() - start: .3f}", "interpolation")

    #DONE: implement the optimization 
    # with graph and regularization constraints
    for _ in range(scKNIFE_iters):
        pos = R_g.T @ (X @ H.T) + A @ W * lam
        neg = R_g.T @ (R_g @ W @ (H @ H.T)) + lam * (D @ W) + beta_W + 1e-9
        W *= (pos / neg)
        
        R_gW = R_g @ W
        pos = (X.T @ R_gW).T
        neg = (R_gW.T @ R_gW) @ H + beta_H + 1e-9
        H *= (pos / neg)
        print(_, "W.max", float(np.max(W)), "H.max", float(np.max(H)))

    #DONE: implement the graph construction
    #TODO: think about integration of spatial structure

    return W, H


start = time.perf_counter()
R_g, D, A, X = get_matrices()
print(f"{time.perf_counter() - start: .3f}", "get matrices")

# START TIMER
start = time.perf_counter()

#NMF(np.ones((100, 100), dtype=int), 10, 200)

#W, H = scKNIFE(np.ones((100,100)), R_g, D, A, 1, 0.1, 0.05, 0.05)

W, H = scKNIFE(X, R_g, D, A, 100, 0.1, 0.05, 0.05, 20, 20)

# END TIMER
print(f"{time.perf_counter() - start : .3f}", "scKNIFE")

np.save("W_matrix.npy", W)
np.save("H_matrix.npy", H)