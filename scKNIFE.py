import time

s = time.perf_counter()

import numpy as np
from numpy.random import default_rng
#import matplotlib.pyplot as plt
#import seaborn as sns

import scipy.sparse as sp
#from scKNIFE_graph.src.scKNIFE_graph.config import PRO_DIR
from matrix_data import get_matrices

print(f"{time.perf_counter() - s : .3f}", "import times")

# R_g maps genes to nodes in unified graph

def NMF(X: sp.csr_matrix, r: int, iters: int, 
        epsilon: float=1e-9, seed: tuple[int, int]=(10,11)):
    # generate matrices of 1's with dimension N*r, r*p
    #W = np.ones((len(X), r), dtype=float)
    #H = np.ones((r, len(X[0])), dtype=float)
    N, M = X.shape
    W = default_rng(seed[0]).random((N, r), dtype=float)
    H = default_rng(seed[1]).random((r, M), dtype=float)

    # for _ in range(iter):
    #     # update values of W and K:
    #     WH = W @ H
    #     for i in range(len(W)):
    #         for k in range(len(W[0])):
    #             pos, neg = 0, 0
    #             for j in range(X.shape[1]):
    #                 pos += H[k][j] * X[i, j] / (WH[i][j] + 1e-9)
    #                 neg += H[k][j] + epsilon
    #             W[i][k] = W[i][k] * pos / neg

    #     WH = W @ H
    #     for k in range(len(H)):
    #         for j in range(len(H[0])):
    #             pos, neg = 0, 0
    #             for i in range(len(X)):
    #                 pos += W[i][k] * X[i, j]   / (WH[i][j] + 1e-9)
    #                 neg += W[i][k] + 1e-9
    #             H[k][j] = H[k][j] * pos / neg 

    # can calculate WH_ij at just nonzero entries X_ij
    coo = X.tocoo()
    rows, cols, x_data = coo.row, coo.col, coo.data

    def build_Q():
        loop = time.perf_counter()



        # wh = np.zeros(len(rows))
        # for c in range(len(rows)):
        #     wh[c] = np.dot(W[rows[c], :], H[:, cols[c]])

        #wh = np.einsum("ik,ki->i", W[rows], H[:, cols])
        #wh = np.sum(W[rows, :] * H[:, cols].T, axis=1)

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



    
    #np.set_printoptions(precision = 2, suppress=True)

    # print(W @ H) # matrix multiplication
    # print(W)
    # sns.heatmap(W, annot=False, fmt=".2f", cmap="coolwarm", cbar=True)
    # plt.title("Matrix Heatmap Visualization")
    # plt.xlabel("Columns")
    # plt.ylabel("Rows")
    # plt.show()


    return W, H

# W, H dense arrays

def scKNIFE(X: sp.csr_matrix, R_g: np.ndarray, D: np.ndarray, A: np.ndarray,
             r: int, lam: float, beta_W: float, beta_H: float):
    L = D - A

    interp = True
    if interp:
        #RgT_X = R_g.T @ X
        #RgT_Rg = R_g.T @ R_g

        #DONE: implement the base NMF algorithm to 
        # determine the values of labeled nodes in W and H
        # 07/06/2026: 
        # read ESL book pages, derived NMF algorithm in P-14.23, 
        # proof from monotone convergence
        # implemented NMF

        W_f, H = NMF(X, r, 100)
        
        # 07/07/2026:
        # (L_uu + \epsilon I) W_u = - L_uf W_f
        # minimizes the Laplacian quadratic form from
        # block matrix [W_u W_f]^T
        #DONE: implement the interpolation

        i_u = L.shape[0] - len(W_f)
        W_u = np.linalg.solve(L[:i_u, :i_u] + np.eye(i_u) * 1e-6, 
                            np.dot(-L[:i_u, i_u:], W_f))
        W = np.vstack((W_u, W_f)) # type: ignore
    


    #DONE: implement the optimization 
    # with graph and regularization constraints
    for _ in range(100):
        pos = np.linalg.multi_dot([R_g.T, X, H.T]) + A @ W * lam
        neg = np.linalg.multi_dot([R_g.T, R_g, W, H, H.T]) \
        + D @ W * lam + beta_W + 1e-9
        W *= pos / neg
        
        pos = np.linalg.multi_dot([W.T, R_g.T, X])
        R_gW = R_g @ W
        neg = np.linalg.multi_dot([R_gW.T, R_gW, H]) + beta_H + 1e-9
        H *= pos / neg

    #TODO: implement the graph construction
    #TODO: think about integration of spatial structure

    return


R_g, D, A, X = get_matrices()

# START TIMER
start = time.perf_counter()

#NMF(np.ones((100, 100), dtype=int), 10, 200)

#W, H = scKNIFE(np.ones((100,100)), R_g, D, A, 1, 0.1, 0.05, 0.05)
W, H = NMF(X, r=5, iters=1000)

# END TIMER
print(f"{time.perf_counter() - start : .3f}", "NMF")

np.save("W_matrix.npy", W)
np.save("H_matrix.npy", H)
