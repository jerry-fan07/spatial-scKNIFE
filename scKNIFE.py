import numpy as np
from numpy.random import default_rng
import seaborn as sns
import matplotlib.pyplot as plt
import time

# R_g maps genes to nodes in unified graph

def NMF(X: np.ndarray, r: int, iter: int):
    # generate matrices of 1's with dimension N*r, r*p
    #W = np.ones((len(X), r), dtype=float)
    #H = np.ones((r, len(X[0])), dtype=float)
    W = default_rng(10).random((len(X), r), dtype=float)
    H = default_rng(11).random((r, len(X[0])), dtype=float)

    for _ in range(iter):
        # update values of W and K:
        WH = np.dot(W, H)
        for i in range(len(W)):
            for k in range(len(W[0])):
                pos, neg = 0, 0
                for j in range(len(X[0])):
                    pos += H[k][j] * X[i][j] / (WH[i][j] + 1e-9)
                    neg += H[k][j] + 1e-9
                W[i][k] = W[i][k] * pos / neg

        WH = np.dot(W, H)
        for k in range(len(H)):
            for j in range(len(H[0])):
                pos, neg = 0, 0
                for i in range(len(X)):
                    pos += W[i][k] * X[i][j] / (WH[i][j] + 1e-9)
                    neg += W[i][k] + 1e-9
                H[k][j] = H[k][j] * pos / neg 
    #np.set_printoptions(precision = 2, suppress=True)

    print(W @ H)
    # print(W)
    # sns.heatmap(W, annot=False, fmt=".2f", cmap="coolwarm", cbar=True)
    # plt.title("Matrix Heatmap Visualization")
    # plt.xlabel("Columns")
    # plt.ylabel("Rows")
    # plt.show()


    return W, H


def scKNIFE(X: np.ndarray, R_g: np.ndarray, D: np.ndarray, A: np.ndarray,
             r: int, lam: float, beta_W: float, beta_H: float):
    L = D - A

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
    i_u = len(L) - len(W_f)
    W_u = np.linalg.solve(L[i_u, :i_u] + np.eye(i_u) * 1e-6, 
                          np.dot(-L[:-i_u, i_u:], W_f))
    W = np.vstack(W_u, W_f) # type: ignore

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

start = time.perf_counter()
NMF(np.ones((100, 100), dtype=int), 10, 200)
print(time.perf_counter() - start)

