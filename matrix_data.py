import scipy.sparse as sp
from scipy.io import mmread
import numpy as np
import rdata
import pyreadr
from scKNIFE_graph.src.scKNIFE_graph.config import PRO_DIR

def get_matrices() -> tuple[sp.csr_matrix, sp.csr_matrix, sp.csr_matrix, sp.csr_matrix]:
    
    L = sp.load_npz(PRO_DIR / "laplacian.npz")
    d = np.asarray(L.sum(axis=1)).ravel()
    D = sp.diags(d, dtype=int)
    A =(D - L).tocsr()

    R_g = sp.load_npz(PRO_DIR / "Rg_matrix.npz")

    
    # bassezA_parsed = rdata.parser.parse_file("/Users/jerryfan/Documents/Github/research2026/spatial-scKNIFE" \
    # "/BASSEZ_A/BassezA_2021_33958794_3patients.rds")
    # bassezA = rdata.conversion.convert(bassezA_parsed)
    # print(type(bassezA))

    X = mmread("bassezA.mtx").tocsr() #type: ignore
    return (R_g, D, A, X)  #type: ignore


get_matrices()