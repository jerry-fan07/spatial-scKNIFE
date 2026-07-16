import json
import pandas as pd
import numpy as np
import scipy.sparse as sp
from .config import SEP_DIR, PRO_DIR
from .etl import utils



def write_genes():
    genes = []
    with open(PRO_DIR / "genes.json", "w") as f:
        json.dump(genes, f, indent=4)
    
def build_Rg():
    complete_map = utils.load_complete_map() 
    with open(PRO_DIR / "genes.json", "r") as f:
        gene_list = json.load(f)
    index = [complete_map[g] for g in gene_list]
    Rg_len = len(index)

    Rg = sp.coo_matrix(
        (np.ones(Rg_len), (range(Rg_len), index)),
        shape = (Rg_len, Rg_len),
        dtype = int
    ).tocsr()
    sp.save_npz(PRO_DIR / "Rg_matrix.npz", Rg)
    return Rg

def build_laplacian():

    # region pull map
    edge_files = [file for file in SEP_DIR.iterdir()
                if file.is_file() and file.suffix == ".tsv"]
    with open(PRO_DIR / "complete_NM.json", "r") as f:
        complete_map = json.load(f)
    hgnc_map = utils.load_hgnc_map()
    # endregion

    # region build laplacian
    srcs = []
    dsts = []
    length = len(complete_map)

    df_list = [pd.read_csv(edge_file, sep='\t', header=0) for edge_file in edge_files]
    for df in df_list:
        # Canonicalize with the same function merge.py
        # so every endpoint resolves to a node id
        src = utils.canonicalize(df["src_id"], hgnc_map).map(complete_map)
        dst = utils.canonicalize(df["dst_id"], hgnc_map).map(complete_map)
        if src.isna().any() or dst.isna().any():
            raise ValueError(
                "Edge endpoint missing from complete_NM.json"
            )
        srcs.append(src.to_numpy(dtype=int))
        dsts.append(dst.to_numpy(dtype=int))

    # symmetrizes the edge list
    # so that the graph is undirected
    rows = np.concatenate(srcs + dsts)
    cols = np.concatenate(dsts + srcs)

    A = sp.coo_matrix(
        (np.ones(len(rows)), (rows, cols)), 
        shape=(length,length), 
        dtype=int
    ).tocsr()
    A.data = np.ones_like(A.data, dtype=int)

    # check for symmetry
    if (A != A.T).nnz != 0:
        raise ValueError(
            "A is not symmetric"
        )
    
    # scipy sparse matrix needs .ravel() since
    # not automatically converted to 1D array
    d = np.asarray(A.sum(axis=1)).ravel()
    D = sp.diags(d, dtype=int)
    L = (D - A).tocsr()
    # endregion

    sp.save_npz(PRO_DIR / "laplacian.npz", L)
    return L

if __name__ == "__main__":
    build_Rg()
    #build_laplacian()