import pandas as pd
import numpy as np
from . import utils
from ..config import RAW_DIR, SEP_DIR

def extract() -> pd.DataFrame:
    df = pd.read_csv(RAW_DIR / "Task_by_Gene.csv", sep=',', index_col=0)

    matrix = df.to_numpy(dtype=int)
    rows = np.array([utils.normalize(t) for t in df.index]) # tasks
    cols = np.array([utils.gene_normalize(g) for g in df.columns]) # genes

    row_indices, col_indices = np.nonzero(matrix)

    df = pd.DataFrame({
        "src_id" : cols[col_indices],
        "dst_id" : rows[row_indices],
        "edge_type" : "task-gene",
        "source" : "Task_by_Gene"
    })

    df.to_csv(SEP_DIR / "metabolic_tasks.tsv", sep='\t', index=False)
    utils.add_nodes(df)
    return df

if __name__ == "__main__":
    extract()

