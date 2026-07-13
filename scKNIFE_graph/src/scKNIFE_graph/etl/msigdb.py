import pandas as pd
from . import utils
from ..config import RAW_DIR

def extract() -> pd.DataFrame:
    edge_list = []
    with open(RAW_DIR / "msigdb_hallmark.gmt", "r") as f:
        edge_list += [(utils.normalize(list[0]), list[i], "gene-gene_set", "msigdb") 
                      for line in f
                      if (list := line.strip('\t'))
                      for i in range(2, len(list))
        ]
    return utils.compile_list(edge_list)