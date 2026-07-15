import pandas as pd
import json
from ..config import RAW_DIR, SEP_DIR

def extract():
    df = pd.read_csv(RAW_DIR / "hgnc_complete_set.txt", 
                     sep='\t', 
                     usecols=[1, 19],
                     header=0)
    code = df["symbol"].to_numpy()
    ens = df["ensembl_gene_id"].to_numpy()
    hgnc_map = {}
    for i in range(len(code)):
        hgnc_map[ens[i]] = code[i]
    with open(SEP_DIR / "hgnc_map.json", "w") as f:
        json.dump(hgnc_map, f, indent=4)
if __name__ == "__main__":      
  extract()
