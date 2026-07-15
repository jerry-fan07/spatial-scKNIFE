"""
PanglaoDB, manual dicts, TNBC panel

"""

import pandas as pd
from . import utils
from ..config import RAW_DIR, SEP_DIR

def panglaoDB() -> pd.DataFrame:
    df = pd.read_csv(RAW_DIR / "PanglaoDB.tsv.gz",
                        sep='\t',
                        header=0)
    df = df[(df["sensitivity_human"] > 0.5) & 
            (df["species"]).str.contains("Hs") &
            (df["canonical marker"] == 1.0)]
    df =  pd.DataFrame({
        "src_id" : df["official gene symbol"],
        "dst_id" : df["cell type"].map(utils.normalize),
        "edge_type": "celltype-marker",
        "source": "panglao"
    })
    return df

# filler function
def TNBC() -> pd.DataFrame:
    return pd.DataFrame()

def extract():
    df_list = []
    for func in (panglaoDB, TNBC):
        df_list.append(func())
    df = pd.concat(df_list)

    df.to_csv(SEP_DIR / "markers.tsv", sep='\t', index=False)
    utils.add_nodes(df)
    return df

if __name__ == "__main__":
    extract()



