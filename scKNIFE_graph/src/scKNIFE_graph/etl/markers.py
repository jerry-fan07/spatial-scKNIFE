"""
PanglaoDB, manual dicts, TNBC panel

"""

import pandas as pd
from . import utils
from ..config import RAW_DIR, SEP_DIR

def panglaoDB():
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
    df.to_csv(SEP_DIR / "panglao.tsv", sep='\t', index=False)
    utils.add_nodes(df, "panglao_NM.json")

    return df

panglaoDB()



