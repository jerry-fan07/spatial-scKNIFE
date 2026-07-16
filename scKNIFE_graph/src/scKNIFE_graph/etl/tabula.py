from pathlib import Path
import pandas as pd
from . import utils
from ..config import RAW_DIR, SEP_DIR

def extract():
    df_list = []
    cluster_map = {} # dict of docts
    
    for file in Path(RAW_DIR / "22_markers").iterdir():
        suffix = file.stem[-7:]
        prefix = file.stem[:-7]
        if (suffix == "classes" and prefix[:4] == "facs"):
            if (prefix in cluster_map):
                raise ValueError("Duplicate file")
            else:
                with open(file, "r") as f:
                    next(f)
                    cluster_map[prefix] = {
                        (int((parts := line.split(','))[0].strip('"'))-1): parts[1].strip().strip('"')
                        for line in f
                        if line.strip()  # skips empty lines and prevents crashes
                    }
                    
    for file in Path(RAW_DIR / "22_markers").iterdir():
        suffix = file.stem[-7:] 
        prefix = file.stem[:-7]
        if (suffix == "markers" and prefix in cluster_map):
            df = pd.read_csv(file, sep=",", usecols = [2, 3, 4, 6, 7], header=0)
            df = df.loc[(df["avg_logFC"] > 1) & (df["pct.1"] - df["pct.2"] > 0.3), ["cluster", "gene"]]
            #print(df["cluster"].unique(), prefix)
            df["cluster"] = df["cluster"].map(cluster_map[prefix])
            df["gene"] = [utils.gene_normalize(g) for g in df["gene"]]
            df_list.append(df)

    df = pd.concat(df_list)
    df = pd.DataFrame({
        "src_id": df["gene"],
        "dst_id": df["cluster"].map(utils.normalize),
        "edge_type": "celltype-marker",
        "source": "tabula"
    })
    df.to_csv(SEP_DIR / "tabula.tsv", sep='\t', index=False)
    utils.add_nodes(df)

    return df

if __name__ == "__main__":
    extract()
