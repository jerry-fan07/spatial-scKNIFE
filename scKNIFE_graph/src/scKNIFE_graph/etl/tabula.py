from pathlib import Path
import pandas as pd
import numpy as np
from . import utils
from ..config import RAW_DIR, SEP_DIR

def extract():
    df_list = []
    cluster_map = {} # dict of docts
    
    for file in Path(RAW_DIR / "22_markers").iterdir():
        suffix = file.stem[-7:]
        prefix = file.stem[:-7]
        if (suffix == "classes"):
            if (prefix in cluster_map):
                raise ValueError("Duplicate file")
            else:
                with open(file, "r") as f:
                    next(f)
                    cluster_map[prefix] = {
                        int((parts := line.split(','))[0].strip('"')): parts[1].strip().strip('"')
                        for line in f
                        if line.strip()  # 👈 Crucial: skips empty lines and prevents crashes!
                    }
                    
    for file in Path(RAW_DIR / "22_markers").iterdir():
        suffix = file.stem[-7:] 
        prefix = file.stem[:-7]
        if (suffix == "markers" and prefix in cluster_map):
            arr = cluster_map[prefix]
            df = pd.read_csv(file, sep=",", usecols = [6, 7], header=0)
            df["cluster"] = df["cluster"].map(cluster_map)
            df_list.append(df)

    #df = pd.concat(df_list)

    print(cluster_map)
    #print(cluster_map["droplet_Bladder_cell_ontology_class_"])
            

    df = pd.concat(df_list)
    df.to_csv(SEP_DIR / "tabula.tsv", sep='\t', index=False)
    utils.compile_list(df)

    return df
    
extract()
