# "Ensembl2Reactome_ALL_Levels.txt" data

# some of column 1 data is not the Ensemble identifier
# filter: $6 == "Homo sapiens" && $1 starts with "ENSG"

import pandas as pd
from . import utils
from ..config import RAW_DIR, SEP_DIR

def extract() -> pd.DataFrame:
    reactome_headers = [
        "ENSG_ID",
        #"Pathway_ID", # 1=Reactome Pathway ID
        "Pathway_Name", # Fuzzy name
        "Species"
    ]

    df = pd.read_csv(RAW_DIR / "Ensembl2Reactome_ALL_Levels.txt", 
                    sep='\t',
                    names=reactome_headers,
                    usecols=[0, 3, 5], # choosing 3=Pathway_Name
                    header=None)

    filter = (df["Species"] == "Homo sapiens") & df["ENSG_ID"].str.startswith("ENSG")
    df = (df[filter]).drop(columns=["Species"])
    df.iloc[:, 1] = df.iloc[:, 1].apply(utils.normalize)
    
    # rename columns, compile to edge list
    df.columns = ["src_id", "dst_id"]
    df["edge_type"] = "gene-pathway" 
    df["source"] = "reactome"

    df.to_csv(SEP_DIR / "reactome.tsv", sep='\t', index=False)

    utils.add_nodes(df)
    # pathway collision with Human-Gem

    return df
if __name__ == "__main__":   
  extract()




