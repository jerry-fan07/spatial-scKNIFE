import pandas as pd
from . import utils
from ..config import RAW_DIR, SEP_DIR

def extract() -> pd.DataFrame:

    trrust_headers = [
        "src_id", # transcription factor
        "dst_id" # target
    ]

    df = pd.read_csv(RAW_DIR / "trrust_rawdata.human.tsv", 
                     sep='\t',
                     names=trrust_headers,
                     usecols=[0, 1],
                     header=None)

    df["edge_type"] = "gene-tf_target"
    df["source"] = "trrust"

    utils.add_nodes(df)
    df.to_csv(SEP_DIR / "trrust.tsv", sep="\t", index=False)

    return df

extract()