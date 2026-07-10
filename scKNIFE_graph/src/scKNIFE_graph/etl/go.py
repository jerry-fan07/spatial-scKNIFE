import pandas as pd
import numpy as np
from . import utils
from ..config import RAW_DIR, SEP_DIR


gaf_columns = [
    "DB",
    "DB_Object_ID",
    "Gene", # ** DB_Object_Symbol
    "Qualifier",
    "GO_ID", # **
    "DB_Reference",
    "Evidence_Code", 
    "With_From", 
    "Aspect", # ** (filtered for P = Biological Processes)
    "DB_Object_Name", 
    "DB_Object_Synonym", 
    "DB_Object_Type", 
    "Taxon", 
    "Date", 
    "Assigned_By", 
    "Annotation_Extension", 
    "Gene_Product_Form_ID"
]

def extract() -> pd.DataFrame:
    """ Extracts data directly from .gaf.gz file to DF, TSV

    Pulls data into Gene and GO_ID dataframe

    Returns:
        DataFrame of src_id, dst_id, edge_type, source

    """
    NODE_MAP = utils.get_map("go_NM.json") # dict
    edge_list = []
    df = pd.read_csv(RAW_DIR / "HUMAN-uniprot.gaf.gz",
                      sep="\t", comment="!", 
                      names=gaf_columns, header=None)
    gene_term_df = df.loc[
        df["Aspect"] == "P",
          ["Gene", "GO_ID"]]
    
    src_array = gene_term_df["Gene"].to_numpy()
    dst_array = gene_term_df["GO_ID"].to_numpy()

    for i in range(len(src_array)):
        utils.add_node(src_array[i], NODE_MAP)
        utils.add_node(dst_array[i], NODE_MAP)
    
    gene_term_df["edge_type"] = "gene-term"
    gene_term_df["source"] = "go"
    gene_term_df = gene_term_df.rename(
        columns = {"Gene": "src_id", 
                   "Aspect": "dst_id"})


    gene_term_df.to_csv(SEP_DIR / "go.tsv", 
                        sep = '\t', index=False)
    
    utils.save_map(NODE_MAP, "go_NM.json")
    return gene_term_df