import pandas as pd
import rdata
from ..config import RAW_DIR, SEP_DIR
from . import utils
import warnings
warnings.filterwarnings("ignore", 
                        category=UserWarning, 
                        module="rdata")

def extract():
    cc_data = rdata.read_rda(RAW_DIR / "CellChatDB.human.rda") # type = dict
    df: pd.DataFrame = cc_data["CellChatDB.human"]["interaction"] # type = pd.DataFrame
    df = df[df["annotation"] != "Non-protein Signaling"]
    edge_list: pd.DataFrame = df[["ligand", "receptor"]]
    edge_list.columns = ["src_id", "dst_id"]
    edge_list["edge_type"] = "ligand-receptor"
    edge_list["source"] = "cellchat"

    edge_list.to_csv(SEP_DIR / "cellchat.tsv", sep='\t', index=False)
    utils.add_nodes(edge_list)

    return edge_list

extract()