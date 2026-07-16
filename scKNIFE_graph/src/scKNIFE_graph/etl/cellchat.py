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

    # expand ligand/receptor complexes into constituent gene symbols
    edge_list = df[["ligand.symbol", "receptor.symbol"]].rename(
        columns={"ligand.symbol": "src_id", "receptor.symbol": "dst_id"})
    edge_list["src_id"] = edge_list["src_id"].str.split(",")
    edge_list["dst_id"] = edge_list["dst_id"].str.split(",")
    edge_list = edge_list.explode("src_id").explode("dst_id")
    edge_list["src_id"] = edge_list["src_id"].str.strip()
    edge_list["dst_id"] = edge_list["dst_id"].str.strip()
    edge_list = edge_list[(edge_list["src_id"] != "") &
                          (edge_list["dst_id"] != "") &
                          (edge_list["src_id"] != edge_list["dst_id"])]
    edge_list["src_id"] = [utils.gene_normalize(s) for s in edge_list["src_id"]]
    edge_list["dst_id"] = [utils.gene_normalize(d) for d in edge_list["dst_id"]]

    edge_list["edge_type"] = "ligand-receptor"
    edge_list["source"] = "cellchat"

    edge_list.to_csv(SEP_DIR / "cellchat.tsv", sep='\t', index=False)

    utils.add_nodes(edge_list)
    return edge_list

if __name__ == "__main__":
    extract()