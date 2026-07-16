import pandas as pd
import inspect
import pathlib
import json
import requests
import time
from ..config import SEP_DIR, PRO_DIR



def compile_list(edge_list : list, file_name: str = "") -> pd.DataFrame:
    """
    Returns a pd.DataFrame object
    Creates a .tsv file from the DF

    edge_list = list of tuples: (src_id, dst_id, edge_type, sources)
    file_name = saved tsv file: "*.tsv"
    
    """
    if not file_name:
        file_name = pathlib.Path(inspect.stack()[1].filename).stem
    df = pd.DataFrame(edge_list, 
                    columns=["src_id",
                            "dst_id",
                            "edge_type",
                            "source"])
    df.to_csv(SEP_DIR / f"{file_name}.tsv", sep='\t', index=False)
    add_nodes(df, f"{file_name}_NM.json")
    return df

# thinking about doing more to minimize duplication
def normalize(name: str) -> str:
    return name.lower().strip().strip('"').replace(' ', "_")

def gene_normalize(name: str) -> str:
    return name.strip().strip('"').upper()

def get_map(name: str) -> dict[str, int]:
    """
    Args:
        name: "*.json"
    Open JSON file as node to ID dictionary
    """
    file = pathlib.Path(SEP_DIR / name)
    if not file.is_file():
        return {}
    with open(SEP_DIR / name, "r") as f:
        return json.load(f)
    
def save_map(map: dict[str, int], name:str):
    """
    Args: 
        map: node to ID dictionary
        name: "*.json"
    Save node to ID dictionary as JSON file
    """
    with open(SEP_DIR / name, "w") as f:
        json.dump(map, f, indent=4)

def add_nodes(edge_list: pd.DataFrame, name: str = ""):
    """ Saves all nodes from a compiled dataframe as "file_NM.json"
    Args:
        edge_list: (src_id, dst_id, edge_type, source) DataFrame
    """
    map_name = pathlib.Path(inspect.stack()[1].filename).stem + "_NM.json"
    if name:
        map_name = name
    NODE_MAP = {}
    for src in edge_list["src_id"].to_numpy(dtype=str):
        if src not in NODE_MAP:
            NODE_MAP[src] = len(NODE_MAP)
    for dst in edge_list["dst_id"].to_numpy(dtype=str):
        if dst not in NODE_MAP:
            NODE_MAP[dst] = len(NODE_MAP)
    save_map(NODE_MAP, map_name)


def add_node(node:str, map:dict[str, int]):
    """ Add node to collective dict
    """
    if node not in map:
        map[node] = len(map)

def load_hgnc_map() -> dict[str, str]:
    """ Ensembl gene id to HGNC symbol lookup 
    """
    file = pathlib.Path(SEP_DIR / "hgnc_map.json")
    if not file.is_file():
        return {}
    with open(file, "r") as f:
        return json.load(f)
    
def load_complete_map() -> dict[str, str]:
    """ Get complete map of nodes
    """
    file = pathlib.Path(PRO_DIR / "complete_NM.json")
    if not file.is_file():
        return {}
    with open(file, "r") as f:
        return json.load(f)


def canonicalize(ids: pd.Series, hgnc_map: dict[str, str]) -> pd.Series:
    """Collapse different ids for the same gene onto one "canonical" node.
    """
    return ids.map(lambda x: hgnc_map.get(x, x))
