import pandas as pd
import inspect
import pathlib
import json
from ..config import SEP_DIR



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
    return name.lower().strip().replace(' ', "_")

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
    NODE_MAP = get_map(map_name)
    for src in edge_list["src_id"].to_numpy(dtype=str):
        add_node(src, NODE_MAP)
    for dst in edge_list["dst_id"].to_numpy(dtype=str):
        add_node(dst, NODE_MAP)
    save_map(NODE_MAP, map_name)


def add_node(node:str, map:dict[str, int]):
    """
    Add node to collective dict
    """
    if node not in map:
        map[node] = len(map)