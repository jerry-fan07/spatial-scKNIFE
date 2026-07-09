import pandas as pd
import pathlib
import json

SEP_DIR = "/Users/jerryfan/Documents/Github/research2026/spatial-scKNIFE/scKNIFE_graph/data/separated/"



def compile_list(edge_list : list, file_name : str):
    """
    Returns a pd.DataFrame object
    Creates a .tsv file from the DF

    edge_list = list of tuples: (src_id, dst_id, edge_type, sources)
    
    """
    df = pd.DataFrame(edge_list, 
                      columns=["src_id", "dst_id", "edge_type", "source"])
    df.to_csv(SEP_DIR + file_name, sep='\t', index=False)
    return df

def normalize(name: str) -> str:
    return name.lower().strip().replace(' ', "_")

def get_map(name: str) -> dict[str, int]:
    """
    Args:
        name: "*.json"
    Open JSON file as node to ID dictionary
    """
    file = pathlib.Path(SEP_DIR + name)
    if not file.is_file():
        return {}
    with open(SEP_DIR + name, "r") as f:
        return json.load(f)
    
def save_map(map: dict[str, int], name):
    """
    Args: 
        map: node to ID dictionary
        name: "*.json"
    Save node to ID dictionary as JSON file
    """
    with open(SEP_DIR + name, "w") as f:
        json.dump(map, f, indent=4)

def add_node(node:str, map:dict[str, int]):
    """
    Add node to collective dict
    """
    if node not in map:
        map[node] = len(map)