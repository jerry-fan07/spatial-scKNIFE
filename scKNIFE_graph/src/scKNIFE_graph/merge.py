""" Merge file

Build the complete node index (complete_NM.json) from the end points of
every edge. Any id that appears in an edge is a node by construction,
and does not drift out of sync with the per-source *_NM.json maps

"""
import json
import pandas as pd
from .config import SEP_DIR, PRO_DIR
from .etl import utils

def build_complete_map():
    edge_files = [file for file in SEP_DIR.iterdir()
                if file.is_file() and file.suffix == ".tsv"]

    hgnc_map = utils.load_hgnc_map()

    complete_map = {}

    for edge_file in edge_files:
        df = pd.read_csv(edge_file, sep='\t', header=0)
        for col in ("src_id", "dst_id"):
            for node in utils.canonicalize(df[col], hgnc_map):
                # this actually makes {source}_NM.json files obsolete
                if node not in complete_map:
                    complete_map[node] = len(complete_map)



    with open(PRO_DIR / "complete_NM.json", "w") as f:
        json.dump(complete_map, f, indent=4)
    print("Merge completed.")
    print("Total Nodes:", len(complete_map))
    return complete_map

if __name__ == "__main__":
    build_complete_map()
