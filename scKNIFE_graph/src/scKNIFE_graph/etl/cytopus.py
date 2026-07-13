import pandas as pd
import cytopus as cp
from . import utils
from ..config import RAW_DIR

def extract() -> pd.DataFrame:
    edge_list = []
    cpG = cp.KnowledgeBase()
    process_map = cpG.processes
    for process in process_map:
        genes = process_map[process]
        edge_list += [(gene, 
                        utils.normalize(process),
                        "gene-biological_process",
                        "cytopus")
                      for gene in genes]
    identity_map = cpG.identities
    for identity in identity_map:
        genes = identity_map[identity]
        edge_list += [(gene, 
                       utils.normalize(identity), 
                       "gene-cellular_identity", 
                       "cytopus")
                      for gene in genes]
    return utils.compile_list(edge_list) # also saves the nodes

extract()