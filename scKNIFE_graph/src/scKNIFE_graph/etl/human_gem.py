import pandas as pd
import libsbml as sbml
from . import utils
from ..config import RAW_DIR

def extract() -> pd.DataFrame:
    """
    Returns: pd.DataFrame
    (src_id, dst_id, edge_type, source)
    
    Extracts raw SBML data from Human-Gem.xml file containing
    reaction data, using libsbml as a parser.
    """
    NODE_MAP = utils.get_map("human_gem_NM.json")
    reader = sbml.SBMLReader()
    
    hg_data = reader.readSBML(RAW_DIR / "Human-Gem.xml")
    model = hg_data.getModel()
    fbc_model = model.getPlugin("fbc")
    edge_list = []

    # GPA tree only holds ENSG IDs, gene symbols held in GeneProduct
    symbol_of = {gp.getId(): gp.getLabel()
                 for gp in fbc_model.getListOfGeneProducts()}
    
    # pathway-reaction edges
    # reaction nodes are added in larger for loop
    for group in model.getPlugin("groups").getListOfGroups():
        edge_list += [(utils.normalize(group.getName()), member.getIdRef(), 
                    "pathway-reaction", "human_gem")
                    for member in group.getListOfMembers()]

    # loop over each reaction
    # --> edges acquired: gene-reaction, reaction-metabolite, gene-metabolite
    for i in range(model.getNumReactions()):
        reaction = model.getReaction(i)
        utils.add_node(reaction.getId(), NODE_MAP)

        def collect_leaves(fbc_assoc: sbml.FbcAssociation) -> list[str]:
            """
            Returns: ENSG strings for genes of the tree

            Recursive method which appends FbcAssociation leaves
            (GeneProductRef objects)
            """
            # ! flattened AND/OR nodes (and/or = all/one gene required)
            if fbc_assoc.isGeneProductRef():
                return [fbc_assoc.getGeneProduct()]
            leaves = []
            for child in fbc_assoc.getListOfAssociations():
                leaves += collect_leaves(child)
            return leaves

        genes = []
        gpa = reaction.getPlugin("fbc").getGeneProductAssociation()
        if gpa is not None:
            genes = collect_leaves(gpa.getAssociation()) # node in fbc tree

        metabolites = []
        
        # gene-reaction edges
        for g in genes:
            gene = symbol_of[g]
            utils.add_node(gene, NODE_MAP)
            edge_list += [(gene, reaction.getId(), 
                           "gene-reaction", "human_gem")]

        # reaction-metabolite edges
        # think about meaning of edges, adding in stoichiometric data later
        for p in range(reaction.getNumProducts()):
            prod = reaction.getProduct(p).getSpecies()
            utils.add_node(prod, NODE_MAP)
            metabolites.append(prod)
            edge_list.append((prod, reaction.getId(),
                               "reaction-metabolite", "human_gem"))
        
        for r in range(reaction.getNumReactants()):
            react = reaction.getReactant(r).getSpecies()
            utils.add_node(react, NODE_MAP)
            metabolites.append(react)
            edge_list.append((react, reaction.getId(), 
                              "reaction-metabolite", "human_gem")) 

        # gene-metabolite edges
        # ! note: gene - reaction already exists, do we need gene - metabolite?
        # no need to update NODE_MAP
        for g in genes:
            gene = symbol_of[g]
            edge_list += [(gene, m, "gene-metabolite", "human_gem") 
                          for m in metabolites]
        
    utils.save_map(NODE_MAP, "human_gem_NM.json")
    # create pd.DataFrame and TSV
    edge_list = utils.compile_list(edge_list)

    return edge_list
