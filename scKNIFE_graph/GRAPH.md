# scKNIFE Graph Construction

## Four Stages:
### 1. Data Acquisition
- <u>acquire.py</u> file holds a dataset with a Source class for each of the 11 datasets, determining how the data is acquired from GitHub or web URLs, as well as pulling metadata into a MANIFEST.json file to record access timestamp, hash, etc.
- Pulled as raw binary data into the <u>/data/raw/</u> directory

### 2. Raw Data ETL
- Raw binary data pulled into a .tsv file through a python file for each source titled "etl_[source].py” and a shared set of "utils.py"
- .tsv files are composed of edge lists between nodes in the subgraph
- **Separated files:** A separate .tsv file is created of each of the data sources and merging occurs in a separate process for stage 3, instead of processing them all together. Advantage comes from updating a single dataset later on instead of all at once.
- **Reasoning for having this separate step:** scipy-sparse matrices are not optimized to be constructed incrementally (adding rows and columns)
  - Loses sources from edge data

**Possible design (not used):**
The raw data of edges from all of the sources are taken one by one in the format (node 1, node 2). We create a list and for each edge data point on the list, also assign an idX to each new node that appears in order (for example, for the first entry, node 1 will have id=1 and node 2 will have id=2) since they are the first distinct nodes. Then in the adjacency matrix we are building, assign the number $\sum_n2^n$ for the nth datasets used instead of a 1 on coordinates (x_i, x_j) and (x_j, x_i) for an edge, and 0 otherwise (then converting back later in the Laplacian).
- Drawback is the issue of not having clean and separated pipelines for each dataset

**Implementation Notes/Issues:**
- *Human-GEM*:
  LIBSBML: https://sbml.org/software/libsbml/5.18.0/docs/formatted/python-api/index.html
  - Gene-reaction exists already, do we need gene-metabolite if reaction is already linked to the metabolites?
    - paper: “To improve traversability between expression-defined genes and metabolite space"
    - Should just solve this problem by incorporating weighted edges to the graph, this is doing essentially the same thing $(f_i - f_j)^2 + (f_j - f_k)^2 + (f_i - f_k)^2 = 2(f_i^2+f_j^2+f_k^2 - f_if_j-f_if_k-f_jf_k)$ 
  - Stoichiometric data should be used, and in addition, for `FbcAssociation` nodes, distinguish between `FbcAnd` versus `FbcOr`
  - Thinking about using fuzzy word matching algorithm or compute words into vectors and calculating distances to check for duplication within pathway names and currently only have normalize function -> `.lower().strip().replace(“ “, “_")`
**Result: 220462, Paper target: 151430 (FLAG)**
- *GO Biological Process*:
  Pulled HUMAN-uniprot.gaf.gz which is decompressed and read through the `pd.read_csv()` function. The columns of interest are `DB_Object_Symbol` (Gene), `GO_ID` (Molecular Function, Biological Process, or Cellular Component), `Aspect` which is filtered to “P” for Biological Process. This therefore maps genes to their biological processes, which are shared by many of them.
From the paper: "To preserve comparability with default Slalom inference while keeping the analysis tractable at atlas scale, we applied a single pre-filter to this prior matrix and retained only gene sets with at least 85 represented genes."
We then also find gene term sizes and filter for >= 85 in GO BP data as well, although this is not mentioned in the paper.
**Result: 205322 edges, Paper target: 197436**

  - Why do we not use GO Molecular Function (gene product activity) and GO Cellular Component (cellular location)?
  - Should use evidence code filtering (EXP, IDA, IMP) -> (IBA, IRD) -> (ISS, ISO) -> (TAS, NAS) -> (IEA) and can implement them for different edge weights in the graph later if used.
- *Reactome*:
  - Gene-pathway edges probably have a lot of collisions with gene-pathway edges in *Human-Gem*, and currently the genes are expressed in Ensembl ID format (decide to convert to gene codes or convert all to ENS id).
**Result: 178495 edges, Paper target: 110530 (FLAG)**

- *TRRUST*:
Pulled gene-tf_gene edges from TRRUST raw data made into a data frame from a simple conversion.
**Result: 9396 edges, Paper target: 12,296 (SLIGHT FLAG)**

- *Cytopus*:
	- cellular processes, cellular identities
    - Dependency pin of 1.3.4 pins the data from this version
    - two versions of data: Cytopus_1.31nc_newcelltypes.txt is not used for graph (contains less entries), and we only pull from the database (this is also what the authors did)
**Result: 9821 + 620 = 10,441 (Exactly matches paper)**

- *Cell Chat*:
  - Used `rdata` module’s `read_rda()` function to extract into a dict of data frames
  - Pulled ligand-receptor edges from the cell chat database, filtering for categories Secreted Signaling, ECM-Receptor, Cell-Cell contact only. Expanded receptor complexes to their symbols (separated by commas) and taking the edges as a complete graph.
**Result: 2203, Paper target: 2287**

- *Hallmark*:
Pulled gene-gene_set edges from the Molecular Signaling Database (msigdb), within the Hallmark collection. Straightforward conversion into a data frame.
**Result: 7371, Paper target: 7322**

- *Tabula*:
  - Used only FACS cells
  - `avg_logFC > 1`, `pct.1 - pct2 > 0.3`
**Issue: 100K+ edges in data, need to collect only 355 edges**