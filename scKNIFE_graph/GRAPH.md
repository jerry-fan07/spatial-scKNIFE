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
- *Human-GEM*: https://sbml.org/software/libsbml/5.18.0/docs/formatted/python-api/index.html
  - Gene-reaction exists already, do we need gene-metabolite if reaction is already linked to the metabolites?
    - paper: “To improve traversability between expression-defined genes and metabolite space"
    - Should just solve this problem by incorporating weighted edges to the graph, this is doing essentially the same thing $(f_i - f_j)^2 + (f_j - f_k)^2 + (f_i - f_k)^2 = 2(f_i^2+f_j^2+f_k^2 - f_if_j-f_if_k-f_jf_k)$ 
  - Stoichiometric data should be used, and in addition, for `FbcAssociation` nodes, distinguish between `FbcAnd` versus `FbcOr`
  - Thinking about using fuzzy word matching algorithm or compute words into vectors and calculating distances to check for duplication within pathway names and currently only have normalize function -> `.lower().strip().replace(“ “, “_")`
