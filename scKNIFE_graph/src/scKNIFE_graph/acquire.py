""" Module for acquiring/pulling data from HTTP links

Holds a list with a Source class for each of the 11 datasets, 
determining how the data is acquired from GitHub or web URLs, 
as well as pulling metadata into a MANIFEST.json
file to record access timestamp, hash, etc.

acquire(src_index: list[int]) pulls data
into /data/raw/

"""



import requests
import hashlib
import pathlib
import json
import os
from datetime import datetime
from .config import RAW_DIR

# source data schema
class Source:
    """ Schema for acquiring data from each source

    Attributes:
        name: dataset name
        method: http / manual
        url: http url given if exists
        raw_path: file path
        notes: notes if manual
    """
    def __init__(self, name: str, url: str | None, method: str, 
                 raw_path:str = None, notes: str | None = None):
        self.name = name
        self.method = method
        self.url = url
        if self.url is None or raw_path is not None:
            self.raw_path = raw_path
        else:
            self.raw_path = self.url.rpartition('/')[-1]
        self.notes = notes
        self.__post_init__()
    
    # requirements for http/manual
    def __post_init__(self):
        """ Validating that source contains info

        Raises: 
            ValueError: http and no url, manual and no notes
        """
        if self.method == "http" and not self.url:
            raise ValueError(
                f"{self.name}: http source needs an URL")
        if self.method == "manual" and not self.notes:
            raise ValueError(f"{self.name}: manual source, need instructions")
    
    # toString
    def __repr__(self):
        return (f"Source(name={self.name!r}), method = {self.method!r}, "
                f"url = {self.url!r}, raw_path = {self.raw_path!r}, "
                f"notes = {self.notes!r}")
    
    # equals
    def __eq__(self, other):
        if other.__class__ is not self.__class__:
            return NotImplemented
        return (self.name, self.method, self.url, self.raw_path, self.notes)\
              == (other.name, other.method, other.url, other.raw_path, other.notes)

# source acquisition list
SOURCES = [
    # 0
    Source("human_gem", 
           "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/refs/heads/main/model/Human-GEM.xml", 
           "http", 
           "Human-Gem.xml"),
    # 1
    Source("go",
           "https://current.geneontology.org/annotations/gaf/HUMAN-uniprot.gaf.gz",
           "http",
           "HUMAN-uniprot.gaf.gz"),
    # 2
    Source("reactome",
           "https://download.reactome.org/97/Ensembl2Reactome_All_Levels.txt",
           "http",
           "Ensembl2Reactome_ALL_Levels.txt"
           ),
    # 3
    Source("trrust",
           "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
           "http"
           ),
    # 4
    Source("cytopus",
           None,
           "dependency",
           None),
    # 5
    Source("msigdb", 
           "https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
           "http_msigdb",
           "msigdb_hallmark.gmt"
           )
]

def manifest(source: Source, sha256_hash: str = ""):
    """ Writes an entry for the given source in MANIFEST.json

    Args:
        source: entry source
    Raises:
        FileNotFoundError: wrong file path
    """
    # retrieve JSON
    # (dict of dicts: source acquisition metadata)
    manifest_path = pathlib.Path(RAW_DIR / "MANIFEST.json")
    if manifest_path.is_file():
        with open(manifest_path, "r") as manifest_file:
            MANIFEST = json.load(manifest_file)
    else:
        MANIFEST = {}
    
    # update JSON file
    if source.name not in MANIFEST:
        MANIFEST[source.name] = {}
    source_dict = MANIFEST[source.name]
    source_dict["url"] = source.url
    source_dict["raw_path"] = source.raw_path
    source_dict["timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    raw_path = pathlib.Path(RAW_DIR / source.raw_path)
    # raise error if file does not exist
    if not raw_path.is_file():
        raise FileNotFoundError

    # compute file hash
    # does not replace if sha256 data exists
    if sha256_hash:
        source_dict["sha256"] = sha256_hash
    else:
        with raw_path.open("rb") as f:
            digest = hashlib.file_digest(f, "sha256")
        source_dict["sha256"] = digest.hexdigest()

    # save source metadata to json:
    with open(RAW_DIR / "MANIFEST.json", "w") as f:
        json.dump(MANIFEST, f, indent=4)
    

def acquire(src_index: list[int]):
    """ Acquires data from the source entry
    
    For each source, if an http link exists, method pulls
    data from the link and downloads into the /data/raw/
    directory. Also computes a hash of the file for
    MANIFEST.json, and writes an entry.

    Args:
        src_index: list of source indices to pull
    
    """

    for i in src_index:
        source = SOURCES[i]
        sha256_hash = hashlib.sha256()
        if source.method == "http":
            # stream to process hash and write to memory
            with requests.get(source.url, stream=True) as r:
                r.raise_for_status() # checks for HTTP error
                with open(RAW_DIR / source.raw_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size = 65536):
                        if chunk: # filter keep-alive empty chunks
                            f.write(chunk)
                            sha256_hash.update(chunk)
            manifest(source, sha256_hash.hexdigest())
        elif source.method == "http_msigdb":
            s = requests.Session()
            s.post("https://www.gsea-msigdb.org/gsea/login",
                   data={"username": "jerry_fan@brown.edu", 
                         "password": "password"})
            with s.get(source.url, stream=True) as r:
                r.raise_for_status() # HTTP error check
                if "login.jsp" in r.url:
                    raise PermissionError("Failed login")
                with open(RAW_DIR / source.raw_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size = 65536):
                        if chunk: # filter keep-alive empty chunks
                            f.write(chunk)
                            sha256_hash.update(chunk)
            manifest(source, sha256_hash.hexdigest())

        else:
            print(f"{source.name} has no http link")


#TODO: write a method that records history and allows retraction of a data pull
# maybe Git is sufficient

# check it is run by python -m scKNIFE_graph/acquire.py
if __name__ == "__main__":
    acquire(src_index=[5])