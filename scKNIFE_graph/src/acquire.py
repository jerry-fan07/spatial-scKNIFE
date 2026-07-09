import requests
import hashlib
import pathlib
import json
from datetime import datetime

# source data schema
class Source:
    def __init__(self, name: str, url: str | None, method: str, 
                 raw_path:str, notes: str | None = None):
        self.name = name
        self.method = method
        self.url = url
        self.raw_path = raw_path
        self.notes = notes
        self.__post_init__()
    
    # requirements for http/manual
    def __post_init__(self):
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

# source acquisition directions
SOURCES = [
    Source("human_gem", 
           "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/refs/heads/main/model/Human-GEM.xml", 
           "http", 
           "Human-Gem.xml")
]

# dict of dicts: source acquisition metadata
MANIFEST = {"human_gem": {}}

# ../data/raw/ directory
RAW_DIR = "/Users/jerryfan/Documents/Github/research2026/spatial-scKNIFE/scKNIFE_graph/data/raw/"
# manifest 
def manifest(source: Source):
    source_dict = MANIFEST[source.name]
    source_dict["url"] = source.url
    source_dict["raw_path"] = source.raw_path
    source_dict["timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    path = pathlib.Path(RAW_DIR + source.raw_path)
    if not path.is_file():
        raise FileNotFoundError
    
    # write source metadata to json:
    with open(RAW_DIR + "MANIFEST.json", "w") as f:
        json.dump(MANIFEST, f, indent=4)

    # does not replace if sha256 data exists
    if source_dict["sha256"]:
        return
    with path.open("rb") as f:
        sha256_hash = hashlib.file_digest(f, "sha256")
    source_dict["sha256"] = sha256_hash.hexdigest()
    

def acquire():
    for source in SOURCES:
        sha256_hash = hashlib.sha256()
        if source.method == "http":
            # stream to process hash and write to memory
            with requests.get(source.url, stream=True) as r:
                r.raise_for_status() # checks for HTTP error
                with open(RAW_DIR + source.raw_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size = 65536):
                        if chunk: # filter keep-alive empty chunks
                            f.write(chunk)
                            sha256_hash.update(chunk)
            MANIFEST[source.name]["sha256"] = sha256_hash.hexdigest()
            manifest(source)
        else:
            print(f"{source.name} has no http link")


acquire()