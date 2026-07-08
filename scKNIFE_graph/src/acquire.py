import requests
import hashlib
import pathlib
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
    
    def __post_init__(self):
        if self.method == "http" and not self.url:
            raise ValueError(
                f"{self.name}: http source needs an URL")
        if self.method == "manual" and not self.notes:
            raise ValueError(f"{self.name}: manual source, need instructions")
    
    def __repr__(self):
        return (f"Source(name={self.name!r}), method = {self.method!r}, "
                f"url = {self.url!r}, raw_path = {self.raw_path!r}, "
                f"notes = {self.notes!r}")
    
    def __eq__(self, other):
        if other.__class__ is not self.__class__:
            return NotImplemented
        return (self.name, self.method, self.url, self.raw_path, self.notes)\
              == (other.name, other.method, other.url, other.raw_path, other.notes)

SOURCES = [
    Source("human_gem", 
           "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/refs/heads/main/model/Human-GEM.xml", 
           "http", 
           "Human-Gem.xml")
]

MANIFEST = {"human_gem": {}}

RAW_DIR = "/Users/jerryfan/Documents/Github/research2026/scKNIFE_graph/data/raw/"

def manifest(source: Source):
    source_dict = MANIFEST[source.name]
    source_dict["url"] = source.url
    source_dict["raw_path"] = source.raw_path
    source_dict["timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    path = pathlib.Path(RAW_DIR + source.raw_path)
    if not path.is_file():
        raise FileNotFoundError
    if source_dict["sha256"]:
        return
    with path.open("rb") as f:
        sha256_hash = hashlib.file_digest(f, "sha256")
    source_dict["sha256"] = sha256_hash.hexdigest()

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

print(MANIFEST)