import argparse

from .etl import hgnc, go, cellchat, human_gem, reactome, trrust, msigdb, markers
from . import acquire
from . import merge, matrix

def main(acquire_list: list[int]):
    acquire.acquire(acquire_list)
    hgnc.extract()
    for script in (go, cellchat, human_gem, reactome, trrust, msigdb): # (add markers)
        script.extract()

    markers.panglaoDB() # TODO: delete later when marker finished

    merge.build_complete_map()
    matrix.build_laplacian()
    
if __name__ == "__main__":
    source_labels = ", ".join(
        f"{i}={source.name}" for i, source in enumerate(acquire.SOURCES))
    parser = argparse.ArgumentParser()
    parser.add_argument("acquire_list", type=int, nargs="+",
                        help=f"IDs to acquire: {source_labels}")

    args = parser.parse_args()
    main(args.acquire_list)


