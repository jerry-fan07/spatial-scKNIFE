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

    parser.add_argument("acquire_list", 
                       type=int, 
                       nargs="*",
                       choices=range(len(acquire.SOURCES)),
                        help=f"IDs to acquire: {source_labels}")
    parser.add_argument("-a", "--all", 
                       action="store_true", 
                       help="Acquire/reacquire all sources")
    args = parser.parse_args()

    if args.all and args.acquire_list:
        parser.error("Passed both list and --all flag")

    if args.all:
        print("Acquiring all sources from HTTP url.")
        main(list(range(len(acquire.SOURCES))))
    elif args.acquire_list:
        main(args.acquire_list)
    else:
        main([])

