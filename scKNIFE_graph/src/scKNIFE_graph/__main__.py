import argparse

from .etl import hgnc, go, cellchat, human_gem, reactome, trrust, msigdb, markers, tabula
from . import acquire
from . import merge, matrix
import time

def main(acquire_list: list[int]):
    t = time.perf_counter()
    acquire.acquire(acquire_list)
    print(f"{time.perf_counter()-t : .3f}")
    t =time.perf_counter()
    
    hgnc.extract()
    for script in (go, 
                   cellchat, 
                   human_gem, 
                   reactome, 
                   trrust, 
                   msigdb, 
                   markers, 
                   tabula):
        script.extract()
        print(f"{time.perf_counter()-t : .3f}", script.__name__)
        t =time.perf_counter()

    merge.build_complete_map()
    print(f"{time.perf_counter()-t : .3f}", "merge")
    t =time.perf_counter()

    matrix.build_laplacian()
    print(f"{time.perf_counter()-t : .3f}", "build laplacian")
    t =time.perf_counter()
    
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

