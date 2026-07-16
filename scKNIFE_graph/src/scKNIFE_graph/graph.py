import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.sparse as sp

from .config import PRO_DIR

LAYOUT_CACHE = PRO_DIR / "layout.npy"


def load_adjacency() -> sp.csr_matrix:
    """Reconstruct the adjacency from the cached Laplacian.

    L = D - A with diag(A) = 0, so A = diag(L) - L. This avoids re-running
    the full ETL just to draw the graph.
    """
    L = sp.load_npz(PRO_DIR / "laplacian.npz")
    A = sp.diags(L.diagonal(), dtype=L.dtype) - L
    return A.tocsr()


def compute_layout(G: nx.Graph) -> np.ndarray:
    """Return an (n, 2) array of node positions ordered by node id.

    The Fruchterman-Reingold layout on ~70k nodes takes minutes, so the
    result is cached to disk and saved *before* drawing: a killed draw
    never loses the expensive part, and reruns are instant.
    """
    if LAYOUT_CACHE.exists():
        return np.load(LAYOUT_CACHE)

    # len(G) >= 500 -> networkx uses the sparse solver (no dense n*n array).
    # iterations capped and seed fixed for a reproducible, bounded layout.
    pos = nx.spring_layout(G, iterations=15, seed=0)
    coords = np.array([pos[i] for i in range(G.number_of_nodes())])
    np.save(LAYOUT_CACHE, coords)
    return coords


def main(draw_edges: bool = False) -> None:
    A = load_adjacency()
    G = nx.from_scipy_sparse_array(A)
    coords = compute_layout(G)

    fig, ax = plt.subplots(figsize=(12, 12))
    if draw_edges:
        # ~470k segments: slow to render and mostly a smear, so opt-in only.
        nx.draw_networkx_edges(
            G,
            pos={i: coords[i] for i in range(len(coords))},
            ax=ax,
            alpha=0.02,
            width=0.1,
        )
    ax.scatter(coords[:, 0], coords[:, 1], s=1, c="lightblue")
    ax.set_axis_off()
    plt.show()


if __name__ == "__main__":
    main(True)
