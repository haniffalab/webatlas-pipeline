import re, ast, pickle
import numpy as np
import pandas as pd
from scipy import sparse
from anndata import AnnData
import scanpy as sc

# ----------------- helpers -----------------
_DIGITS = re.compile(r"(\d+)")
def assign_new_features(
    adata,
    X_new,                           # (n_obs, old_vars + k)
    new_feature_names,               # list/Index length = old_vars + k
    feature_types=None               # optional pd.Series indexed by new_feature_names
    ):
    """
    Safely assign widened X and matching .var (always set var first).
    """
    # Ensure shapes match n_obs
    if X_new.shape[0] != adata.n_obs:
        raise ValueError(f"Row mismatch: X_new has {X_new.shape[0]} rows, adata has {adata.n_obs} obs")

    # 1) set VAR FIRST
    adata.var = pd.DataFrame(index=pd.Index(new_feature_names, name=adata.var.index.name))
    if feature_types is not None:
        # align and set
        adata.var["feature_type"] = feature_types.reindex(adata.var.index).astype("category")
    adata.var_names = adata.var.index

    # 2) then set X
    if not sparse.isspmatrix(X_new):
        # keep memory sensible if large
        from scipy import sparse as _sp
        X_new = _sp.csr_matrix(X_new)
    adata.X = X_new

    # sanity
    assert adata.n_vars == X_new.shape[1], \
        f"n_vars ({adata.n_vars}) != X_new.shape[1] ({X_new.shape[1]})"
    
def _normalize_obs_label(x):
    """Normalize obs_id entry to a string like 'cell_123'."""
    # unwrap 1-element containers
    if isinstance(x, (list, tuple, np.ndarray)) and len(x) == 1:
        x = x[0]
    # bytes → str
    if isinstance(x, (bytes, np.bytes_)):
        x = x.decode("utf-8", errors="ignore")
    # numpy scalar → int
    if isinstance(x, (np.integer,)):
        return f"cell_{int(x)}"
    if isinstance(x, (np.floating,)) and np.isfinite(x):
        return f"cell_{int(x)}"
    # string that might look like "['cell_123']"
    if isinstance(x, str):
        xs = x.strip()
        if (xs.startswith("[") and xs.endswith(("]", ")"))):
            try:
                v = ast.literal_eval(xs)
                return _normalize_obs_label(v)
            except Exception:
                return xs
        return xs
    return str(x)

def _parse_cell_positions(obs_ids):
    """Map 'cell_k' → (zero-based position = k-1) → row index in pickle matrix."""
    pos2row = {}
    for row_idx, raw in enumerate(obs_ids):
        s = _normalize_obs_label(raw)
        m = _DIGITS.findall(s)
        if not m:
            raise ValueError(f"Could not find an integer in obs id '{raw}' → '{s}'.")
        one_based = int(m[-1])
        if one_based <= 0:
            raise ValueError(f"Cell numbering must be 1-based; got {one_based} from '{raw}'.")
        pos2row[one_based - 1] = row_idx
    return pos2row

def _normalize_feature_name(x):
    """Normalize factor names to plain strings (hashable)."""
    if isinstance(x, (list, tuple, np.ndarray)):
        if len(x) == 1:
            return _normalize_feature_name(x[0])
        return "_".join(map(_normalize_feature_name, list(x)))
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode("utf-8", errors="ignore")
    if isinstance(x, (np.integer,)):
        return str(int(x))
    if isinstance(x, (np.floating,)):
        return str(int(x)) if np.isfinite(x) else "nan"
    return str(x)

# ----------------- main function -----------------
def add_factors_to_adata_by_position(
    adata: AnnData,
    factors_pkl: str,
    obs_ids_pkl: str,
    matrix_pkl: str,
    *,
    missing: str = "drop",        # "drop" or "fill_zeros"
    factor_prefix: str | None = None,
    allow_oob: str = "ignore",    # "ignore" or "error"
) -> AnnData:
    # load pickles
    with open(factors_pkl, "rb") as f:
        raw_factor_names = list(pickle.load(f))
    with open(obs_ids_pkl, "rb") as f:
        obs_ids = list(pickle.load(f))
    with open(matrix_pkl, "rb") as f:
        M = pickle.load(f)

    if not sparse.isspmatrix(M):
        raise TypeError("matrix_pkl must contain a scipy.sparse matrix")
    M = M.tocsr(copy=False)

    n_obs_factors, n_factors = M.shape
    factor_names = [_normalize_feature_name(x) for x in raw_factor_names]
    if factor_prefix:
        factor_names = [f"{factor_prefix}{nm}" for nm in factor_names]

    if len(factor_names) != n_factors:
        raise ValueError(f"len(factor_names)={len(factor_names)} != M.shape[1]={n_factors}")
    if len(obs_ids) != n_obs_factors:
        raise ValueError(f"len(obs_ids)={len(obs_ids)} != M.shape[0]={n_obs_factors}")

    n_obs_adata = adata.n_obs

    pos2row_full = _parse_cell_positions(obs_ids)
    bad = [p for p in pos2row_full if p < 0 or p >= n_obs_adata]
    if bad and allow_oob == "error":
        raise ValueError(f"{len(bad)} positions exceed AnnData rows (0..{n_obs_adata-1}). Examples: {bad[:5]}")

    # drop invalid
    pos2row = {p: j for p, j in pos2row_full.items() if 0 <= p < n_obs_adata}

    if missing not in {"drop", "fill_zeros"}:
        raise ValueError("missing must be 'drop' or 'fill_zeros'")

    if missing == "drop":
        keep_positions_sorted = np.array(sorted(pos2row.keys()), dtype=int)
        adata._inplace_subset_obs(keep_positions_sorted)
        n_obs_adata = adata.n_obs
        rows = np.arange(n_obs_adata, dtype=int)
        cols = np.array([pos2row[p] for p in keep_positions_sorted], dtype=int)
        data = np.ones_like(rows, dtype=float)
        S = sparse.csr_matrix((data, (rows, cols)), shape=(n_obs_adata, n_obs_factors))
    else:
        rows, cols, data = [], [], []
        for i in range(n_obs_adata):
            j = pos2row.get(i)
            if j is not None:
                rows.append(i); cols.append(j); data.append(1.0)
        S = sparse.csr_matrix((data, (rows, cols)), shape=(n_obs_adata, n_obs_factors))

    M_aligned = S @ M

    # concatenate features
    X0 = adata.X
    if sparse.isspmatrix(X0):
        X0 = X0.tocsr(copy=False)
    else:
        X0 = sparse.csr_matrix(X0)
    X_new = sparse.hstack([X0, M_aligned], format="csr")

    # rebuild .var
    gene_names = list(adata.var_names)
    combined_names = pd.Index(gene_names + factor_names)
    if combined_names.has_duplicates:
        counts = {}
        new_names = []
        for nm in combined_names:
            k = counts.get(nm, 0)
            new_names.append(nm if k == 0 else f"{nm}.{k}")
            counts[nm] = k + 1
        combined_names = pd.Index(new_names)

    feature_type = pd.Series(
        ["gene"] * len(gene_names) + ["factor"] * len(factor_names),
        index=combined_names, name="feature_type", dtype="category",
    )

    print("adata.n_obs:", adata.n_obs)
    print("old n_vars:", len(gene_names))
    print("n_factors:", len(factor_names))
    print("X_new shape:", X_new.shape)
    print("len(combined_names):", len(combined_names))
    adata = rebuild_with_new_features(
        adata,
        X_new=X_new,
        feature_names=combined_names,
        feature_type=feature_type
        )

    # optional sanity check
    #assert adata.var.shape[0] == X_new.shape[1]

    return adata

def rebuild_with_new_features(adata, X_new, feature_names, feature_type=None):
    """
    Return a *new* AnnData with X_new (n_obs × new_n_vars) and a fresh .var,
    preserving obs/obsm/uns/layers/obsp/varm from the input.
    """
    if X_new.shape[0] != adata.n_obs:
        raise ValueError(f"Row mismatch: X_new has {X_new.shape[0]} rows, adata has {adata.n_obs} obs")

    # Ensure names are an Index and lengths match
    feature_names = pd.Index(feature_names, name=(adata.var.index.name if adata.var is not None else None))
    if X_new.shape[1] != len(feature_names):
        raise ValueError(f"Col mismatch: X_new has {X_new.shape[1]} cols, feature_names has {len(feature_names)}")

    # Build new .var
    new_var = pd.DataFrame(index=feature_names)
    if feature_type is not None:
        ft = pd.Series(feature_type, index=feature_names, dtype="category")
        new_var["feature_type"] = ft

    # Ensure sparse CSR for large data
    if not sparse.isspmatrix(X_new):
        from scipy import sparse as _sp
        X_new = _sp.csr_matrix(X_new)

    # Create the new AnnData and copy sidecars
    new = sc.AnnData(X_new, obs=adata.obs.copy(), var=new_var)
    new.obsm = adata.obsm.copy()
    new.layers = dict(adata.layers)
    new.uns = adata.uns.copy()
    new.obsp = adata.obsp.copy()
    new.varm = adata.varm.copy()
    return new