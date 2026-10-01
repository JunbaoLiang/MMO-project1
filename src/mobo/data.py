"""Load, clean and split the raw data. Section numbers refer to docs/benchmark_spec.md."""
import hashlib
import math
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDConfig, RDLogger

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402  (RDKit ships SAscore as a contrib script, not a module)

RDLogger.DisableLog("rdApp.*")  # invalid SMILES are expected; don't spam stderr


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_raw(cfg):
    """§1-2: verify the file hash, then return raw_row, smiles and the internal (maximized) objectives."""
    actual = sha256(cfg["raw_path"])
    if actual != cfg["raw_sha256"]:
        raise ValueError(f"sha256 mismatch for {cfg['raw_path']}: got {actual}")

    raw = pd.read_csv(cfg["raw_path"], float_precision="round_trip")  # labels exactly as written
    df = pd.DataFrame({"raw_row": raw.index, "smiles": raw[cfg["smiles_column"]]})
    for obj in cfg["objectives"]:
        df[obj["name"]] = obj["sign"] * raw[obj["column"]]
    return df


def canonical(smiles):
    mol = Chem.MolFromSmiles(smiles)  # includes sanitization
    return None if mol is None else Chem.MolToSmiles(mol)


def sa_score(smiles):
    try:
        return sascorer.calculateScore(Chem.MolFromSmiles(smiles))
    except Exception:  # §3 step 5: any failure means the molecule is dropped
        return math.nan


def clean(df, objective_names, sa_max):
    """§3: apply the six filters in order.

    Returns the cleaned frame (sorted by canonical SMILES, with mol_id = row number)
    and the row count after each step.
    """
    counts = {"raw": len(df)}

    df = df.assign(canonical_smiles=df["smiles"].map(canonical))
    df = df[df["canonical_smiles"].notna()]
    counts["parsed"] = len(df)

    df = df[~df["canonical_smiles"].duplicated(keep=False)]  # drop whole duplicate groups
    counts["unique"] = len(df)

    df = df[np.isfinite(df[objective_names]).all(axis=1)]
    counts["finite_objectives"] = len(df)

    df = df.assign(sascore=df["canonical_smiles"].map(sa_score))
    df = df[df["sascore"].notna()]
    counts["sascore_ok"] = len(df)

    df = df[df["sascore"] <= sa_max]
    counts["sa_feasible"] = len(df)

    df = df.sort_values("canonical_smiles").reset_index(drop=True)
    df.insert(0, "mol_id", df.index)
    return df[["mol_id", "raw_row", "canonical_smiles", "sascore", *objective_names]], counts


def split(n, pool_size, pool_seed, dev_size, dev_seed):
    """§4: draw the pool from mol_ids 0..n-1, then dev from the rest. Both returned sorted."""
    pool = np.random.default_rng(pool_seed).choice(n, pool_size, replace=False)
    rest = np.setdiff1d(np.arange(n), pool)  # sorted
    dev = np.random.default_rng(dev_seed).choice(rest, dev_size, replace=False)
    return np.sort(pool), np.sort(dev)
