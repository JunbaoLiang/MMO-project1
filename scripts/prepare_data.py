"""Week 1: raw CSV -> cleaned, split molecules + dev-derived constants.

Spec §1-4 (clean and split), §8 (scaling), §9 (reference point).
Run from the repo root:  python scripts/prepare_data.py
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import yaml

from mobo.data import clean, load_raw, split

CONFIG = Path("configs/benchmark.yaml")
OUT = Path("data/processed")


def digest(strings):
    """Identity of an ordered list of strings (e.g. the pool's canonical SMILES)."""
    return hashlib.sha256("\n".join(strings).encode()).hexdigest()


def main():
    cfg = yaml.safe_load(CONFIG.read_text())
    names = [obj["name"] for obj in cfg["data"]["objectives"]]

    df, counts = clean(load_raw(cfg["data"]), names, cfg["data"]["sa_max"])
    s = cfg["split"]
    pool, dev = split(len(df), s["pool_size"], s["pool_seed"], s["dev_size"], s["dev_seed"])

    # df's index equals mol_id, so .loc with id arrays selects those molecules
    df["split"] = "unused"
    df.loc[pool, "split"] = "pool"
    df.loc[dev, "split"] = "dev"

    OUT.mkdir(parents=True, exist_ok=True)
    df[["mol_id", "raw_row", "canonical_smiles", "sascore", "split"]].to_parquet(
        OUT / "molecules.parquet", index=False)
    df[["mol_id", *names]].to_parquet(OUT / "labels.parquet", index=False)

    # §8-9: constants from dev labels only (numpy std defaults to ddof=0)
    y_dev = df.loc[dev, names].to_numpy()
    mean, std = y_dev.mean(axis=0), y_dev.std(axis=0)
    ref = np.quantile((y_dev - mean) / std, cfg["reference_point"]["quantile"], axis=0)

    # Tracked in Git: counts, identities and dev-only statistics. Never pool labels.
    report = {
        "spec_version": cfg["spec_version"],
        "raw_sha256": cfg["data"]["raw_sha256"],
        "counts": counts,
        "split_sizes": {k: int((df["split"] == k).sum()) for k in ["pool", "dev", "unused"]},
        "pool_smiles_sha256": digest(df.loc[pool, "canonical_smiles"]),
        "dev_smiles_sha256": digest(df.loc[dev, "canonical_smiles"]),
        "dev_constants": {
            "objectives": names,
            "mean": mean.tolist(),
            "std": std.tolist(),
            "reference_point": ref.tolist(),
        },
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
