import numpy as np
import pandas as pd
import pytest

from mobo.data import clean, load_raw, sha256, split

OBJ = ["neg_delta_est", "f_osc"]


def toy():
    return pd.DataFrame({
        "raw_row": range(7),
        "smiles": [
            "CCO",                  # kept
            "not_a_smiles",         # step 1: parse fails
            "c1ccccc1",             # step 3: same molecule as the next row -> both dropped
            "C1=CC=CC=C1",
            "CCN",                  # step 4: NaN objective
            "CC(=O)O",              # kept
            "OC1C2OC3C(O)C1OC23",   # step 6: SAscore 4.97 > 4.5
        ],
        "neg_delta_est": [-0.5, -0.1, -0.2, -0.3, np.nan, -0.4, -0.6],
        "f_osc": [0.1, 0.2, 0.3, 0.3, 0.1, 0.05, 0.2],
    })


def test_clean_applies_each_step():
    df, counts = clean(toy(), OBJ, sa_max=4.5)

    assert counts == {"raw": 7, "parsed": 6, "unique": 4, "finite_objectives": 3,
                      "sascore_ok": 3, "sa_feasible": 2}
    assert df["canonical_smiles"].tolist() == ["CC(=O)O", "CCO"]  # sorted
    assert df["mol_id"].tolist() == [0, 1]
    assert df["raw_row"].tolist() == [5, 0]


def test_mol_id_ignores_raw_row_order():
    shuffled = toy().sample(frac=1, random_state=0)
    a, _ = clean(toy(), OBJ, sa_max=4.5)
    b, _ = clean(shuffled, OBJ, sa_max=4.5)
    pd.testing.assert_frame_equal(a, b)


def test_split_is_disjoint_sorted_and_deterministic():
    pool, dev = split(100, pool_size=10, pool_seed=42, dev_size=20, dev_seed=7)

    assert len(pool) == 10 and len(dev) == 20
    assert not set(pool) & set(dev)
    assert (np.diff(pool) > 0).all() and (np.diff(dev) > 0).all()

    pool2, dev2 = split(100, pool_size=10, pool_seed=42, dev_size=20, dev_seed=7)
    assert (pool == pool2).all() and (dev == dev2).all()


def test_load_raw_checks_hash_and_signs(tmp_path):
    path = tmp_path / "raw.csv"
    path.write_text("smiles,gap,osc\nCCO,0.25,0.1\n")
    cfg = {
        "raw_path": path,
        "raw_sha256": sha256(path),
        "smiles_column": "smiles",
        "objectives": [
            {"name": "neg_delta_est", "column": "gap", "sign": -1},
            {"name": "f_osc", "column": "osc", "sign": 1},
        ],
    }

    df = load_raw(cfg)
    assert df.loc[0, "neg_delta_est"] == -0.25 and df.loc[0, "f_osc"] == 0.1

    with pytest.raises(ValueError, match="sha256 mismatch"):
        load_raw({**cfg, "raw_sha256": "0" * 64})
