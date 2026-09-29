# Sample-Efficient Multi-Objective Molecular Discovery

A finite-pool, bi-objective active search benchmark on Tartarus organic emitters
(maximize −ΔE_ST and f_osc under SAscore ≤ 4.5), comparing Random, Greedy-HVI,
ParEGO and qLogEHVI under a 500-evaluation offline-oracle budget.

Status: Project Start Gate passed; next is Week 1 data cleaning. Rules: [docs/benchmark_spec.md](docs/benchmark_spec.md). Plan: [docs/plan.md](docs/plan.md).

## Setup

```bash
conda env create -f environment.yml   # or environment.lock.yml for exact versions
conda activate mmo
python scripts/smoke_env.py
```
