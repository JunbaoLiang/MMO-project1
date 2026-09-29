# Raw data

Raw files are not tracked in Git. Identity and provenance live in [manifest.yaml](manifest.yaml).

## Fetch and verify

```bash
curl -sSfL -o data/raw/gdb13.csv \
  https://raw.githubusercontent.com/the-matter-lab/Tartarus/5accbf7217e586095760db45ecb94078e1089584/datasets/gdb13.csv
shasum -a 256 data/raw/gdb13.csv   # must equal manifest sha256
chmod 444 data/raw/gdb13.csv
```

Never edit files in this directory. All cleaning writes to `data/processed/`.
