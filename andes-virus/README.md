# Andes virus assets

This directory contains Andes virus reference assets for Loculus and Pathoplexus experiments.

Regenerate the directory with:

```bash
python andes-virus/scripts/generate_andes_virus.py --output-dir andes-virus
```

The script requires `biopython`. Minimizer generation uses
[`loculus-project/nextclade-sort-minimizers`](https://github.com/loculus-project/nextclade-sort-minimizers)
and requires that script's dependencies: `click`, `numpy`, and `pyyaml`.
