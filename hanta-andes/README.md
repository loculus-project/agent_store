# Andes Virus [Hantavirus] asset generator

This directory contains the generator for Andes Virus reference assets used by Pathoplexus.
The generated assets are served from the Pathoplexus reference store:

https://pathoplexus.github.io/reference_store/hanta-andes/

Regenerate the directory with:

```bash
python hanta-andes/scripts/generate_andes_virus.py --output-dir /tmp/hanta-andes
```

The script requires `biopython`. Minimizer generation uses
[`loculus-project/nextclade-sort-minimizers`](https://github.com/loculus-project/nextclade-sort-minimizers)
and requires that script's dependencies: `click`, `numpy`, and `pyyaml`.
