#!/usr/bin/env python3
"""Generate and verify BDBV reference assets.

The assets in ../reference.fasta and ../genes/*.fasta are headerless sequence
files consumed by Pathoplexus through [[URL:...]] references. This script
regenerates them from the official Nextclade BDBV dataset and checks the
checked-in files byte-for-byte.
"""

from __future__ import annotations

import argparse
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path


DATASET_BASE_URL = (
    "https://raw.githubusercontent.com/nextstrain/nextclade_data/"
    "master/data_output/nextstrain/orthoebolavirus/bdbv/unreleased"
)
GENES = ("NP", "VP35", "VP40", "GP", "GP_003", "VP30", "VP24", "L")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass(frozen=True)
class CdsFeature:
    name: str
    start: int
    end: int
    strand: str
    phase: int


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-base-url", default=DATASET_BASE_URL)
    parser.add_argument(
        "--write",
        action="store_true",
        help="Overwrite the checked-in bdbv reference files with regenerated output.",
    )
    args = parser.parse_args()

    repo_bdbv_dir = Path(__file__).resolve().parents[1]
    generated = generate_assets(args.dataset_base_url)

    if args.write:
        write_assets(repo_bdbv_dir, generated)

    mismatches = compare_assets(repo_bdbv_dir, generated)
    if mismatches:
        for mismatch in mismatches:
            print(mismatch)
        raise SystemExit(f"{len(mismatches)} BDBV reference asset(s) differ")

    print("All BDBV reference assets match byte-for-byte.")
    for relative_path, contents in generated.items():
        print(f"{relative_path}: {len(contents)} bytes")


def generate_assets(dataset_base_url: str) -> dict[Path, bytes]:
    reference_fasta = download_text(f"{dataset_base_url.rstrip('/')}/reference.fasta")
    genome_annotation = download_text(f"{dataset_base_url.rstrip('/')}/genome_annotation.gff3")

    reference = parse_fasta(reference_fasta)
    features = parse_cds_features(genome_annotation)

    assets = {Path("reference.fasta"): f"{reference}\n".encode()}
    for gene in GENES:
        cds_parts = features.get(gene)
        if not cds_parts:
            raise SystemExit(f"Could not find CDS feature for {gene}")
        protein = translate_cds(reference, cds_parts)
        assets[Path("genes") / f"{gene}.fasta"] = f"{protein}\n".encode()
    return assets


def download_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode()


def parse_fasta(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">")).upper()


def parse_cds_features(gff: str) -> dict[str, list[CdsFeature]]:
    features: dict[str, list[CdsFeature]] = {}
    for line in gff.splitlines():
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) != 9 or columns[2] != "CDS":
            continue
        attrs = parse_gff_attributes(columns[8])
        name = attrs.get("Name")
        if name is None:
            raise SystemExit(f"CDS feature has no Name attribute: {line}")
        features.setdefault(name, []).append(
            CdsFeature(
                name=name,
                start=int(columns[3]),
                end=int(columns[4]),
                strand=columns[6],
                phase=int(columns[7]),
            )
        )

    for cds_parts in features.values():
        cds_parts.sort(key=lambda feature: feature.start)
    return features


def parse_gff_attributes(attributes: str) -> dict[str, str]:
    parsed = {}
    for part in attributes.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        parsed[key] = urllib.parse.unquote(value)
    return parsed


def translate_cds(reference: str, features: list[CdsFeature]) -> str:
    if len({feature.strand for feature in features}) != 1:
        raise SystemExit(f"CDS has mixed strands: {[feature.name for feature in features]}")

    cds = "".join(
        reference[feature.start - 1 + feature.phase : feature.end]
        for feature in features
    )
    if features[0].strand == "-":
        cds = reverse_complement(cds)
    elif features[0].strand != "+":
        raise SystemExit(f"Unsupported CDS strand: {features[0].strand}")

    if len(cds) % 3:
        raise SystemExit(f"CDS length is not divisible by 3 for {features[0].name}")

    protein = "".join(CODON_TABLE[cds[index : index + 3]] for index in range(0, len(cds), 3))
    return protein[:-1] if protein.endswith("*") else protein


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_assets(base_dir: Path, assets: dict[Path, bytes]) -> None:
    for relative_path, contents in assets.items():
        path = base_dir / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(contents)


def compare_assets(base_dir: Path, generated: dict[Path, bytes]) -> list[str]:
    mismatches = []
    for relative_path, generated_bytes in generated.items():
        checked_in_path = base_dir / relative_path
        if not checked_in_path.exists():
            mismatches.append(f"Missing checked-in file: {relative_path}")
            continue
        checked_in_bytes = checked_in_path.read_bytes()
        if checked_in_bytes != generated_bytes:
            mismatches.append(
                f"Mismatch for {relative_path}: "
                f"checked-in={len(checked_in_bytes)} bytes generated={len(generated_bytes)} bytes"
            )
    return mismatches


if __name__ == "__main__":
    main()
