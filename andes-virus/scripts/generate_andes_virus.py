#!/usr/bin/env python3
"""Generate the Andes virus assets stored in this repository.

The script fetches the three RefSeq records used for Andes virus, writes
headerless reference and gene sequence files, creates minimal Nextclade
datasets, builds a segment minimizer, and writes the dataset-server index.
"""

from __future__ import annotations

import argparse
import copy
import json
import shutil
import sys
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path
from textwrap import wrap

try:
    import numpy as np
    from Bio import Entrez, SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError as error:  # pragma: no cover
    raise SystemExit(
        "Missing dependencies. Install with: python -m pip install biopython numpy"
    ) from error


TAXON_ID = 1980456
DEFAULT_EMAIL = "loculus-agent@example.com"
DEFAULT_IMAGE_URL = "http://theo:8765/copper_r.png"
DATASET_COMPATIBILITY = {"cli": "3.0.0-alpha.0", "web": "3.0.0-alpha.0"}
MINIMIZER_ALGO_VERSION = "1"
MINIMIZER_JSON_SCHEMA_VERSION = "3.0.0"
MINIMIZER_K = 17


@dataclass(frozen=True)
class GeneSpec:
    name: str
    protein_id: str
    gff_attributes: tuple[str, ...]
    genbank_qualifiers: tuple[str, ...]


@dataclass(frozen=True)
class SegmentSpec:
    segment: str
    accession: str
    reference_description: str
    genes: tuple[GeneSpec, ...]


SEGMENTS: tuple[SegmentSpec, ...] = (
    SegmentSpec(
        "L",
        "NC_003468.2",
        "Andes virus segment L",
        (
            GeneSpec(
                "RdRp",
                "NP_604473.1",
                ("Name", "gbkey", "locus_tag", "protein_id", "product", "ID", "Dbxref"),
                ("ID", "Dbxref", "Name", "gbkey", "locus_tag", "product", "protein_id"),
            ),
        ),
    ),
    SegmentSpec(
        "M",
        "NC_003467.2",
        "Andes virus segment M",
        (
            GeneSpec(
                "GPC",
                "NP_604472.1",
                ("Name", "gbkey", "locus_tag", "protein_id", "ID", "Dbxref", "product"),
                ("ID", "Dbxref", "Name", "gbkey", "locus_tag", "product", "protein_id"),
            ),
        ),
    ),
    SegmentSpec(
        "S",
        "NC_003466.1",
        "Andes virus segment S",
        (
            GeneSpec(
                "N",
                "NP_604471.1",
                ("Name", "gbkey", "locus_tag", "protein_id", "ID", "product", "Dbxref"),
                ("ID", "Dbxref", "Name", "gbkey", "locus_tag", "product", "protein_id"),
            ),
            GeneSpec(
                "NSs",
                "YP_004928151.1",
                ("Name", "gbkey", "locus_tag", "protein_id", "ID", "product", "Dbxref"),
                ("ID", "Dbxref", "Name", "gbkey", "locus_tag", "product", "protein_id"),
            ),
        ),
    ),
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="andes-virus", type=Path)
    parser.add_argument("--email", default=DEFAULT_EMAIL)
    parser.add_argument("--image-url", default=DEFAULT_IMAGE_URL)
    parser.add_argument("--skip-image", action="store_true")
    parser.add_argument("--threshold", default=32, type=int)
    args = parser.parse_args()

    if not 1 <= args.threshold <= 32:
        raise SystemExit("--threshold must be between 1 and 32")

    Entrez.email = args.email
    output_dir = args.output_dir
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)
    (output_dir / "genes").mkdir()

    segment_records = {}
    for spec in SEGMENTS:
        record = fetch_genbank_record(spec.accession)
        segment_records[spec.segment] = record
        write_text(output_dir / f"{spec.accession}.fasta", str(record.seq).upper())
        write_nextclade_dataset(output_dir, spec, record)
        for gene in spec.genes:
            write_text(output_dir / "genes" / f"{gene.name}.fasta", extract_cds_translation(record, gene))

    write_nextclade_index(output_dir / "nextclade_data" / "index.json")
    write_minimizer(output_dir / "segment-minimizer.json", segment_records, args.threshold)
    if not args.skip_image:
        image_path = output_dir / "images" / "copper_r.png"
        image_path.parent.mkdir()
        download_file(args.image_url, image_path)


def fetch_genbank_record(accession: str) -> SeqRecord:
    with Entrez.efetch(db="nuccore", id=accession, rettype="gbwithparts", retmode="text") as handle:
        return SeqIO.read(handle, "genbank")


def write_nextclade_dataset(output_dir: Path, spec: SegmentSpec, record: SeqRecord) -> None:
    dataset_dir = output_dir / "nextclade" / spec.segment
    dataset_dir.mkdir(parents=True)
    fasta = format_fasta(spec.accession, record.description, str(record.seq).upper(), line_width=60)
    files = {
        "CHANGELOG.md": "# Changelog\n\n## 2026-05-08\n\nInitial minimal dataset generated from NCBI GenBank reference annotation.\n",
        "README.md": f"# Andes virus {spec.segment} segment\n\nMinimal Nextclade dataset for Andes virus {spec.segment} segment using Andes virus reference.\n",
        "genome_annotation.gff3": genbank_to_gff3(record),
        "pathogen.json": json.dumps(pathogen_json(spec), indent=2) + "\n",
        "reference.fasta": fasta,
        "reference.gb": nextclade_genbank(record, spec),
        "sequences.fasta": fasta,
    }
    for filename, content in files.items():
        write_text(dataset_dir / filename, content)

    release_dir = output_dir / "nextclade_data" / "andes-virus" / spec.segment / "unreleased"
    release_dir.mkdir(parents=True)
    with zipfile.ZipFile(release_dir / "dataset.zip", "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for filename in (
            "CHANGELOG.md",
            "README.md",
            "genome_annotation.gff3",
            "pathogen.json",
            "reference.fasta",
            "sequences.fasta",
        ):
            archive.write(dataset_dir / filename, arcname=filename)


def pathogen_json(spec: SegmentSpec) -> dict:
    return {
        "$schema": "https://raw.githubusercontent.com/nextstrain/nextclade/refs/heads/release/packages/nextclade-schemas/input-pathogen-json.schema.json",
        "schemaVersion": "3.0.0",
        "files": {
            "reference": "reference.fasta",
            "pathogenJson": "pathogen.json",
            "genomeAnnotation": "genome_annotation.gff3",
            "examples": "sequences.fasta",
            "readme": "README.md",
            "changelog": "CHANGELOG.md",
        },
        "attributes": {
            "name": f"Andes virus {spec.segment}",
            "reference name": spec.reference_description,
            "reference accession": spec.accession,
        },
        "qc": {
            "missingData": {"enabled": True, "missingDataThreshold": 2000, "scoreBias": 500},
            "snpClusters": {"enabled": False},
            "mixedSites": {"enabled": True, "mixedSitesThreshold": 15},
            "frameShifts": {"enabled": True, "scoreWeight": 20},
            "stopCodons": {"enabled": True, "scoreWeight": 50},
        },
    }


def write_nextclade_index(path: Path) -> None:
    datasets = []
    for spec in SEGMENTS:
        datasets.append(
            {
                "path": f"andes-virus/{spec.segment}",
                "shortcuts": [f"andes-virus/{spec.segment}"],
                "enabled": True,
                "attributes": {
                    "name": f"Andes virus {spec.segment} segment",
                    "reference name": spec.reference_description,
                    "reference accession": spec.accession,
                },
                "files": {
                    "changelog": "CHANGELOG.md",
                    "examples": "sequences.fasta",
                    "genomeAnnotation": "genome_annotation.gff3",
                    "pathogenJson": "pathogen.json",
                    "readme": "README.md",
                    "reference": "reference.fasta",
                },
                "capabilities": {"qc": ["frameShifts", "missingData", "mixedSites", "stopCodons"]},
                "versions": [{"tag": "unreleased", "compatibility": DATASET_COMPATIBILITY}],
                "version": {"tag": "unreleased", "compatibility": DATASET_COMPATIBILITY},
            }
        )
    write_json(
        path,
        {
            "schemaVersion": "3.0.0",
            "collections": [
                {
                    "meta": {
                        "id": "loculus-agent-store",
                        "title": "Loculus agent store datasets",
                        "description": "Minimal Nextclade datasets staged by automation for Loculus testing.",
                        "maintainers": [{"name": "Loculus project", "url": "https://github.com/loculus-project"}],
                        "urls": [{"name": "source", "url": "https://github.com/loculus-project/agent_store"}],
                    },
                    "datasets": datasets,
                }
            ],
        },
    )


def genbank_to_gff3(record: SeqRecord) -> str:
    lines = [
        "##gff-version 3",
        "#!gff-spec-version 1.21",
        "#!processor NCBI annotwriter",
        f"##sequence-region {record.id} 1 {len(record.seq)}",
        f"##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={TAXON_ID}",
    ]
    for feature in record.features:
        if feature.type not in {"source", "CDS"}:
            continue
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        strand = "+" if feature.location.strand != -1 else "-"
        phase = "0" if feature.type == "CDS" else "."
        qualifiers = feature.qualifiers
        attributes = []
        if feature.type == "source":
            segment = first(qualifiers, "segment")
            attributes.extend(
                [
                    f"ID={record.id}:1..{len(record.seq)}",
                    f"Dbxref=taxon:{TAXON_ID}",
                    f"Name={segment}",
                    "gbkey=Src",
                    "genome=genomic",
                ]
            )
            add_attribute(attributes, "mol_type", qualifiers)
            attributes.append("old-name=Andes virus")
            for key in ("segment", "strain"):
                add_attribute(attributes, key, qualifiers)
        else:
            gene = gene_for_feature(feature)
            values = gff_cds_values(feature, gene)
            for key in gene.gff_attributes:
                add_value(attributes, key, values.get(key))
        lines.append(
            "\t".join(
                [
                    record.id,
                    "RefSeq",
                    "region" if feature.type == "source" else "CDS",
                    str(start),
                    str(end),
                    ".",
                    strand,
                    phase,
                    ";".join(attributes),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def extract_cds_translation(record: SeqRecord, gene: GeneSpec) -> str:
    for feature in record.features:
        if feature.type == "CDS" and gene.protein_id in feature.qualifiers.get("protein_id", []):
            translation = first(feature.qualifiers, "translation")
            if translation:
                return translation + "*"
            return str(feature.extract(record.seq).translate(to_stop=False))
    raise ValueError(f"Could not find CDS {gene.name} ({gene.protein_id}) in {record.id}")


def nextclade_genbank(record: SeqRecord, spec: SegmentSpec) -> str:
    record = copy.deepcopy(record)
    cds_features = []
    for feature in record.features:
        if feature.type == "source":
            cds_features.append(feature)
        elif feature.type == "CDS":
            gene = gene_for_feature(feature)
            values = genbank_cds_values(feature, gene)
            feature.qualifiers.clear()
            for key in gene.genbank_qualifiers:
                value = values.get(key)
                if value is not None:
                    feature.qualifiers[key] = [value]
            cds_features.append(feature)
    record.features = cds_features
    return record.format("genbank")


def gene_for_feature(feature) -> GeneSpec:
    protein_id = first(feature.qualifiers, "protein_id")
    for spec in SEGMENTS:
        for gene in spec.genes:
            if gene.protein_id == protein_id:
                return gene
    raise ValueError(f"Unexpected CDS protein_id {protein_id}")


def gff_cds_values(feature, gene: GeneSpec) -> dict[str, str]:
    product = first(feature.qualifiers, "product") or ""
    gene_id = gene_id_from_feature(feature)
    return {
        "Name": gene.name,
        "gbkey": "CDS",
        "locus_tag": first(feature.qualifiers, "locus_tag") or "",
        "protein_id": gene.protein_id,
        "product": product,
        "ID": f"cds-{gene.protein_id}",
        "Dbxref": f"GenBank:{gene.protein_id},{gene_id}",
    }


def genbank_cds_values(feature, gene: GeneSpec) -> dict[str, str]:
    values = gff_cds_values(feature, gene)
    values["locus_tag"] = gene.name
    return values


def gene_id_from_feature(feature) -> str:
    for db_xref in feature.qualifiers.get("db_xref", []):
        if db_xref.startswith("GeneID:"):
            return db_xref
    raise ValueError(f"Missing GeneID for {first(feature.qualifiers, 'protein_id')}")


def write_minimizer(path: Path, records_by_segment: dict[str, SeqRecord], threshold: int) -> None:
    cutoff = 1 << threshold
    refs = {
        segment: SeqRecord(Seq(str(records_by_segment[segment].seq).upper()), id=segment, description="")
        for segment in ("L", "M", "S")
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(serialize_ref_search_index(make_ref_search_index(refs, cutoff))), encoding="utf-8")


def make_ref_search_index(refs: dict[str, SeqRecord], cutoff: int) -> dict:
    minimizers_by_reference = []
    for name, ref in sorted(refs.items()):
        minimizers = get_ref_search_minimizers(ref, cutoff)
        minimizers_by_reference.append(
            {
                "minimizers": minimizers,
                "meta": {"name": name, "length": len(ref.seq), "nMinimizers": len(minimizers)},
            }
        )
    index = {"minimizers": {}, "references": []}
    for ref_index, minimizer_set in enumerate(minimizers_by_reference):
        for minimizer in minimizer_set["minimizers"]:
            index["minimizers"].setdefault(minimizer, []).append(ref_index)
        index["references"].append(minimizer_set["meta"])
    normalization = np.array([x["length"] / x["nMinimizers"] for x in index["references"]])
    return {
        "schemaVersion": MINIMIZER_JSON_SCHEMA_VERSION,
        "version": MINIMIZER_ALGO_VERSION,
        "params": {"k": MINIMIZER_K, "cutoff": cutoff},
        **index,
        "normalization": normalization,
    }


def get_ref_search_minimizers(seq: SeqRecord, cutoff: int) -> np.ndarray:
    seq_str = str(seq.seq).upper().replace("-", "")
    minimizers = []
    for i in range(len(seq_str) - MINIMIZER_K):
        minimizer_hash = get_hash(seq_str[i : i + MINIMIZER_K], cutoff)
        if minimizer_hash < cutoff:
            minimizers.append(minimizer_hash)
    return np.unique(minimizers)


def get_hash(kmer: str, cutoff: int) -> int:
    value = 0
    offset = 0
    for index, nuc in enumerate(kmer):
        if index % 3 == 2:
            continue
        if nuc not in "ACGT":
            return cutoff + 1
        if nuc in "AC":
            value += 1 << offset
        if nuc in "AT":
            value += 1 << (offset + 1)
        offset += 2
    return invertible_hash(value)


def invertible_hash(value: int) -> int:
    mask = (1 << 32) - 1
    value = (~value + (value << 21)) & mask
    value = value ^ (value >> 24)
    value = (value + (value << 3) + (value << 8)) & mask
    value = value ^ (value >> 14)
    value = (value + (value << 2) + (value << 4)) & mask
    value = value ^ (value >> 28)
    value = (value + (value << 31)) & mask
    return value


def serialize_ref_search_index(index: dict) -> dict:
    serialized = copy.deepcopy(index)
    serialized["minimizers"] = {str(key): value for key, value in serialized["minimizers"].items()}
    serialized["normalization"] = serialized["normalization"].tolist()
    return serialized


def format_fasta(identifier: str, description: str, sequence: str, line_width: int | None = None) -> str:
    if line_width:
        sequence = "\n".join(wrap(sequence, line_width))
    return f">{identifier} {description}\n{sequence}\n"


def first(qualifiers: dict[str, list[str]], key: str) -> str | None:
    values = qualifiers.get(key)
    return values[0] if values else None


def add_attribute(attributes: list[str], key: str, qualifiers: dict[str, list[str]]) -> None:
    for value in qualifiers.get(key, []):
        add_value(attributes, key, value)


def add_value(attributes: list[str], key: str, value: str | None) -> None:
    if value:
        attributes.append(f"{key}={escape_gff3(value)}")


def escape_gff3(value: str) -> str:
    return (
        value.replace("%", "%25")
        .replace(";", "%3B")
        .replace("=", "%3D")
        .replace("&", "%26")
        .replace("\t", "%09")
        .replace("\n", "%0A")
    )


def download_file(url: str, path: Path) -> None:
    with urllib.request.urlopen(url, timeout=60) as response:  # noqa: S310
        path.write_bytes(response.read())


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def write_json(path: Path, content: dict) -> None:
    write_text(path, json.dumps(content, indent=2) + "\n")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
