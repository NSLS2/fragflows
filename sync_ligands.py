#!/usr/bin/env python
"""Sync ligand PDB/CIF files from a models directory into matching
pandda_analyse processed_datasets/<sample>/ligand_files directories.

For each subdirectory name shared between the models directory and
pandda_analyse/processed_datasets, this clears out the corresponding
ligand_files directory and copies over any files matching
*_[0-9][0-9].pdb or *_[0-9][0-9].cif from the models subdirectory.
"""
import argparse
import shutil
from pathlib import Path

LIGAND_PATTERNS = ("*_[0-9][0-9].pdb", "*_[0-9][0-9].cif")


def find_ligand_files(sample_dir: Path):
    for pattern in LIGAND_PATTERNS:
        yield from sample_dir.glob(pattern)


def sync_ligands(models_dir: Path, pandda_analyse_dir: Path, dry_run: bool = False):
    processed_datasets_dir = pandda_analyse_dir / "processed_datasets"
    if not processed_datasets_dir.is_dir():
        raise FileNotFoundError(
            f"processed_datasets directory not found: {processed_datasets_dir}"
        )

    sample_names = sorted(
        p.name for p in models_dir.iterdir() if p.is_dir()
    )

    for sample_name in sample_names:
        model_sample_dir = models_dir / sample_name
        dataset_dir = processed_datasets_dir / sample_name
        if not dataset_dir.is_dir():
            continue

        ligand_files_dir = dataset_dir / "ligand_files"

        ligand_files = sorted(find_ligand_files(model_sample_dir))
        if not ligand_files:
            continue

        print(f"[{sample_name}] syncing {len(ligand_files)} file(s) -> {ligand_files_dir}")

        if dry_run:
            for existing in sorted(ligand_files_dir.glob("*")) if ligand_files_dir.is_dir() else []:
                print(f"  would remove {existing}")
            for src in ligand_files:
                print(f"  would copy {src} -> {ligand_files_dir / src.name}")
            continue

        # only reached when matching ligand files exist in model_sample_dir, so it's safe to clear
        if ligand_files_dir.is_dir():
            for existing in ligand_files_dir.iterdir():
                if existing.is_file() or existing.is_symlink():
                    existing.unlink()
                else:
                    shutil.rmtree(existing)
        else:
            ligand_files_dir.mkdir(parents=True)

        for src in ligand_files:
            shutil.copy2(src, ligand_files_dir / src.name)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("models_dir", type=Path, help="Directory containing per-sample model subdirectories")
    parser.add_argument(
        "pandda_analyse_dir",
        type=Path,
        help="pandda_analyse directory containing processed_datasets/<sample>/ligand_files",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be done without modifying anything",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    sync_ligands(args.models_dir, args.pandda_analyse_dir, dry_run=args.dry_run)


if __name__ == "__main__":
    main()