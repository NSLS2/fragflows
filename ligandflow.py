#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 14:13:19 2023

@author: dkreitler
"""

import os
import subprocess
from prefect import task, flow
from prefect.task_runners import ConcurrentTaskRunner
from pathlib import Path
import yaml
import argparse
import pandas as pd
from ligands.utils import enumerate_stereo_from_smiles

args = argparse.ArgumentParser()

args.add_argument(
    "--datasets",
    type=str,
    default=None,
    help="comma separated list of xtal_ids to process"
)

with open("config.yaml", "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

MODELS_DIRECTORY = config["ligandflow"]["models_directory"]
LIGAND_CSV = config["ligandflow"]["ligand_csv"]


@task(name="locate_models_dir", tags="acedrg_job")
def find_sample_path(sample_dict: dict):
    for r, d, f in os.walk(MODELS_DIRECTORY):
        for d_ in d:
            if sample_dict["xtal_id"] == d_:
                return Path(r, d_)  # sample_dir
    return None


@task(name="validate_sample_path", tags=["acedrg_job"])
def validate_sample_dir(sample_path: Path):
    sample_name = sample_path.parts()[-1]
    return all(
        [
            (sample_path / Path(f"{sample_name}{extension}")).exists()
            for extension in (".dimple.pdb", ".dimple.mtz")
        ]
    )


@task(name="generate_acedrg_params", tags=["acedrg_job"])
def generate_acedrg_params(sample_dict: dict):
    """do some processing on the sample_dict to make a dict that contains
    acedrg parameters"""
    smiles_list = sample_dict["smiles"].split(".")
    smiles_list.sort(key=len)
    smiles = smiles_list[-1]  # largest molecule

    acedrg_params = {
        "smiles": smiles,
        "catalog_id": sample_dict["catalog_id"],
    }
    
    return acedrg_params


@task(name="run_acedrg", tags=["acedrg_job"])
def run_acedrg(acedrg_params: dict, sample_path: Path):
    stereo_isomers = enumerate_stereo_from_smiles(acedrg_params["smiles"])
    if len(stereo_isomers) > 16:
        raise ValueError(f"Too many stereoisomers for SMILES: {acedrg_params['smiles']}", acedrg_params)

    isomer_count = 0
    for stereo_isomer in stereo_isomers:
        acedrg_params["smiles"] = stereo_isomer
        acedrg_params["isomer_count"] = str(isomer_count).zfill(2)
        isomer_count += 1
        cmd = "acedrg --smi {smiles} -o {catalog_id}_{isomer_count} -r {isomer_count}".format(**acedrg_params)
        acedrg_process = subprocess.Popen(
            cmd.split(),
            cwd=sample_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        acedrg_process.communicate()


@flow(name="acedrg_flow", task_runner=ConcurrentTaskRunner)
def acedrg_flow(samples: list, **kwargs):
    sample_paths = [find_sample_path(sample_dict) for sample_dict in samples]

    run_acedrg.map(
        [
            generate_acedrg_params(samples[idx])
            for idx, sample_dir in enumerate(sample_paths)
            if sample_dir != None
        ],
        [path for path in sample_paths if path != None],
    )


if __name__ == "__main__":
    parsed_args = args.parse_args()

    ligand_df = pd.read_csv(LIGAND_CSV)
    if parsed_args.datasets:
        dataset_list = [d.strip() for d in parsed_args.datasets.split(",") if d.strip()]
        ligand_df = ligand_df[ligand_df["xtal_id"].isin(dataset_list)]

    samples = ligand_df.to_dict("records")

    sample_chunks = [samples[i : i + 20] for i in range(0, len(samples), 20)]
    for chunk in sample_chunks:
        acedrg_flow(chunk)
