# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import subprocess
import csv
from prefect import task, flow
from prefect.task_runners import ConcurrentTaskRunner
from pathlib import Path
import pandas
import multiprocessing
import shlex
import yaml
from datetime import datetime
import gemmi
import time
import argparse
from typing import Optional

with open("config.yaml", "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

EXPORT_DATA_DIRECTORY = config["refineflow"]["export_data_directory"]
LIGAND_CSV = config["refineflow"]["ligand_csv"]
ligand_df = pandas.read_csv(LIGAND_CSV)
BASENAME = config["refineflow"]["basename"]
BASENAME_HKL = config["refineflow"]["basename_hkl"]
REFINEMENT_PROGRAM = config["refineflow"]["refinement_program"]


def _parse_selected_datasets(arg_str: Optional[str]):
    if not arg_str:
        return None
    # support both comma-separated and accidental spaces
    items = [x.strip() for x in arg_str.split(",")]
    items = [x for x in items if x]
    return set(items) if items else set()


def build_jobs_list(selected: Optional[set[str]])->list:

    jobs_list = []

    for d in os.listdir(EXPORT_DATA_DIRECTORY):
        ligand_row = ligand_df.loc[ligand_df["xtal_id"] == d, "catalog_id"]
        print(ligand_row)
        ligand = ligand_row.iloc[0] if not ligand_row.empty else None
        print(ligand)
        if ligand == "DMSO":
            ligand = None

        job_dict = {
            "xyzin": f"{d}-{BASENAME}.pdb",
            "xyzout": f"{d}-{BASENAME}_refine",
            "hklin": f"{d}-{BASENAME_HKL}.mtz",
            "hklout": f"{d}-{BASENAME}_refine.mtz",
            "hklout_cif": f"{d}-{BASENAME}_refine.reflections.cif",
            "restraints": f"{d}-{BASENAME}.restraints-{REFINEMENT_PROGRAM}.params",
            "ligand": None,
            "sample_dir": str(Path(EXPORT_DATA_DIRECTORY) / Path(d)),
        }

        if ligand is not None:
            job_dict["ligand"] = f"{ligand}.cif"

        if selected:
            if d in selected:
                jobs_list.append(job_dict)
        else:
            jobs_list.append(job_dict)
    return jobs_list


@task(name="run_phenix", tags=["phenix_job"])
def run_phenix(phenix_params: dict):
    cmd = 'phenix.refine {xyzin} {hklin} {ligand} {restraints} write_reflection_cif_file=True strategy="*individual_sites *individual_adp *occupancies" --overwrite'.format(
        **phenix_params
    )
    print(f"running: {cmd}\nin {phenix_params['sample_dir']}")
    phenix_process = subprocess.Popen(
        shlex.split(cmd),
        cwd=phenix_params["sample_dir"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    phenix_process.communicate()


@task(name="run_refmac", tags=["refmac_job"])
def run_refmac(refine_params: dict):
    with open(
        f"{refine_params['sample_dir']}/refmac.{datetime.now().strftime('%Y%m%d%H%M%S')}.log",
        "w",
    ) as log_file:
        cmd = "refmac5 hklin {hklin} hklout {hklout} xyzin {xyzin} xyzout {xyzout} libin {ligand}".format(
            **refine_params
        )
        with subprocess.Popen(
            cmd.split(),
            stdin=subprocess.PIPE,
            stdout=log_file,
            cwd=refine_params["sample_dir"],
        ) as proc:
            print(f"running: {cmd}\nin {refine_params['sample_dir']}")
            # include custom occupancy group restraints
            # additional refmac config parameters for fine tuning
            stdout = proc.communicate(
                input=f"@{refine_params['restraints']}\n@{os.getcwd()}/refmac.params".encode()
            )
    return refine_params


@task(name="mtz2cif", tags=["mtz2cif_job"])
def mtz2cif(refine_params: dict):
    print(
        "converting {sample_dir}/{hklout} to {sample_dir}/{hklout_cif}".format(
            **refine_params
        )
    )
    try:
        mtz = gemmi.read_mtz_file("{sample_dir}/{hklout}".format(**refine_params))
    except RuntimeError as e:
        print(
            f"Caught {e} for {refine_params['hklout']}, waiting 5 sec and re-trying..."
        )
        time.sleep(5)
        mtz = gemmi.read_mtz_file("{sample_dir}/{hklout}".format(**refine_params))

    mtz_to_cif = gemmi.MtzToCif()
    doc = gemmi.cif.Document()
    doc.parse_string(mtz_to_cif.write_cif_to_string(mtz))
    for block in doc:
        if block.name == "mtz":
            block.name = refine_params["xyzin"]
    doc.write_file("{sample_dir}/{hklout_cif}".format(**refine_params))
    return refine_params


@task(name="refmac_mmcif_block_rename", tags=["refmac_mmcif_block_rename_job"])
def refmac_mmcif_block_rename(refine_params: dict):
    print(
        "renaming {sample_dir}/{xyzout}.mmcif cif block to denote ensembled".format(
            **refine_params
        )
    )
    doc = gemmi.cif.read_file("{sample_dir}/{xyzout}.mmcif".format(**refine_params))
    for block in doc:
        if block.name == "xxxx":  # default given by refmac
            block.name = refine_params["xyzout"]
            break  # assume that it is the first block
    # overwrite refmac-generated file with our updated version
    doc.write_file("{sample_dir}/{xyzout}.mmcif".format(**refine_params))


@flow(name="phenix_flow", task_runner=ConcurrentTaskRunner)
def phenix_flow(jobs, **kwargs):
    run_phenix.map(jobs)


@flow(name="refmac_flow", task_runner=ConcurrentTaskRunner)
def refmac_flow(jobs, **kwargs):
    output1 = run_refmac.map(jobs)
    output2 = mtz2cif.map(output1)
    refmac_mmcif_block_rename.map(output2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Refinement runner")
    parser.add_argument(
        "--select_datasets", "--select-datasets",
        dest="select_datasets",
        default=None,
        help='Comma-separated dataset IDs to refine, e.g. --select_datasets "xtal-001, xtal-002". '
             'If omitted, refine all datasets found.'
    )
    args, _ = parser.parse_known_args()
    selected = _parse_selected_datasets(args.select_datasets)
    jobs_list = build_jobs_list(selected)
    n_cpus = multiprocessing.cpu_count()
    if n_cpus < 30:
        n_chunks = n_cpus
    else:
        n_chunks = 30

    job_chunks = [
        jobs_list[i : i + n_chunks] for i in range(0, len(jobs_list), n_chunks)
    ]
    for chunk in job_chunks:
        if REFINEMENT_PROGRAM == "phenix":
            phenix_flow(chunk)
        elif REFINEMENT_PROGRAM == "refmac":
            refmac_flow(chunk)
        else:
            raise Exception(f"Refinement program {REFINEMENT_PROGRAM} not implemented")
