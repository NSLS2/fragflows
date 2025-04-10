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

with open("config.yaml","r") as yaml_file:
    config = yaml.safe_load(yaml_file)

EXPORT_DATA_DIRECTORY = config['refineflow']['export_data_directory']
LIGAND_CSV = config['refineflow']['ligand_csv']
ligand_df = pandas.read_csv(ligand_csv)
BASENAME = config['refineflow']['basename']
BASENAME_HKL = config['refineflow']['basename_hkl']

jobs_list = []

for d in os.listdir(EXPORT_DATA_DIRECTORY):
    ligand_row = ligand_df.loc[ligand_df['xtal_id'] == d, 'catalog_id']
    print(ligand_row)
    ligand = ligand_row.iloc[0] if not ligand_row.empty else None
    print(ligand)
    if ligand == "DMSO":
        ligand = None

    job_dict = {
            "xyzin": f"{d}-{BASENAME}.pdb",
            "hklin": f"{d}-{BASENAME_HKL}.mtz",
            "restraints": f"{d}-{BASENAME}.restraints-phenix.params",
            "ligand": None,
            "sample_dir": str(Path(EXPORT_DATA_DIRECTORY) / Path(d)),
        }
 
    if ligand is not None:
        job_dict['ligand'] = f"{ligand}.cif"
        
    jobs_list.append(job_dict)

@task(name="run_phenix", tags=["phenix_job"])
def run_phenix(phenix_params: dict):
    cmd = 'phenix.refine {xyzin} {hklin} {ligand} write_reflection_cif_file=True strategy="*individual_sites *individual_adp *occupancies" --overwrite'.format(
        **phenix_params
    )
    print(f"running: {cmd}\nin {phenix_params['sample_dir']}")
    phenix_process = subprocess.Popen(
        shlex.split(cmd),
        cwd=phenix_params['sample_dir'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    phenix_process.communicate()


@flow(name="phenix_flow", task_runner=ConcurrentTaskRunner)
def phenix_flow(jobs, **kwargs):
    run_phenix.map(jobs)


if __name__ == "__main__":
    n_cpus = multiprocessing.cpu_count()
    if n_cpus < 30:
        n_chunks = n_cpus
    else:
        n_chunks = 30

    job_chunks = [jobs_list[i : i + n_chunks] for i in range(0, len(jobs_list), n_chunks)]
    for chunk in job_chunks:
        phenix_flow(chunk)
