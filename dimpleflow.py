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
import yaml

with open("config.yaml","r") as yaml_file:
    config = yaml.safe_load(yaml_file)

DATA_DIRECTORY = config['gather_xray_data']['data_directory']
SAMPLE_NAME = config['gather_xray_data']['sample_name']
PROCESSING_DATA_DIRECTORY = ['dimpleflow']['processing_data_directory']
MODELS_DIRECTORY = ['dimpleflow']['models_directory']
REFERENCE_PDB = ['dimpleflow']['reference_pdb']
FILTERED_XRAY_CSV = ['dimpleflow']['filtered_xray_csv']

jobs_df = pandas.read_csv(FILTERED_XRAY_CSV)
jobs_list = []

for index, row in jobs_df.iterrows():
    jobs_list.append(
        {
            "hklout": f"{row['xtal_id']}.dimple.mtz",
            "xyzout": f"{row['xtal_id']}.dimple.pdb",
            "xyzin": reference_model,
            "hklin": row["filepath"],
            "sample_dir": str(models_dir / Path(f"{row['xtal_id']}")),
            "xtal_id": row["xtal_id"],
        }
    )


@task(name="run_dimple", tags=["dimple_job"])
def run_dimple(dimple_params: dict):
    cmd = "dimple --hklout {hklout} --xyzout {xyzout} {xyzin} {hklin} {sample_dir}".format(
        **dimple_params
    )
    print(cmd)
    dimple_process = subprocess.Popen(
        cmd.split(),
        cwd=models_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    dimple_process.communicate()


@flow(name="dimple_flow", task_runner=ConcurrentTaskRunner)
def dimple_flow(jobs, **kwargs):
    run_dimple.map(jobs)


if __name__ == "__main__":
    n_cpus = multiprocessing.cpu_count()
    if n_cpus < 30:
        n_chunks = n_cpus
    else:
        n_chunks = 30

    job_chunks = [jobs_list[i : i + n_chunks] for i in range(0, len(jobs_list), n_chunks)]
    for chunk in job_chunks:
        dimple_flow(chunk)
