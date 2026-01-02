# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import shutil
from prefect import task, flow, get_run_logger
from prefect.task_runners import ConcurrentTaskRunner
from pathlib import Path
import multiprocessing
import yaml
from ensemble_models import sync_solvent, ensemble_merge
import gemmi
import numpy as np

with open("config.yaml", "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

EXPORT_DATA_DIRECTORY = config["exportflow"]["export_data_directory"]
PANDDA_ANALYSE_DIRECTORY = config["exportflow"]["pandda_analyse_directory"]

EXPORT_FILE_LIST = config["exportflow"]["export_list"]
with open(EXPORT_FILE_LIST, 'r') as f:
    export_list = [line.strip() for line in f]

GROUND_BASENAME = config["exportflow"]["ground_basename"]
CHANGED_BASENAME = config["exportflow"]["changed_basename"]
ENSEMBLE_BASENAME = config["exportflow"]["ensemble_basename"]
VALIDATE_ONLY = config["exportflow"]["validate_only"]

if Path(PANDDA_ANALYSE_DIRECTORY).parts[:6] != Path(EXPORT_DATA_DIRECTORY).parts[:6]:
    raise Exception('Double check export and analyse directories, they should match.')

file_pattern_manifest = [
    '*export*.ccp4',
    f'*{GROUND_BASENAME}*.pdb',
    f'*{CHANGED_BASENAME}*.pdb',
    'ligand_files/*.cif',
    'ligand_files/*.pdb',
]

dir_dict = {
    sample: {
        'source_dir': Path(PANDDA_ANALYSE_DIRECTORY) / Path(sample),
        'export_dir': Path(EXPORT_DATA_DIRECTORY) / Path(sample),
        'ensemble': Path(EXPORT_DATA_DIRECTORY) / Path(f'{sample}/{sample}-{ENSEMBLE_BASENAME}.pdb'),
        'restraints': Path(EXPORT_DATA_DIRECTORY) / Path(f'{sample}/{sample}-{ENSEMBLE_BASENAME}.ff-refmac.params')
    }
    for sample in export_list
}

source_dir_list = [v['source_dir'].is_dir() for v in dir_dict.values()]

if not all(source_dir_list):
    missing_dirs = [v['source_dir'] for v in dir_dict.values() if not v['source_dir'].is_dir()]
    raise Exception(f'no analyse directories found for {missing_dirs}')

for sample, dirs in dir_dict.items():
    d = dirs['source_dir']

    event_maps = list(d.glob('*event*.ccp4'))
    ligand_files = list(d.glob('ligand_files/*.cif')) + list(d.glob('ligand_files/*.pdb'))
    
    # symlinks
    ground_state = list(d.glob(f'*{GROUND_BASENAME}*.pdb'))
    ground_state_reflections = list(d.glob(f'*{GROUND_BASENAME}*.mtz'))
    changed_state = list(d.glob(f'modelled_structures/*{CHANGED_BASENAME}*.pdb'))

    dir_dict[sample].update({
        'xtal_id': sample,
        'event_maps': event_maps,
        'ligand_files': ligand_files,
        'ground_state': ground_state,
        'ground_state_reflections': ground_state_reflections,
        'changed_state': changed_state,
    })

import pprint
pprint.pprint(dir_dict)


@task(name="validate", tags=["validate_job"])
def validate(dir_dict: dict):
    logger = get_run_logger()
    try:

        print(dir_dict['ground_state'])
        ground = gemmi.read_structure(str(dir_dict['ground_state'][0]))
        changed = gemmi.read_structure(str(dir_dict['changed_state'][0]))

        # water clash
        sync_solvent.check_for_solvent_clash(ground)
        sync_solvent.check_for_solvent_clash(changed)

        # ion/lig/nonpolymer clash
        sync_solvent.check_for_nonpolymer_clashes(ground)
        sync_solvent.check_for_nonpolymer_clashes(changed)

        ground_atoms = [a for m in ground for c in m for r in c for a in r]
        changed_atoms = [a for m in changed for c in m for r in c for a in r]

        # check for b_aniso
        if any([b_eq > 0 for b_eq in [b.b_eq() for b in ground_atoms]]):
            raise Exception(f'Detected aniso adp in ground state for {dir_dict}')
        
        if np.mean(np.array([a.b_iso for a in ground_atoms])) < 5:
            raise Exception(f'suspiciously low B-factor for {dir_dict}')
        return dir_dict
    
    except Exception:
        logger.exception(
            f"VALIDATE TASK FAILED xtal_id: {dir_dict['xtal_id']}"
        )
        raise

@task(name="make_export_dir", tags=["make_export_dir_job"])
def make_export_dir(dir_dict: dict):
    Path(dir_dict['export_dir']).mkdir(parents=True, exist_ok=True)
    return dir_dict

@task(name="merge_ensemble", tags=["merge_ensemble_job"])
def merge_ensemble_dir(dir_dict: dict):
    logger = get_run_logger()
    try:
        ground = gemmi.read_structure(str(dir_dict['ground_state'][0]))
        changed = gemmi.read_structure(str(dir_dict['changed_state'][0]))
        em = ensemble_merge.EnsembleMerger(
            ground,
            changed,
            dir_dict['xtal_id'],
            occupancy_kwargs={"eps": 3.3, "min_samples":1}
        )
        em.run()

        # write model/restraint files to avoid passing anything to next task
        em.acceptor.write_pdb(str(dir_dict['ensemble']))
        with open(dir_dict['restraints'], 'w') as f:
            f.write(em.refmac_restraints_to_string())

        return dir_dict

    except Exception:
        logger.exception(
            f"ENSEMBLE MERGE TASK FAILED xtal_id: {dir_dict['xtal_id']}"
        )
        raise

@task(name="copy_files", tags=["copy_files_job"])
def copy_files(dir_dict: dict):
    s = dir_dict['source_dir']
    d = dir_dict['export_dir']
    files = dir_dict['ground_state'] + dir_dict['changed_state'] + dir_dict['event_maps'] + dir_dict['ligand_files'] + dir_dict['ground_state_reflections']
    for f in files:
        shutil.copy2(f, d, follow_symlinks=True)

    return dir_dict

@flow(name="export_flow", task_runner=ConcurrentTaskRunner)
def export_flow(jobs, **kwargs):
    validate_output = validate.map(jobs)
    make_dir_output = make_export_dir.map(validate_output)
    merge_output = merge_ensemble_dir.map(make_dir_output)
    copy_output = copy_files.map(merge_output)


if __name__ == "__main__":
    n_cpus = multiprocessing.cpu_count()
    if n_cpus < 30:
        n_chunks = n_cpus - 2
    else:
        n_chunks = 5
    jobs_list = list(dir_dict.values())
    job_chunks = [
        jobs_list[i : i + n_chunks] for i in range(0, len(jobs_list), n_chunks)
    ]

    for chunk in job_chunks:
        export_flow(chunk)