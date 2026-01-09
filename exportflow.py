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
from prefect.context import get_run_context
from datetime import datetime
import traceback
import json

with open("config.yaml", "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

EXPORT_DATA_DIRECTORY = config["exportflow"]["export_data_directory"]
PANDDA_ANALYSE_DIRECTORY = config["exportflow"]["pandda_analyse_directory"]

EXPORT_FILE_LIST = config["exportflow"]["export_list"]
with open(EXPORT_FILE_LIST, "r") as f:
    export_list = [line.strip() for line in f]

GROUND_BASENAME = config["exportflow"]["ground_basename"]
CHANGED_BASENAME = config["exportflow"]["changed_basename"]
ENSEMBLE_BASENAME = config["exportflow"]["ensemble_basename"]
VALIDATE_ONLY = config["exportflow"]["validate_only"]

if Path(PANDDA_ANALYSE_DIRECTORY).parts[:6] != Path(EXPORT_DATA_DIRECTORY).parts[:6]:
    raise Exception("Double check export and analyse directories, they should match.")

dir_dict = {
    sample: {
        "source_dir": Path(PANDDA_ANALYSE_DIRECTORY) / Path(sample),
        "export_dir": Path(EXPORT_DATA_DIRECTORY) / Path(sample),
        "ensemble": Path(EXPORT_DATA_DIRECTORY)
        / Path(f"{sample}/{sample}-{ENSEMBLE_BASENAME}.pdb"),
        "restraints": Path(EXPORT_DATA_DIRECTORY)
        / Path(f"{sample}/{sample}-{ENSEMBLE_BASENAME}.ff-refmac.params"),
    }
    for sample in export_list
}

source_dir_list = [v["source_dir"].is_dir() for v in dir_dict.values()]

if not all(source_dir_list):
    missing_dirs = [
        v["source_dir"] for v in dir_dict.values() if not v["source_dir"].is_dir()
    ]
    raise Exception(f"no analyse directories found for {missing_dirs}")

for sample, dirs in dir_dict.items():
    d = dirs["source_dir"]

    event_maps = list(d.glob("*event*.ccp4"))
    average_map = list(d.glob("*average*.ccp4"))
    ligand_files = list(d.glob("ligand_files/*.cif")) + list(
        d.glob("ligand_files/*.pdb")
    )

    # symlinks
    ground_state = list(d.glob(f"*{GROUND_BASENAME}*.pdb"))
    ground_state_reflections = list(d.glob(f"*{GROUND_BASENAME}*.mtz"))
    changed_state = list(d.glob(f"modelled_structures/*{CHANGED_BASENAME}*.pdb"))

    dir_dict[sample].update(
        {
            "xtal_id": sample,
            "event_maps": event_maps,
            "average_map": average_map,
            "ligand_files": ligand_files,
            "ground_state": ground_state,
            "ground_state_reflections": ground_state_reflections,
            "changed_state": changed_state,
            "messages": [],

        }
    )

def record_message(
        dir_dict: dict,
        stage: str,
        level: str,
        message: str,
        exc: Exception=None,
):
    entry = {
        'time': datetime.utcnow().isoformat(),
        'stage': stage,
        'level': level,
        'message': message,
    }
    if exc is not None:
        entry['exception'] = repr(exc)
        entry['traceback'] = traceback.format_exc()
    
    dir_dict['messages'].append(entry)


@task(name="validate", tags=["validate_job"])
def validate(dir_dict: dict, lig_label: str='UNL'):
    logger = get_run_logger()
    try:

        print(dir_dict["ground_state"])
        ground = gemmi.read_structure(str(dir_dict["ground_state"][0]))
        changed = gemmi.read_structure(str(dir_dict["changed_state"][0]))

        # chain mismatch check
        ground_chains = set([c.name for m in ground for c in m])
        changed_chains = set([c.name for m in changed for c in m])

        if ground_chains != changed_chains:
            raise Exception(f"Mismatching chain names between donor/acceptor states for {dir_dict}")

        # check for clash with placed fragments/ligands
        ns = gemmi.NeighborSearch(changed, 2.2, 0) # within 2.2A of model 0
        for nc, chain in enumerate(changed[0]):
            for nr, res in enumerate(chain):
                for na, atom in enumerate(res):
                    if res.name != lig_label:
                        ns.add_atom(atom, nc, nr, na)
        lig_atoms = [a for c in changed[0] for r in c for a in r if r.name == lig_label]
        for atom in lig_atoms:
            neighbors = ns.find_atoms(atom.pos)
            if neighbors:
                raise Exception(f'found clashing atoms near {atom} for {dir_dict}')


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
            raise Exception(f"Detected aniso adp in ground state for {dir_dict}")

        if np.mean(np.array([a.b_iso for a in ground_atoms])) < 5:
            raise Exception(f"suspiciously low B-factor for {dir_dict}")
        
        record_message(
            dir_dict,
            stage=get_run_context().task_run.name,
            level='info',
            message='validation passed'
        )
        dir_dict['validate_ok'] = True

    except Exception as e:
        print(f"CAUGHT VALIDATION EXCEPTION {e} FOR {dir_dict['xtal_id']}")
        record_message(
            dir_dict,
            stage=get_run_context().task_run.name,
            level='error',
            message=str(e),
            exc=e,
        )
        dir_dict['validate_ok'] = False

    return dir_dict


@task(name="make_export_dir", tags=["make_export_dir_job"])
def make_export_dir(dir_dict: dict):
    Path(dir_dict["export_dir"]).mkdir(parents=True, exist_ok=True)
    return dir_dict


@task(name="merge_ensemble", tags=["merge_ensemble_job"])
def merge_ensemble(dir_dict: dict, write_files=True):
    logger = get_run_logger()

    if dir_dict['validate_ok'] == False:
        record_message(
            dir_dict,
            stage=get_run_context().task_run.name,
            level='info',
            message='skipped due to validation fail'
        )
        return dir_dict

    try:
        ground = gemmi.read_structure(str(dir_dict["ground_state"][0]))
        changed = gemmi.read_structure(str(dir_dict["changed_state"][0]))
        em = ensemble_merge.EnsembleMerger(
            ground,
            changed,
            dir_dict["xtal_id"],
            occupancy_kwargs={"eps": 3.3, "min_samples": 1},
        )
        em.run()

        if write_files:
            # write model/restraint files to avoid passing anything to next task
            em.acceptor.write_pdb(str(dir_dict["ensemble"]))
            with open(dir_dict["restraints"], "w") as f:
                f.write(em.refmac_restraints_to_string())

        record_message(
            dir_dict,
            stage=get_run_context().task_run.name,
            level='info',
            message='ensemble merge passed',
        )

        dir_dict['ensemble_merge_ok'] = True

    except Exception as e:
        print(f"ENSEMBLE MERGE TASK FAILED xtal_id {dir_dict['xtal_id']} {e}")
        record_message(
            dir_dict,
            stage=get_run_context().task_run.name,
            level='error',
            message=str(e),
            exc=e,
        )

        dir_dict['ensemble_merge_ok'] = False

    return dir_dict


@task(name="copy_files", tags=["copy_files_job"])
def copy_files(dir_dict: dict):

    if dir_dict['ensemble_merge_ok']:
        s = dir_dict["source_dir"]
        d = dir_dict["export_dir"]
        files = (
            dir_dict["ground_state"]
            + dir_dict["changed_state"]
            + dir_dict["event_maps"]
            + dir_dict["ligand_files"]
            + dir_dict["ground_state_reflections"]
            + dir_dict["average_map"]
        )
        for f in files:
            shutil.copy2(f, d, follow_symlinks=True)

        return dir_dict
    



@flow(name="export_flow", task_runner=ConcurrentTaskRunner)
def export_flow(jobs, **kwargs):

    if VALIDATE_ONLY:
        validate_output = validate.map(jobs)
        ensemble_output = merge_ensemble.map(validate_output, write_files=False)
        return ensemble_output
    else:
        validate_output = validate.map(jobs)
        make_dir_output = make_export_dir.map(validate_output)
        merge_output = merge_ensemble.map(make_dir_output)
        copy_output = copy_files.map(merge_output)
        return copy_output


if __name__ == "__main__":
    n_cpus = multiprocessing.cpu_count()
    if n_cpus < 30:
        n_chunks = n_cpus - 2
    else:
        n_chunks = 30
    jobs_list = list(dir_dict.values())
    job_chunks = [
        jobs_list[i : i + n_chunks] for i in range(0, len(jobs_list), n_chunks)
    ]

    all_results = []
    for chunk in job_chunks:
        chunk_results = export_flow(chunk)
        all_results.extend([r.result() for r in chunk_results])

    with open("validate.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
