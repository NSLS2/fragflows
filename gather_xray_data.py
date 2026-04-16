#!/nsls2/conda/envs/2023-1.1-py39/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 09:00:15 2023

@author: dkreitler
"""

import pandas
import datetime
import yaml
import numpy as np
import argparse
import os
from fragflows_db.ingest import ingest_mx_processing_path
from fragflows_db.utils import get_xml_paths
from fragflows_db.database_init import init_db, session_scope
from fragflows_db.crud import get_mx_processing_results_df, extend_dataframe_mx_stats, reflection_file_to_df
import re

with open("config.yaml", "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

DATA_DIRECTORIES = config["gather_xray_data"]["data_directory"]
SAMPLE_NAME = config["gather_xray_data"]["sample_name"]

parser = argparse.ArgumentParser()
parser.add_argument(
    "--update",
    action="store_true",
    help="if set, will update the database with new files instead of creating a new one",
)
parser.add_argument(
    "--db",
    type=str,
    default="fragflows.db",
    help="path to SQLite database file",
)
parser.add_argument(
    "--r_mrg_threshold",
    type=float,
    default=0.45,
)
parser.add_argument(
    "--symm",
    type=str,
    help="a fuzzy selection filter for space group, e.g. C222 will pick up C2221"
)
parser.add_argument(
    "--symm_strict",
    type=str,
    help="if set, will only consider files with this exact space group."
)
args = parser.parse_args()
init_db(args.db)

if args.update:
    if type(DATA_DIRECTORIES) != list:
        DATA_DIRECTORIES = [DATA_DIRECTORIES]

    xml_paths = []
    for data_dir in DATA_DIRECTORIES:
        if not os.path.exists(data_dir):
            raise Exception(f"data directory {data_dir} does not exist")
        
        xml_paths.extend(get_xml_paths(data_dir))

    xml_paths = [xp for xp in xml_paths if SAMPLE_NAME in xp]
    print(f"found {len(xml_paths)} xml files in {DATA_DIRECTORIES} matching sample name {SAMPLE_NAME}")

    for xp in xml_paths:
        print(f"found {len(xp)} xml files in {DATA_DIRECTORIES}")
        with session_scope() as session:
            try:
                ingest_mx_processing_path(session, xp)
                print(f"successfully ingested {xp}")
            except Exception as e:
                print(f"failed to ingest {xp}: {e}")

with session_scope() as session:
    df = get_mx_processing_results_df(session)
    df = extend_dataframe_mx_stats(df)
    df = reflection_file_to_df(df)

def normalize_space_group(sg: str) -> str:
    """Normalize space group string for fuzzy matching."""
    return "".join(str(sg).split()).upper()

def fuzzy_pattern(symm: str) -> str:
    return '.*'.join(re.escape(ch) for ch in normalize_space_group(symm))

df["symm"] = df["symm"].fillna("").apply(normalize_space_group)

if args.symm and args.symm_strict:
    raise Exception("cannot use both symm and symm_strict arguments together")
elif args.symm:
    pattern = fuzzy_pattern(args.symm)
    df = df[df["symm"].str.contains(pattern, regex=True)]
elif args.symm_strict:
    df = df[df["symm"] == normalize_space_group(args.symm_strict)]

def pick_row(group: pandas.DataFrame)->pandas.Series:

    acceptable_r_mrg = group.loc[group["r_mrg"] <= args.r_mrg_threshold]
    if not acceptable_r_mrg.empty:
        return group.loc[acceptable_r_mrg["hi"].idxmin()]

    less_acceptable_r_mrg = group.loc[group["r_mrg"] > args.r_mrg_threshold]
    if not less_acceptable_r_mrg.empty:
        return group.loc[less_acceptable_r_mrg["r_mrg"].idxmin()]

    return group.iloc[0]

df["r_mrg"] = df["r_mrg"].astype(np.float64)
df = df.groupby("xtal_id", group_keys=False).apply(pick_row, include_groups=False)

print(df)
df.to_csv(
    f"{SAMPLE_NAME}.{datetime.datetime.now().strftime('%Y%m%d')}.filtered.csv",
    index=False
)
