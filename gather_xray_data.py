#!/nsls2/conda/envs/2023-1.1-py39/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 09:00:15 2023

@author: dkreitler
"""

import os
import pandas
import subprocess
import datetime

"""
INSTRUCTIONS:

Edit DATA_DIRECTORY to include path to xray data, e.g.
/nsls2/data/amx/proposals/2023-3/pass-313937/313937-20231105-dtime

Edit SAMPLE_NAME to include unique sample id in master.h5 filename
for including in the dataset
"""

with open("config.yaml","r") as yaml_file:
    config = yaml.safe_load(yaml_file)

DATA_DIRECTORY = config['gather_xray_data']['data_directory']
SAMPLE_NAME = config['gather_xray_data']['sample_name']

df = pandas.DataFrame()
find_cmd = ["find", f"{DATA_DIRECTORY}", "-maxdepth", "4", "-name", "*summary.csv"]
find_output = subprocess.check_output(find_cmd, universal_newlines=True)
csv_files = find_output.splitlines()

df = pandas.DataFrame()

for f in csv_files:
    df_to_add = pandas.read_csv(f)
    df_to_add["pipeline"] = ""
    for index, row in df_to_add.iterrows():
        df_to_add.at[index, "pipeline"] = f.split("/")[-1].split(".")[0]
    df = pandas.concat([df, df_to_add], ignore_index=True)

# fix autoPROC -> autoProc
df["pipeline"][df["pipeline"] == "autoPROC"] = "autoProc"

# add xtal_id
df["xtal_id"] = ""
for index, row in df.iterrows():
    df.at[index, "xtal_id"] = row["Sample_Path"]

df_filtered = df[df["xtal_id"].str.contains(f"{SAMPLE_NAME}")]
final_df = df_filtered.groupby("xtal_id").apply(
    lambda group: group.loc[group["Hi"].idxmin()]
)

# get list of all reflection files
find_cmd = [
    "find",
    f"{DATA_DIRECTORY}",
    "-name",
    "truncate-unique.mtz",
]
reflection_files = subprocess.check_output(
    find_cmd, universal_newlines=True
).splitlines()

print(reflection_files)

# find file paths
final_df["filepath"] = ""
for index, row in final_df.iterrows():
    for f in reflection_files:
        if (row["Sample_Path"] in f):
            final_df.at[index, "filepath"] = f

print(final_df)
final_df.to_csv(
    f"{SAMPLE_NAME}.{datetime.datetime.now().strftime('%Y%m%d')}.filtered.csv"
)
