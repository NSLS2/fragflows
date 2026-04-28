from pathlib import Path
import pandas as pd
from fragflows_db.crud import mtz_from_xml
from fragflows_db.utils import sha256sum

# some helpful functions for combining dataframes from different
# fragflows versions



# 20260421
# helper function for extending the original csv/dataframe used for 
# dimpleflow with ispyb xml and hdf5 info needed for group deposition.

def merge_dfs_for_group_dep(df1: pd.DataFrame, df2: pd.DataFrame, columns_to_add: list = None) -> pd.DataFrame:
    if columns_to_add is None:
        columns_to_add = ["xml_path", "data_collection_date", "wavelength", "det_description", "det_serial_no"]
    
    for idx, row in df1.iterrows():
        mtz_filepath1 = row["filepath"]
        if not Path(mtz_filepath1).exists():
            raise Exception(f"File {mtz_filepath1} does not exist.")
        mtz_checksum = sha256sum(mtz_filepath1)
        print(f"Computed checksum for {mtz_filepath1}: {mtz_checksum}")
        df1.at[idx, "mtz_checksum"] = mtz_checksum

    for idx, row in df2.iterrows():
        df2.at[idx, "mtz_checksum"] = sha256sum(mtz_from_xml(row["xml_path"]))

    for col in columns_to_add:
        mapping = df2.dropna(subset=["mtz_checksum", col]).set_index("mtz_checksum")[col]
        df1[col] = df1["mtz_checksum"].map(mapping)
    
    return df1.copy()