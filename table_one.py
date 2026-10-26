import pandas
import yaml
from deposition.load import cif_to_payload
from deposition.models import ReflectionStats, GroupDepCif
import os
import gemmi
import argparse

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

groupdep_directory = config["groupdepflow"]["groupdep_directory"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv",
        type=str,
        default="table_one.csv",
        help="Path to the output CSV file",
    )
    parser.add_argument(
        "--groupdep_directory",
        type=str,
        default=None,
        help="Path to the groupdep directory",
    )

    parser.add_argument(
        "--ground_model_id","--ground-model-id",
        type=str,
        default=None,
        help="ID of the ground model structure factor cif",
    )

    args = parser.parse_args()
    groupdep_directory = args.groupdep_directory if args.groupdep_directory else groupdep_directory
    ground_model_id = args.ground_model_id if args.ground_model_id else "ground-sf"


    models = []

    for file in os.listdir(groupdep_directory):

        # changed state cif files
        if file.endswith(".cif") and "-sf" not in file:
            cif_file_path = os.path.join(groupdep_directory, file)
            cif_block = gemmi.cif.read_file(cif_file_path)[0]
            payload = cif_to_payload(cif_block)
            rs = ReflectionStats(**payload)
            groupdep_cif = GroupDepCif(reflection_stats=rs, **payload)
            models.append(groupdep_cif)

        elif file.endswith("ground-sf.cif") or file.endswith(f"{ground_model_id}-sf.cif"):
            cif_file_path = os.path.join(groupdep_directory, file)
            for block in gemmi.cif.read_file(cif_file_path):
                payload = cif_to_payload(block)
                rs = ReflectionStats(**payload)
                groupdep_cif = GroupDepCif(reflection_stats=rs, **payload)
                models.append(groupdep_cif)

    df = pandas.DataFrame([model.model_dump_flat(mode="json") for model in models])
    # Drop columns where all entries are "?"
    df = df.loc[:, (df != "?").any()]
    df = df.drop(columns=['reflns_pdbx_diffrn_id','reflns_shell_pdbx_diffrn_id'])
    print(df)
    df.to_csv(args.csv, index=False)

        

