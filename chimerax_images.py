import gc
from chimerax.core.commands import run as cx
import sys

sys.path.append(
    "/nsls2/software/mx/conda_envs/2024-2.0-py311/lib/python3.11/site-packages"
)
import argparse
import gemmi
import tempfile
import shutil
import os
import numpy as np
from PIL import ImageDraw, ImageFont, Image
import time
import json
import re
import yaml


def get_events(sf):
    doc = gemmi.cif.read_file(sf)
    events = []
    for block, rblock in zip(doc, gemmi.as_refln_blocks(gemmi.cif.read_file(sf))):
        try:
            details = block.find_value("_diffrn.details")
            if details is None:
                continue

            s = gemmi.cif.as_string(details).strip()
            # _diffrn.details may be stored with an extra wrapping quote layer.
            if len(s) >= 2 and s[0] == s[-1] and s[0] in ("'", '"'):
                s = s[1:-1]

            parsed_s = json.loads(s)
            if "event_map" in str(parsed_s.get("ligand_evidence", "")):
                events.append(
                    {
                        "event": rblock,
                        "x": float(parsed_s["event_site_x"]),
                        "y": float(parsed_s["event_site_y"]),
                        "z": float(parsed_s["event_site_z"]),
                    }
                )

        except Exception as e:
            print(f"caught {e}")

    return events


# find list of ligands


def get_ligands(st, events):
    ligands = []
    for m in gemmi.read_structure(st):
        for c in m:
            for r in c:
                if r.name == "UNL":
                    lig_position = np.add.reduce([a.pos for a in r]) / len(r)
                    lig_position = np.array(lig_position.tolist())
                    dist_to_events = [
                        lig_position
                        - np.array(
                            [event["x"], event["y"], event["z"]], dtype=np.float64
                        )
                        for event in events
                    ]
                    min_ind, _ = min(
                        enumerate(dist_to_events), key=lambda x: np.linalg.norm(x[1])
                    )
                    ligands.append(
                        {
                            "chain": c.name,
                            "seqid": r.seqid.num,
                            "position": lig_position,
                            "event": events[min_ind],
                        }
                    )

    print(ligands)
    return ligands


def add_text_to_img(image: str, text: str):

    # Open your PNG
    img = Image.open(image)

    # Create a drawing object
    draw = ImageDraw.Draw(img)

    # Choose font (defaults to a basic one if you don’t have a ttf file)
    # You can specify a font file (like Arial.ttf) if you want nicer text
    font = ImageFont.load_default()

    draw.multiline_text(
        (20, 20),
        text,
        fill="red",
        font=ImageFont.truetype(
            "/usr/share/fonts/nimbus-sans-l/NimbusSanL-Bold.ttf", 16
        ),
    )

    # Save result
    img.save(image)




def write_image(ligand: dict, dataset: str, group_dir: str, out_dir=os.getcwd()):
    with tempfile.TemporaryDirectory() as tmpdir:
        st = f"{group_dir}/{dataset}.cif"
        sf = f"{group_dir}/{dataset}-sf.cif"
        doc = gemmi.cif.read_file(f"{sf}")
        rblocks = gemmi.as_refln_blocks(doc)
        g = ligand["event"]["event"].transform_f_phi_to_map("pdbx_FWT", "pdbx_PHWT")
        print(g)
        cm = gemmi.Ccp4Map()
        cm.grid = g
        cm.update_ccp4_header()
        cm.write_ccp4_map(f"{tmpdir}/fwt.ccp4")
        del g
        del cm
        gc.collect()
        cx(session, f"open {st}")
        cx(session, "hide #1 atoms")
        cx(session, f"show /{ligand['chain']}:{ligand['seqid']} atoms")
        cx(session, f"select /{ligand['chain']}:{ligand['seqid']}")
        cx(session, "view sel")
        cx(session, f"color /{ligand['chain']}:{ligand['seqid']} byelement")
        cx(session, "set bgColor white")
        cx(session, f"open {tmpdir}/fwt.ccp4")
        cx(session, f"vop cover #1 #2 atomBox sel")
        cx(session, "volume #3 level 0.25")
        cx(session, "color #3 blue")
        cx(session, "transparency 50 surfaces")
        cx(session, "volume zone #3 near sel newMap false")
        cx(session, "close #2")
        cx(session, "hide cartoons")
        img_name = f"{dataset}_{ligand['chain']}{ligand['seqid']}.png"
        cx(
            session,
            f"save {out_dir}/{img_name} transparentBackground false width 512 height 512",
        )
        add_text_to_img(
            f"{out_dir}/{img_name}",
            f"{dataset}\nevent pos.: {ligand['event']['x']},{ligand['event']['y']},{ligand['event']['z']}\nlig. pos.: {ligand['position']}\n{ligand['chain']}{ligand['seqid']}",
        )
        cx(session, "turn y 90")
        img_name_rot90 = f"{dataset}_{ligand['chain']}{ligand['seqid']}_rot90.png"
        cx(
            session,
            f"save {out_dir}/{img_name_rot90} transparentBackground false width 512 height 512",
        )
        add_text_to_img(
            f"{out_dir}/{img_name_rot90}",
            f"{dataset}\nevent pos.: {ligand['event']['x']},{ligand['event']['y']},{ligand['event']['z']}\nlig. pos.: {ligand['position']}\n{ligand['chain']}{ligand['seqid']}\nrotated 90 deg",
        )
        # shutil.move("test4.png", cwd)
        cx(session, "close #3")
        cx(session, "close #1")


def write_images(dataset, group_dir, out_dir=os.getcwd()):
    st = f"{group_dir}/{dataset}.cif"
    sf = f"{group_dir}/{dataset}-sf.cif"
    events = get_events(sf)
    ligands = get_ligands(st, events)
    for ligand in ligands:
        write_image(ligand, dataset, group_dir, out_dir=out_dir)


parser = argparse.ArgumentParser(description="Generate ChimeraX ligand images from group deposition CIFs")
parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
parser.add_argument("--out-dir", default=os.getcwd(), help="Output directory for images")
parser.add_argument("--exclude-pattern", default="ground", help="Regex pattern to exclude datasets by name")
args, _ = parser.parse_known_args()

with open(args.config) as f:
    config = yaml.safe_load(f)

group_dir = config["groupdepflow"]["groupdep_directory"]
out_dir = args.out_dir
os.makedirs(out_dir, exist_ok=True)

exclude_re = re.compile(args.exclude_pattern)
datasets = [x[:-7] for x in os.listdir(group_dir) if x.endswith("sf.cif")]
print(datasets)
for d in datasets:
    if exclude_re.search(d):
        print(f"Skipping {d} (matches exclude pattern '{args.exclude_pattern}')")
        continue
    print(d)
    write_images(d, group_dir, out_dir=out_dir)
cx(session, "exit")
