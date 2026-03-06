from coot import *
from coot_utils import *

import yaml
from pathlib import Path
import os
import argparse

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk

# This script when loaded in coot 1.1.xx will load
# autorefine structure factors, models, and event maps
# from an exportflow directory that has been autorefined
# with refmac5. A navigation browser gtk4 window will
# pop up when loaded. Tested with Coot 1.1.19.

def script_args(script_filename="coot_autorefine_inspect.py"):
    for idx, a in enumerate(sys.argv):
        if os.path.basename(a) == script_filename:
            return sys.argv[idx + 1:]
    
    if "--" in sys.argv:
        return sys.argv[sys.argv.index("--") + 1:]
    
    return []

args = script_args()

if len(args) > 1:
    raise Exception('too many input args')

start_n = 0
if len(args) == 1:
    start_n = int(args[0]) - 1


with open('config.yaml','r') as f:
    data = yaml.safe_load(f)

EXPORT_DIR = data['refineflow']['export_data_directory']
DATASETS = []

for dataset in os.listdir(EXPORT_DIR):
    DATASETS.append(
        (f"{EXPORT_DIR}/{dataset}/{dataset}-ensemble-model_refine.pdb",
         f"{EXPORT_DIR}/{dataset}/{dataset}-ensemble-model_refine.mtz",
         [str(d) for d in Path(f"{EXPORT_DIR}/{dataset}").glob("*event*.ccp4")],
         dataset)
    )

N_DATASETS = len(DATASETS)
print(DATASETS)


# MTZ columns labels assume refmac5
F_COL    = "FWT"
PHI_COL  = "PHWT"
DF_COL   = "DELFWT"
DPHI_COL = "PHDELWT"

# --- 2) State (current index + currently loaded molecule IDs) ---
_state = {
    "i": start_n,
    "imol_model": None,
    "imol_2fofc": None,
    "imol_fofc": None,
    "event_maps": [],
    "event_displayed_idx": None,
    "xtal_id": None
}

def _close_current():
    """Close currently loaded molecules/maps (optional but nice)."""
    for k in ("imol_model", "imol_2fofc", "imol_fofc"):
        imol = _state.get(k)
        if isinstance(imol, int) and imol >= 0:
            try:
                close_molecule(imol)          # mark closed
            except Exception:
                pass
            _state[k] = None

    event_maps = _state.get("event_maps")
    if event_maps:
        for imol in event_maps:
            if isinstance(imol, int) and imol >=0:
                try:
                    close_molecule(imol)
                except Exception:
                    pass
    _state["event_maps"] = []
    _state["event_displayed_idx"] = None
    _state["xtal_id"] = None

    # Some builds actually delete closed molecules only after this call:
    try:
        end_delete_closed_molecules()
    except Exception:
        pass

def _load_dataset(i):
    model_file, mtz_file, event_maps, xtal_id = DATASETS[i]

    _close_current()

    # Load model
    imol_model = read_pdb(model_file)

    # Load maps from MTZ
    imol_2fofc = make_and_draw_map(mtz_file, F_COL, PHI_COL, "", 0, 0)
    set_map_colour(imol_2fofc,0,0,0.9)

    # Load difference map from MTZ
    imol_fofc = make_and_draw_map(mtz_file, DF_COL, DPHI_COL, "", 0, 1)

    _state["imol_model"] = imol_model
    _state["imol_2fofc"] = imol_2fofc
    _state["imol_fofc"] = imol_fofc

    for event_map in event_maps:
        imol_event_map = read_ccp4_map(event_map, False)
        _state["event_maps"].append(imol_event_map)
        set_map_displayed(imol_event_map, 0)
        set_map_colour(imol_event_map, 1, 1, 1)

    _state["xtal_id"] = xtal_id

    # Optional: tell Coot which map to refine against (2Fo-Fc)
    try:
        set_imol_refinement_map(imol_2fofc)
    except Exception:
        pass

    add_status_bar_text(f"Loaded dataset {i+1}/{len(DATASETS)}")
    update_tool_window_labels()

def next_dataset():
    print('next pressed')
    if not DATASETS:
        add_status_bar_text("No datasets configured")
        return
    _state["i"] = (_state["i"] + 1) % len(DATASETS)
    _load_dataset(_state["i"])

def prev_dataset():
    if not DATASETS:
        add_status_bar_text("No datasets configured")
        return
    _state["i"] = (_state["i"] - 1) % len(DATASETS)
    _load_dataset(_state["i"])

def next_event():
    print('next event pressed')
    if _state["event_maps"]:
        print('next_event maps detected')
        if _state["event_displayed_idx"] is not None:
            imol = _state["event_maps"][_state["event_displayed_idx"]]
            set_map_displayed(imol, 0)
            _state["event_displayed_idx"] = (_state["event_displayed_idx"] + 1) % len(_state["event_maps"])
        else:
            _state["event_displayed_idx"] = 0
        imol = _state["event_maps"][_state["event_displayed_idx"]]
        set_map_displayed(imol, 1)
    else:
        print("no event maps to display")
        
def toggle_2fofc():
    print('toggle 2fofc pressed')
    imol = _state["imol_2fofc"]
    set_map_displayed(imol, 0 if map_is_displayed(imol) else 1)

def toggle_fofc():
    print('toggle fofc pressed')
    imol = _state["imol_fofc"]
    set_map_displayed(imol, 0 if map_is_displayed(imol) else 1)

def toggle_event():
    print('toggle event map pressed')
    if len(_state["event_maps"]) == 0:
        print("no event maps to toggle")
        return
    elif _state["event_displayed_idx"] is None:
        print("no event map displayed")
        return
    else:
        imol = _state["event_maps"][_state["event_displayed_idx"]]
        set_map_displayed(imol, 0 if map_is_displayed(imol) else 1)


# --- 3) Bind keys ---

def _next_cb(*args):
    add_status_bar_text("NEXT fired")
    next_dataset()
    return 1  # some key handlers like a truthy return

def _next_event_cb(*args):
    add_status_bar_text("NEXT event")
    next_event()
    return 1

def _prev_cb(*args):
    add_status_bar_text("PREV fired")
    prev_dataset()
    return 1

def _toggle_2fofc_cb(*args):
    add_status_bar_text("TOGGLE 2FOFC fired")
    toggle_2fofc()
    return 1

def _toggle_fofc_cb(*args):
    add_status_bar_text("TOGGLE FOFC fired")
    toggle_fofc()
    return 1

def _toggle_event_cb(*args):
    add_status_bar_text("TOGGLE event")
    toggle_event()
    return 1

# Keep hard references so the callback cannot be GC'd
NEXT_CB = _next_cb
PREV_CB = _prev_cb
NEXT_EVENT_CB = _next_event_cb
TOGGLE_2FOFC_CB = _toggle_2fofc_cb
TOGGLE_FOFC_CB = _toggle_fofc_cb
TOGGLE_EVENT_CB = _toggle_event_cb

add_key_binding_gtk4_py("v", 0, NEXT_CB, "Dataset: Next")
add_key_binding_gtk4_py("c", 0, PREV_CB, "Dataset: Prev")
add_key_binding_gtk4_py("t", 0, NEXT_EVENT_CB, "Dataset_event: Next")
add_key_binding_gtk4_py("g", 0, TOGGLE_2FOFC_CB, "toggle 2fofc map on/off")
add_key_binding_gtk4_py("h", 0, TOGGLE_EVENT_CB, "toggle event map on/off")
add_key_binding_gtk4_py("d", 0, TOGGLE_FOFC_CB, "toggle fofc map on/off")



_TOOLWIN = None  # keep a hard reference so it doesn't get GC'd
_DATASET_LABEL = None

def show_dataset_tool_window():
    global _TOOLWIN, _DATASET_LABEL
    if _TOOLWIN is not None:
        try:
            _TOOLWIN.present()
            return
        except Exception:
            _TOOLWIN = None  # if it was destroyed

    # If Coot is running as a Gtk.Application, parent our window to it
    app = None
    try:
        app = Gtk.Application.get_default()
    except Exception:
        pass

    win = Gtk.Window(title="Dataset Browser", application=app)
    win.set_default_size(320, 180)

    box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    box.set_margin_top(10); box.set_margin_bottom(10)
    box.set_margin_start(10); box.set_margin_end(10)

    # --- status area ---
    _DATASET_LABEL = Gtk.Label(xalign=0.0)
    box.append(_DATASET_LABEL)

    # --- row 1: dataset navigation ---
    row1 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    btn_prev = Gtk.Button(label="Prev dataset")
    btn_next = Gtk.Button(label="Next dataset")
    row1.append(btn_prev)
    row1.append(btn_next)

    # clicked handlers receive the button as first arg → ignore it
    btn_prev.connect("clicked", lambda *_: prev_dataset())
    btn_next.connect("clicked", lambda *_: next_dataset())

    # --- row 2: event maps ---
    row2 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    btn_next_evt = Gtk.Button(label="Next event map")
    btn_tog_evt  = Gtk.Button(label="Toggle event")
    row2.append(btn_next_evt)
    row2.append(btn_tog_evt)

    btn_next_evt.connect("clicked", lambda *_: next_event())
    btn_tog_evt.connect("clicked",  lambda *_: toggle_event())

    # --- row 3: map toggles ---
    row3 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    btn_2fo = Gtk.Button(label="Toggle 2Fo-Fc")
    btn_df  = Gtk.Button(label="Toggle Fo-Fc")
    row3.append(btn_2fo)
    row3.append(btn_df)

    btn_2fo.connect("clicked", lambda *_: toggle_2fofc())
    btn_df.connect("clicked",  lambda *_: toggle_fofc())

    box.append(row1)
    box.append(row2)
    box.append(row3)

    win.set_child(box)

    # Keep reference + show
    _TOOLWIN = win
    win.present()

def update_tool_window_labels():
    # If window not created yet, nothing to update
    global _DATASET_LABEL
    if _DATASET_LABEL is None:
        return

    dataset_name = _state.get("xtal_id", "N/A")
    dataset_num = _state.get("i") + 1 # shift index
    _DATASET_LABEL.set_text(f"Autorefine review: {dataset_name} ({dataset_num}/{N_DATASETS})")

# Call once to show it:
show_dataset_tool_window()

# Optional: load the first one immediately
_load_dataset(_state["i"])
