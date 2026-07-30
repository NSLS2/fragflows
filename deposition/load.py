import xmltodict
from .utils import flatten_dict
from typing import Any, Union
from .mapping import ALL_SHELL_TO_MODEL, HIGH_SHELL_TO_MODEL
from .models import ReflectionStats
import gemmi
from deposition.utils import PathResolver

_path_resolver = PathResolver()

def ispyb_xml_to_dict(xmlpath: str) -> dict:
    xmlpath = _path_resolver.resolve(xmlpath)
    with open(xmlpath, "rb") as f:
        data = xmltodict.parse(f)
    return flatten_dict(data)


def resolution_shells_to_payload(
    shells: list[dict[str, Any]], desired_shell: str = None
) -> dict:
    for shell in shells:
        if shell["scalingStatisticsType"] == desired_shell:
            chosen_shell = shell
            break
    if chosen_shell is None:
        raise ValueError("resolution shell not found")
    result: dict[str, Union[str, float]] = {}
    if desired_shell == "overall":
        key_map = ALL_SHELL_TO_MODEL
    elif desired_shell == "outerShell":
        key_map = HIGH_SHELL_TO_MODEL
    else:
        raise ValueError("Unable to find matching shell")
    for k, v in chosen_shell.items():
        if k in key_map:
            # add leading 0 to string if missing, artifact of autoPROC.xml
            if v.split('.')[0] == '':
                v = ''.join(['0',v])
            result[key_map[k]] = v
    return result


def ispyb_xml_to_cif_block(
    xmlpath: str, blockname="reflection_stats"
) -> gemmi.cif.Block:
    flat_ispyb_xml_dict = ispyb_xml_to_dict(xmlpath)
    shells = flat_ispyb_xml_dict[
        "AutoProcContainer.AutoProcScalingContainer.AutoProcScalingStatistics"
    ]
    overall_shell = resolution_shells_to_payload(shells, desired_shell="overall")
    
    high_shell = resolution_shells_to_payload(shells, desired_shell="outerShell")
    reflection_stats_dict = overall_shell | high_shell

    if hasattr(ReflectionStats, "model_validate"):
        reflection_stats = ReflectionStats.model_validate(reflection_stats_dict)
    else:
        reflection_stats = ReflectionStats(**reflection_stats_dict)
    out_block = gemmi.cif.Block(blockname)
    if hasattr(reflection_stats, "model_dump"):
        payload = reflection_stats.model_dump(by_alias=True)
    else:
        payload = reflection_stats.dict(by_alias=True)

    for k, v in payload.items():
        out_block.set_pair(k, v)
    return out_block

def cif_to_payload(cif_block: gemmi.cif.Block, ignored_categories: list[str]=['_atom_site','_entity_poly_seq','_refln']) -> dict:
    payload = {}
    for item in cif_block:
        if item.pair:
            # Simple key-value pair
            key, value = item.pair
            if key in payload:
                raise ValueError(f"Duplicate key found in CIF block: {key}")
            payload[key] = value
            
        elif item.loop:
            # Loop with multiple columns (tags) and rows
            loop = item.loop
            tags = loop.tags  # List of column headers
            if any(ignored_category == tags[0].split('.')[0] for ignored_category in ignored_categories):
                # Skip ignored categories
                continue
            values = loop.values  # Flat list of all values
            
            # Chunk values by number of tags
            num_cols = len(tags)
            rows = [values[i:i + num_cols] for i in range(0, len(values), num_cols)]
            
            # Store as list of dicts, one dict per row
            loop_data = [dict(zip(tags, row)) for row in rows]
            
            # You can store as single dict if one row, or list if multiple
            if len(loop_data) == 1:
                # Single row: merge tags as separate keys or nest under loop name
                payload.update(loop_data[0])
            else:
                # Multiple rows: use loop name (first tag prefix) as key
                loop_name = tags[0].rsplit('.', 1)[0]  # e.g., "_atom_site.label" → "_atom_site"
                payload[loop_name] = loop_data

    resolution_shells = payload.get('_refine_ls_shell', [])

    if isinstance(resolution_shells, list) and len(resolution_shells) > 0:
        high_resolution_shell = min(resolution_shells, key=lambda x: float(x.get('_refine_ls_shell.d_res_low')))
        for key, value in high_resolution_shell.items():
            if key in payload:
                raise ValueError(f"Duplicate key found in CIF block: {key}")
            payload[key] = value

    refln_stats = payload.get('_reflns', [])
    if isinstance(refln_stats, list) and len(refln_stats) == 1:
        overall_refln_stats = refln_stats[0]
        for key, value in overall_refln_stats.items():
            if key in payload:
                raise ValueError(f"Duplicate key found in CIF block: {key}")
            payload[key] = value

    refln_shell_stats = payload.get('_reflns_shell', [])
    if isinstance(refln_shell_stats, list) and len(refln_shell_stats) == 1:
        overall_refln_shell_stats = refln_shell_stats[0]
        for key, value in overall_refln_shell_stats.items():
            if key in payload:
                raise ValueError(f"Duplicate key found in CIF block: {key}")
            payload[key] = value

    return payload