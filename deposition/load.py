import xmltodict
from .utils import flatten_dict
from typing import Any, Union
from .mapping import ALL_SHELL_TO_MODEL, HIGH_SHELL_TO_MODEL
from .models import ReflectionStats
import gemmi


def ispyb_xml_to_dict(xmlpath: str) -> dict:
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
    print(overall_shell)
    high_shell = resolution_shells_to_payload(shells, desired_shell="outerShell")
    reflection_stats_dict = overall_shell | high_shell
    print(reflection_stats_dict)
    reflection_stats = ReflectionStats(**reflection_stats_dict)
    out_block = gemmi.cif.Block(blockname)
    print(reflection_stats.dict())
    for k, v in reflection_stats.dict(by_alias=True).items():
        out_block.set_pair(k, v)
    return out_block
