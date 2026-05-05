import gemmi
from deposition.load import ispyb_xml_to_dict

PIPELINE_PARAMETERS = {
    "autoPROC 1.0.5 (20240123 - see: http://www.globalphasing.com/autoproc/)": [
        {
            "name": "autoPROC",
            "classification": "\"data reduction\"",
            "version": "1.0.5",
            "location": "http://www.globalphasing.com/autoproc/",
        },
        {"name": "Aimless", "classification": "\"data scaling\""},
    ],
    "fast_dp": [
        {"name": "FAST_DP", "classification": "\"data reduction\""},
        {"name": "Aimless", "classification": "\"data scaling\""},
        {"name": "XDS", "classification": "\"data reduction\""},
    ],
    "default": [
        {"name": "DIMPLE", "classification": "\"phasing\"", "type": "program"},
        {"name": "PHASER", "classification": "\"phasing\""},
        {"name": "Coot", "classification": "\"model building\""},
        {
            "name": "gemmi",
            "classification": "\"data extraction\"",
            "version": f"{gemmi.__version__}",
            "location": "https://github.com/project-gemmi/gemmi",
        },
    ],
}


def xml_path_to_pipeline_programs(xml_path: str) -> str:
    ispyb_dict = ispyb_xml_to_dict(xml_path)
    program = ispyb_dict[
        "AutoProcContainer.AutoProcProgramContainer.AutoProcProgram.processingPrograms"
    ]
    # program is a string which we need to clean and extract to use as a key
    return [*PIPELINE_PARAMETERS.get(program, []), *PIPELINE_PARAMETERS['default']]
