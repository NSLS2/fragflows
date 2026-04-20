from pathlib import Path
from lxml import etree
from deposition import load
from fragflows_db.crud import get_or_create_xtal, get_or_create_hdf5_file, get_or_create_mx_processing_result
from fragflows_db.utils import sha256sum



def ingest_mx_processing_path(session, filepath: str):

    if Path(filepath).is_symlink():
        return  # silently ignore autoPROC symlinks

    filepath = Path(filepath).resolve()
    
    if not filepath.exists():
        raise Exception(f"XML file {filepath} does not exist")
    
    filepath = str(filepath)
    
    filename = Path(filepath).name.lower().split('.')[0]

    # Get the directory containing the XSD schema files
    schema_dir = Path(__file__).parent / "xml_schemas"

    if filename == 'autoproc':
        xsd_path = schema_dir / "autoproc.xsd"
        xsd_doc = etree.parse(str(xsd_path))
    elif filename == 'fast_dp':
        xsd_path = schema_dir / "fast_dp.xsd"
        xsd_doc = etree.parse(str(xsd_path))
    else:
        raise Exception(f'unknown processing pipeline xml: {filepath}')

    schema = etree.XMLSchema(xsd_doc)

    xml_doc = etree.parse(filepath)

    if not schema.validate(xml_doc):
        print(f"Invalid schema for {filepath}")
        for error in schema.error_log:
            print(error.message)
        return

    mx_processing_dict = load.ispyb_xml_to_dict(filepath)
    hdf5_directory = mx_processing_dict['AutoProcContainer.AutoProcScalingContainer.AutoProcIntegrationContainer.Image.fileLocation']
    hdf5_filename = mx_processing_dict['AutoProcContainer.AutoProcScalingContainer.AutoProcIntegrationContainer.Image.fileName']
    hdf5_path = f'{hdf5_directory}/{hdf5_filename}'
    hdf5_path = Path(hdf5_path).resolve()
    if not hdf5_path.exists():
        raise Exception(f"HDF5 file {hdf5_path} does not exist for XML file {filepath}")
    hdf5_path = str(hdf5_path)
    xtal_name = Path(mx_processing_dict['AutoProcContainer.AutoProcScalingContainer.AutoProcIntegrationContainer.Image.fileLocation']).parts[-3]

    xtal = get_or_create_xtal(session, xtal_name)
    hdf5_file = get_or_create_hdf5_file(session, xtal, hdf5_path, hdf5_filename)
    mx_processing_result = get_or_create_mx_processing_result(
                                session,
                                hdf5 = hdf5_file,
                                xml_path = filepath,
                                pipeline = filename,
                                checksum = sha256sum(filepath)
    )
