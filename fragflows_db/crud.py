from fragflows_db.data_models import Xtal, HDF5File, MXProcessingResult
import pandas as pd
from deposition import load
import h5py

def get_or_create_xtal(session, name: str) -> Xtal:
    xtal = session.query(Xtal).filter_by(name=name).one_or_none()
    if xtal:
        return xtal

    xtal = Xtal(name=name)
    session.add(xtal)
    return xtal


def get_or_create_hdf5_file(session, xtal: Xtal, path: str, filename: str) -> HDF5File:
    hdf5 = session.query(HDF5File).filter_by(
        path=path,
        filename=filename
    ).one_or_none()

    if hdf5:
        return hdf5

    hdf5 = HDF5File(
        xtal=xtal,
        path=path,
        filename=filename,
    )
    session.add(hdf5)
    return hdf5


def get_or_create_mx_processing_result(
    session,
    hdf5: HDF5File,
    xml_path: str,
    pipeline: str,
    checksum: str,
    **kwargs,
) -> MXProcessingResult:

    result = session.query(MXProcessingResult).filter_by(
        xml_path=xml_path
    ).one_or_none()

    if result:
        return result

    result = MXProcessingResult(
        hdf5_file=hdf5,
        xml_path=xml_path,
        pipeline=pipeline,
        checksum=checksum,
        **kwargs,
    )

    session.add(result)
    return result

import pandas as pd

def get_mx_processing_results_df(session):
    query = (
        session.query(
            MXProcessingResult.id.label("mx_processing_result_id"),
            MXProcessingResult.xml_path,
            MXProcessingResult.pipeline,
            MXProcessingResult.imported_at,
            MXProcessingResult.checksum,

            HDF5File.id.label("hdf5_file_id"),
            HDF5File.path.label("hdf5_path"),
            HDF5File.filename.label("hdf5_filename"),

            Xtal.id.label("xtal_id"),
            Xtal.name.label("xtal_name"),
        )
        .join(MXProcessingResult.hdf5_file)
        .join(HDF5File.xtal)
    )

    return pd.read_sql(query.statement, session.bind)


def get_stats_from_xml(xml_filepath: str):
    d = load.ispyb_xml_to_dict(xml_filepath)
    symm = d['AutoProcContainer.AutoProc.spaceGroup']
    a = d['AutoProcContainer.AutoProc.refinedCell_a']
    b = d['AutoProcContainer.AutoProc.refinedCell_b']
    c = d['AutoProcContainer.AutoProc.refinedCell_c']
    alpha = d['AutoProcContainer.AutoProc.refinedCell_alpha']
    beta = d['AutoProcContainer.AutoProc.refinedCell_beta']
    gamma = d['AutoProcContainer.AutoProc.refinedCell_gamma']
    overall = [
        d_ for d_ in d['AutoProcContainer.AutoProcScalingContainer.AutoProcScalingStatistics'] if d_['scalingStatisticsType'] == 'overall'
    ][0]
    overall_hi = overall['resolutionLimitHigh']
    overall_lo = overall['resolutionLimitLow']
    overall_r_mrg = overall['rMerge']
    overall_comp = overall['completeness']
    overall_mult = overall['multiplicity']
    outer = [
        d_ for d_ in d['AutoProcContainer.AutoProcScalingContainer.AutoProcScalingStatistics'] if d_['scalingStatisticsType'] == 'outerShell'
    ][0]
    outer_hi = outer['resolutionLimitHigh']
    outer_lo = outer['resolutionLimitLow']
    outer_r_mrg = outer['rMerge']
    outer_comp = outer['completeness']
    outer_mult = outer['multiplicity']
    outer_cchalf = outer['ccHalf']
    return {
        'hi': overall_hi,
        'lo': overall_lo,
        'r_mrg': overall_r_mrg,
        'overall_comp': overall_comp,
        'overall_mult': overall_mult,
        'outer_hi': outer_hi,
        'outer_lo': outer_lo,
        'outer_comp': outer_comp,
        'outer_mult': outer_mult,
        'outer_cchalf': outer_cchalf,
        'symm': symm,
        'a': a,
        'b': b,
        'c': c,
        'alpha': alpha,
        'beta': beta,
        'gamma': gamma,
    }

def extend_dataframe_mx_stats(df: pd.DataFrame):
    for idx, row in df.iterrows():
        d = get_stats_from_xml(row['xml_path'])
        for key, value in d.items():
            df.loc[idx, key] = value
    return df

def reflection_file_to_df(df: pd.DataFrame):
    for idx, row in df.iterrows():
        d = load.ispyb_xml_to_dict(row['xml_path'])
        r = [m for m in d['AutoProcContainer.AutoProcProgramContainer.AutoProcProgramAttachment'] if m['fileType'] == 'Result' and m['fileName'].endswith('.mtz')][0]
        df.loc[idx, 'mtz_path'] = f"{r['filePath']}/{r['fileName']}"
    return df


def extract_collection_info_from_hdf5(hdf5_filepath):
     f = h5py.File(hdf5_filepath)
     return {
         'data_collection_date': f['/entry/instrument/detector/detectorSpecific/data_collection_date'][()].decode(),
         'det_description': f['/entry/instrument/detector/description'][()].decode(),
         'det_serial_no': f['/entry/instrument/detector/detector_number'][()].decode(),
         'wavelength': f['/entry/instrument/beam/incident_wavelength'][()],
     }

def extend_dataframe_collection_info(df: pd.DataFrame):
    for idx, row in df.iterrows():
        print(f"fetching collection info from {row['hdf5_path']}...")
        d = extract_collection_info_from_hdf5(row['hdf5_path'])
        for key, value in d.items():
            df.loc[idx, key] = value
    return df