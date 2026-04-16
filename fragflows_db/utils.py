import os
from typing import List
from pathlib import Path
import hashlib

def sha256sum(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()

def get_xml_paths(
        datapath: str,
        xml_files: List=['autoPROC.xml', 'fast_dp.xml'],
        dirs_to_avoid: List=['dozor']
    ) -> List[str]:
    xml_paths = []
    for r,d,f in os.walk(datapath):
        d = [d_ for d_ in d if d_ not in dirs_to_avoid]
        for f_ in f:
            if f_ in xml_files:
                print(f'{r}/{f_}')
                xml_paths.append(f'{r}/{f_}')
    return xml_paths

def get_hdf5_paths(
        datapath: str,
        hdf5_tag: List=['master.h5'],
        include_rasters: bool = False,
        dirs_to_avoid: List=['dozor', 'autoProcOutput']
)-> List[str]:
    hdf5_paths = []
    for r,d,f in os.walk(datapath):
        d = [d_ for d_ in d if d_ not in dirs_to_avoid]
        for f_ in f:
            if include_rasters:
                if f_.endswith(hdf5_tag):
                    print(f'{r}/{f_}')
                    hdf5_paths.append(f'{r}/{f_}')
            else:
                if f_.endswith(hdf5_tag) and 'Raster' not in f_:
                    print(f'{r}/{f_}')
                    hdf5_paths.append(f'{r}/{f_}')
    return hdf5_paths
