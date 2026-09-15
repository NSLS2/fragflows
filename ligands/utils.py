from pathlib import Path
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions
)
from deposition.cif_blocks import deduplicate_cif_loops
import gemmi

def enumerate_stereo_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    opts = StereoEnumerationOptions(
        onlyUnassigned=True,   # preserve stereo that's already specified
        unique=True,           # remove duplicates
        maxIsomers=1024,       # safety limit
        tryEmbedding=False
    )
    if mol is None:
        raise ValueError(f"invalid smiles {smiles}")
    isomers = [Chem.MolToSmiles(isomer, isomericSmiles=True) for isomer in EnumerateStereoisomers(mol, options=opts)]
    return isomers

def consolidate_lig_restraints(filelist: list[Path])->gemmi.cif.Document:

    print(filelist)
    doclist = [gemmi.cif.read_file(str(f)) for f in filelist if f.suffix == '.cif']
    for doc in doclist:
        if len(doc) != 2:
            raise ValueError(f"Expected 2 blocks in CIF file, found {len(doc)}")
        
    block_dict = {}
    for f in doclist:
        for b in f:
            if block_dict.get(b.name) is None:
                block_dict[b.name] = b
            elif block_dict.get(b.name) is not None and b.name != 'comp_list':
                raise ValueError(f"Duplicate block found for {b.name} which is not 'comp_list'")
            else:
                # merged rows land on the target (b); the returned source copy has the loop stripped and is discarded
                _, block_dict[b.name] = deduplicate_cif_loops(block_dict.get(b.name), b, '_chem_comp.id')

    doc = gemmi.cif.Document()
    for block in block_dict.values():
        doc.add_copied_block(block)

    return doc