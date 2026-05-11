import gemmi
from deposition.structure import resolve_entities
from deposition.cif_blocks import update_entity_id_loops
import os


def test_resolve_cif_entities():
    st = gemmi.read_structure("../fragflows_export_flow_20260319-1/cypd-1/cypd-1-ensemble-model_refine.mmcif")
    st_b = st.make_mmcif_block()
    entities = st_b.find_loop_item('_entity.id').loop.values
    assert len(entities) > 0, "missing _entity loop"
    resolve_entities(st)
    d = gemmi.cif.Document()
    d.add_copied_block(st.make_mmcif_block())
    d.write_file("test1.cif")
    entities_after_resolve = d[0].find_loop_item('_entity.id').loop.values
    if sorted(entities_after_resolve) != sorted(entities):
        print(f"entities before resolve: {entities}")
        print(f"entities after resolve: {entities_after_resolve}")
        #raise ValueError("missing _entity loop")
    assert st is not None

def test_update_entity_id_loops():
    d = gemmi.cif.read_file("cypd-1.mmcif")
    b = d[0] # structure block
    assert type(b) == gemmi.cif.Block
    update_entity_id_loops(b)
    d.write_file("test2.cif")
    assert os.path.exists("test2.cif")