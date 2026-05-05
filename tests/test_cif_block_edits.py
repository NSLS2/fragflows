from deposition.cif_blocks import convert_cif_pairs_to_loop
import gemmi

def test_convert_cif_pairs_to_loop():
    block = gemmi.cif.Block('test')
    block.set_pair("_refine.pdbx_ls_cross_validation_method", "THROUGHOUT")
    block.set_pair("_refine.pdbx_method_to_determine_struct", "\"MOLECULAR REPLACEMENT\"")
    new_block = convert_cif_pairs_to_loop(block, "_refine")
    assert new_block is not None
    assert new_block.find_loop_item("_refine.pdbx_ls_cross_validation_method") is not None
    assert new_block.find_loop_item("_refine.pdbx_method_to_determine_struct") is not None
    assert new_block.find_loop("_refine.pdbx_ls_cross_validation_method")[0] == "THROUGHOUT"