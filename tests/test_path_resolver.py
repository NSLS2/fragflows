from deposition.utils import PathResolver
import os
import yaml


def test_path_resolver():
    # Create a temporary config file for testing
    test_config = {
        "path_mappings": [
            "/old/root/path",
            "/new/root/path"
        ]
    }
    
    with open("test_config.yaml", "w") as f:
        yaml.dump(test_config, f)
    
    resolver = PathResolver(mapping_config="test_config.yaml")
    
    # Test that the mapping is loaded correctly
    assert len(resolver.mappings) == 1
    assert resolver.mappings[0]['old_root'] == "/old/root/path"
    assert resolver.mappings[0]['new_root'] == "/new/root/path"
    
    # Test path resolution
    old_path = "/old/root/path/subdir/file.txt"
    expected_new_path = "/new/root/path/subdir/file.txt"
    resolved_path = resolver.resolve(old_path)
    assert resolved_path == expected_new_path
    
    # Clean up the temporary config file
    os.remove("test_config.yaml")

def test_path_resolver():
    resolver = PathResolver()
    assert resolver is not None
    old_path = "/mxprocessing/DLS_processed_2022_20241122/CypD-x0246/autoPROC.xml"
    new_path = "/nsls2/data/amx/proposals/2024-3/pass-316947/DLS_processed_2022_20241122/CypD-x0246/autoPROC.xml"
    resolved_path = resolver.resolve(old_path)
    assert resolved_path == new_path

    old_path = "/nsls2/data/amx/proposals/2024-3/pass-316947/DLS_processed_2022_20241122/CypD-x0246/autoPROC.xml"
    new_path = "/nsls2/data/amx/proposals/2024-3/pass-316947/DLS_processed_2022_20241122/CypD-x0246/autoPROC.xml"
    resolved_path = resolver.resolve(old_path)
    assert resolved_path == new_path