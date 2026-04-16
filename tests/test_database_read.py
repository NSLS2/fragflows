import pytest
import os
import fragflows_db.database_init as db_init
from sqlalchemy import inspect
from fragflows_db.ingest import ingest_mx_processing_path
from fragflows_db.data_models import MXProcessingResult
from fragflows_db.utils import get_xml_paths
from fragflows_db.crud import get_mx_processing_results_df
from fragflows_db.crud import extend_dataframe_mx_stats

@pytest.fixture(autouse=True)
def cleanup_test_db():
    """Automatically cleanup test database after test runs."""
    yield  # Run the test
    # Cleanup after test
    if os.path.exists("test_read.db"):
        os.unlink("test_read.db")

def test_database_read():
    db_init.init_db("test_read.db")
    DATA_DIRECTORY = ""

    xml_paths = get_xml_paths(DATA_DIRECTORY)
    assert len(xml_paths) > 0

    for xp in xml_paths:
        print(f"found {len(xp)} xml files in {DATA_DIRECTORY}")
        with db_init.session_scope() as session:
            try:
                ingest_mx_processing_path(session, xp)
                print(f"successfully ingested {xp}")
            except Exception as e:
                print(f"failed to ingest {xp}: {e}")

        #test reading as df
    with db_init.session_scope() as session:
        df = get_mx_processing_results_df(session)
        print(df)
        df = extend_dataframe_mx_stats(df)
        print(df)
        assert len(df) > 0
        assert "pipeline" in df.columns
        assert "xml_path" in df.columns
