import pytest
import os
import fragflows_db.database_init as db_init
from sqlalchemy import inspect
from fragflows_db.ingest import ingest_mx_processing_path
from fragflows_db.data_models import MXProcessingResult


@pytest.fixture(autouse=True)
def cleanup_test_db():
    """Automatically cleanup test database after test runs."""
    yield  # Run the test
    # Cleanup after test
    if os.path.exists("test_update.db"):
        os.unlink("test_update.db")

def test_database_update():
    db_init.init_db("test_update.db")
    with db_init.session_scope() as session:
        # Just test that we can create a session and it has the right tables
        test_path = ""
        ingest_mx_processing_path(session, test_path)
        # Check that the data was ingested
        result = session.query(MXProcessingResult).filter_by(xml_path=test_path).one_or_none()
        assert result is not None
        assert result.pipeline == "autoproc"

    # make sure that re-ingesting the same file doesn't create a duplicate entry
    with db_init.session_scope() as session:
        # Just test that we can create a session and it has the right tables
        test_path = ""
        ingest_mx_processing_path(session, test_path)
        result = session.query(MXProcessingResult).filter_by(xml_path=test_path).one_or_none()
        assert result is not None
        assert result.pipeline == "autoproc"

    # test fast_dp xml validation and ingestion
    with db_init.session_scope() as session:
        # Just test that we can create a session and it has the right tables
        test_path = ""
        ingest_mx_processing_path(session, test_path)
        result = session.query(MXProcessingResult).filter_by(xml_path=test_path).one_or_none()
        assert result is not None
        assert result.pipeline == "fast_dp"