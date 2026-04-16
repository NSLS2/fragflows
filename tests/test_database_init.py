import pytest
import os
import fragflows_db.database_init as db_init
from sqlalchemy import inspect

@pytest.fixture(autouse=True)
def cleanup_test_db():
    """Automatically cleanup test database after test runs."""
    yield  # Run the test
    # Cleanup after test
    if os.path.exists("test.db"):
        os.unlink("test.db")

def test_database_init():
    db_init.init_db("test.db")
    with db_init.session_scope() as session:
        # Just test that we can create a session and it has the right tables
        engine = session.get_bind()
        inspector = inspect(engine)
        tables_names = inspector.get_table_names()
        assert "xtal" in tables_names
        assert "hdf5_file" in tables_names
        assert "mx_processing_result" in tables_names
