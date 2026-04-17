from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from contextlib import contextmanager
import os

_engine = None
SessionLocal = sessionmaker(future=True)

def init_db(db_path="fragflows.db"):
    global _engine

    _engine = create_engine(
        f"sqlite:///{db_path}",
        echo=False,
        future=True
    )

    from fragflows_db.data_models import Base

    Base.metadata.create_all(_engine)
    SessionLocal.configure(bind=_engine)

@contextmanager
def session_scope():
    """Provide a transactional scope around a series of operations."""
    session = SessionLocal()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()