from sqlalchemy import (
    Column,
    Integer,
    String,
    ForeignKey,
    UniqueConstraint,
    DateTime,
)
from sqlalchemy.orm import declarative_base, relationship
from datetime import datetime, timezone

Base = declarative_base()


class Xtal(Base):
    __tablename__ = "xtal"

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True, index=True, nullable=False)

    hdf5_files = relationship(
        "HDF5File",
        back_populates="xtal",
        cascade="all, delete-orphan",
    )
    __table_args__ = (
        UniqueConstraint("name", name="uq_xtal_name"),
    )

    def __repr__(self):
        return f"<Xtal(name={self.name})>"


class HDF5File(Base):
    __tablename__ = "hdf5_file"

    id = Column(Integer, primary_key=True)

    xtal_id = Column(Integer, ForeignKey("xtal.id"), index=True, nullable=False)

    path = Column(String, nullable=False)
    filename = Column(String, nullable=False)

    xtal = relationship("Xtal", back_populates="hdf5_files")

    processing_results = relationship(
        "MXProcessingResult",
        back_populates="hdf5_file",
        cascade="all, delete-orphan",
    )

    __table_args__ = (
        UniqueConstraint("path", "filename", name="uq_hdf5_file"),
    )

    def __repr__(self):
        return f"<HDF5File({self.path}/{self.filename})>"


class MXProcessingResult(Base):
    __tablename__ = "mx_processing_result"

    id = Column(Integer, primary_key=True)

    hdf5_file_id = Column(
        Integer,
        ForeignKey("hdf5_file.id"),
        index=True,
        nullable=False,
    )

    xml_path = Column(String, unique=True, nullable=False)
    pipeline = Column(String)
    imported_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))
    checksum = Column(String, nullable=False)
    hdf5_file = relationship("HDF5File", back_populates="processing_results")

    def __repr__(self):
        return f"<MXProcessingResult(xml={self.xml_path})>"
    
  