"""GUI integration with the shipped reader, independent of pygacode's age.

pygacode still supplies the ordinary field/flux reader. Do not monkey-patch it
or fall back to an arbitrary installed FULLT codec when validation fails.
"""
from cgyro_fullt_reader import (
    fullt_reader_for_case, reader_for_case, validate_loaded_case,
)


def open_fullt(case, suffix=".cgyro.fullt_asym", data=None):
    """Use one checked, slice-capable reader for raw MPI and FTZ layouts."""
    reader = fullt_reader_for_case(case, suffix)
    if reader is not None and data is not None:
        try:
            validate_loaded_case(reader, data)
        except Exception:
            reader.close()
            raise
    return reader


def open_compressed(case):
    """Select FTZ explicitly, including Fortran-style generated input values."""
    return reader_for_case(case)
