"""Utilities for resolving repository-relative paths.

This module centralises filesystem locations used throughout the
:mod:`sbtn_leaf` package so that callers can reliably locate bundled data
files regardless of the current working directory.  All helpers return
:class:`pathlib.Path` objects for convenient interoperability with the
standard library and third-party packages.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Union, overload

PathLike = Union[str, Path]


def _resolve_root() -> Path:
    """Return the repository root directory.

    The project structure follows the :pep:`517` recommendation of placing
    package sources under ``src/``.  This helper walks two parents up from
    the module file (``.../src/sbtn_leaf/paths.py``) to land on the
    repository root directory.
    """

    return Path(__file__).resolve().parent.parent.parent


_PROJECT_ROOT = _resolve_root()
_DATA_DIR = _PROJECT_ROOT / "data"
_DOCUMENTATION_DIR = _PROJECT_ROOT / "documentation"
_EXAMPLES_DIR = _PROJECT_ROOT / "examples"
_LEAFS_DIR = _PROJECT_ROOT / "LEAFs"


def project_root() -> Path:
    """Return the absolute path to the repository root."""

    return _PROJECT_ROOT


def data_dir() -> Path:
    """Return the absolute path to the repository ``data`` directory."""

    return _DATA_DIR


def documentation_dir() -> Path:
    """Return the absolute path to the ``documentation`` directory."""

    return _DOCUMENTATION_DIR


def examples_dir() -> Path:
    """Return the absolute path to the ``examples`` directory."""

    return _EXAMPLES_DIR

def LEAFs_dir() -> Path:
    """Return the absolute path to the ``LEAFs`` directory."""

    return _LEAFS_DIR


@overload
def data_path(*parts: PathLike) -> Path:
    ...


@overload
def data_path(parts: Iterable[PathLike]) -> Path:
    ...


def data_path(*parts: PathLike | Iterable[PathLike]) -> Path:
    """Return a path inside the repository ``data`` directory.

    Parameters
    ----------
    parts:
        Optional path segments that will be joined beneath ``data``.  The
        helper accepts either variadic positional arguments or a single
        iterable for convenience.  ``Path`` objects are returned to maximise
        compatibility with consumers that accept ``os.PathLike`` values.
    """

    if len(parts) == 1 and isinstance(parts[0], Iterable) and not isinstance(parts[0], (str, bytes, Path)):
        parts_iterable = parts[0]
    else:
        parts_iterable = parts

    return _DATA_DIR.joinpath(*map(Path, parts_iterable))


@overload
def leaf_path(*parts: PathLike) -> Path:
    ...


@overload
def leaf_path(parts: Iterable[PathLike]) -> Path:
    ...


def leaf_path(*parts: PathLike | Iterable[PathLike]) -> Path:
    """Return a path inside the repository ``LEAFs`` directory.

    Parameters
    ----------
    parts:
        Optional path segments that will be joined beneath ``LEAFs``.  The
        helper accepts either variadic positional arguments or a single
        iterable for convenience.  ``Path`` objects are returned to maximise
        compatibility with consumers that accept ``os.PathLike`` values.
    """

    if len(parts) == 1 and isinstance(parts[0], Iterable) and not isinstance(parts[0], (str, bytes, Path)):
        parts_iterable = parts[0]
    else:
        parts_iterable = parts

    return _LEAFS_DIR.joinpath(*map(Path, parts_iterable))


def project_path(*parts: PathLike) -> Path:
    """Return a path anchored at the repository root."""

    return _PROJECT_ROOT.joinpath(*map(Path, parts))
