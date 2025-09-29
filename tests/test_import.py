import importlib


def test_import_package():
    """Basic smoke test: can import the package and check version."""
    pkg = importlib.import_module("sbtn_leaf")

    assert hasattr(pkg, "__version__")
    assert pkg.__version__.startswith("0.")


def test_import_functions():
    """Check that key functions are exposed at top-level."""
    import sbtn_leaf as sl

    # adjust these if you rename your public API
    assert hasattr(sl, "run_model")
    assert hasattr(sl, "ModelConfig")