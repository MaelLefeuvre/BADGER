import subprocess
import pytest
import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--extensive", action="store_true", default=False, help="test BADGER extensively, by running every module on a dataset."
    )
    parser.addoption(
        "--network", action="store_true", default=False, help="Test BADGER rules that require a network connection"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "extensive: mark test as computationally intensive")
    config.addinivalue_line("markers", "network: mark test as requiring a network connection")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--network"):
        return
    if config.getoption("--extensive"):
        return
    skip_extensive = pytest.mark.skip(reason="need --extensive option to run")
    skip_network   = pytest.mark.skip(reason="need --network option to run")
    for item in items:
        if "extensive" in item.keywords:
            item.add_marker(skip_extensive)
    
        if "network" in item.keywords:
            item.add_marker(skip_network)
