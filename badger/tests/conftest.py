import subprocess
import pytest
import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--extensive", action="store_true", default=False, help="Test BADGER extensively, by running every module on a dataset. This can take a very long time..."
    )
    parser.addoption(
        "--network", action="store_true", default=False, help="Test BADGER rules that require a network connection"
    )
    parser.addoption(
        "--urls", action="store_true", default=False, help="Test the validity and accessibility of the URLs used by BADGER"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "extensive: mark test as computationally intensive")
    config.addinivalue_line("markers", "network: mark test as requiring a network connection")
    config.addinivalue_line("markers", "urls: mark test which checks the validity and accessibility of every URL used by BADGER")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--extensive"):
        return
    #if config.getoption("--network"):
    #    return
    #if config.getoption("--urls"):
    #    return
    skip_extensive = pytest.mark.skip(reason="need --extensive option to run")
    skip_network   = pytest.mark.skip(reason="need --network option to run")
    skip_urls      = pytest.mark.skip(reason="need --urls option to run")
    for item in items:
        if "urls" in item.keywords and not config.getoption("--urls"):
            item.add_marker(skip_urls)
            continue
        if "network" in item.keywords and not config.getoption("--network"):
            item.add_marker(skip_network)
            continue
        if "extensive" in item.keywords and not config.getoption("--extensive"):
            item.add_marker(skip_extensive)
            continue

