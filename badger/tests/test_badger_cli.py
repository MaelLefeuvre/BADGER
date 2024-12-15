import shlex, sys, os, tarfile
from pathlib import Path
import pytest
from badger import badger

TEST_CONFIG_ARGS="--config gargammel=\"{'coverage': '0.005', 'comp-endo': '0.99', 'comp-cont': '0.01', 'compt_bact': '0.00'}\""
SMK_SYMLINK_DIRS=["config", "resources", "workflow"]
MIN_REQUIRED_CORES="18"
MIN_REQUIRED_MEM="32000"

def _run_badger_main(command: str = ""):
    print(f"Running badger.main with the following arguments: {command}")
    if " -- " not in command:
        command = f"{command} --"
    command = f"{command} {TEST_CONFIG_ARGS}"
    return badger.main([sys.argv[0]] + shlex.split(command))


@pytest.fixture(scope = 'module')
def setup_env(tmp_path_factory):
    cwd    = os.getcwd()
    tmpdir =  tmp_path_factory.mktemp("test-badger")
    for smk_dir in SMK_SYMLINK_DIRS:
        symlink_path = os.path.join(tmpdir, smk_dir)
        symlink_target = os.path.relpath(smk_dir, tmpdir)
        Path(symlink_path).symlink_to(symlink_target)
    os.chdir(tmpdir)
    yield tmpdir
    os.chdir(cwd)

def test_dry_run_setup(setup_env):
    _run_badger_main("setup -- --forcerun fetch_data --dry-run")

def test_dry_run_all(setup_env):
    _run_badger_main(f"run --cores {MIN_REQUIRED_CORES} --mem-mb {MIN_REQUIRED_MEM} -- --forcerun fetch_data --dry-run")

def test_dry_run_archive(setup_env):
    _run_badger_main("archive -- --forcerun fetch_data --dry-run")

def test_dry_run_loop_pipeline(setup_env):
    _run_badger_main(f"loop-pipeline --cores {MIN_REQUIRED_CORES} --mem-mb {MIN_REQUIRED_MEM} -i 3 -- --forcerun fetch_data --dry-run")

@pytest.mark.dependency(name="create-conda-envs")
@pytest.mark.extensive
@pytest.mark.network
def test_create_conda_envs(setup_env):
    command = "setup --no-fetch-data -- --forcerun fetch_data"
    print(command)
    badger.main([sys.argv[0]] + shlex.split(command))

@pytest.mark.dependency(name="fetch_data", depends=["create-conda-envs"])
@pytest.mark.extensive
@pytest.mark.network
def test_fetch_data(setup_env):
    command = "setup --cores 1 --no-create-envs"
    print(command)
    badger.main([sys.argv[0]] + shlex.split(command))

@pytest.mark.dependency(name="run", depends = ["fetch_data", "create-conda-envs"])
@pytest.mark.extensive
def test_run_module(setup_env):
    command = "run"
    print(command)
    badger.main([sys.argv[0]] + shlex.split(command))

@pytest.mark.dependency(name="archive", depends = ["run"])
@pytest.mark.extensive
def test_archive_module(setup_env):
    command = "archive"
    print(command)
    badger.main([sys.argv[0]] + shlex.split(command))
