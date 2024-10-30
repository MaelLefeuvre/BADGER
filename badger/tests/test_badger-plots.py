import shlex, sys, tarfile, subprocess, pytest
from os import path, chdir, mkdir

TEST_TARBALL="test-data/test-archives.tar.xz"

def _shell_cmd(command: str = "", stdout = subprocess.PIPE, stderr = subprocess.PIPE):
    print(f"\nRunning command {command}")
    return subprocess.run(shlex.split(command), stdout=stdout, stderr=stderr)

@pytest.fixture(scope = 'session')
def setup_plot_testdata(tmp_path_factory):
    tmpdir       = tmp_path_factory.mktemp("test-badger-plots")
    test_tarball = path.join(path.dirname(__file__), TEST_TARBALL)
    tar = tarfile.open(test_tarball, mode='r')
    tar.extractall(tmpdir, filter = "data")
    chdir(tmpdir)
    mkdir("plots")
    yield tmpdir

@pytest.mark.dependency(name="make-input", depends=[])
def test_badger_plots_make_input(setup_plot_testdata):
    command = 'badger-plots make-input -d test-archives -s 0.02X,0.04X,0.06X'
    result = _shell_cmd(command, open("plots/input.yml", 'w+'))
    assert(result.returncode==0)

@pytest.mark.dependency(name="template", depends=["make-input"])
def test_badger_plots_template(setup_plot_testdata):
    command = 'badger-plots template -i plots/input.yml -P test-archives/pedigree-codes.txt --output-dir plots'
    result = _shell_cmd(command, open("plots/params.yml", 'w+'))
    assert(result.returncode==0)

@pytest.mark.dependency(name="plot", depends=["make-input", "template"])
def test_badger_plots_plot(setup_plot_testdata):
    command = 'badger-plots plot --yaml plots/params.yml --threads 1'
    result = _shell_cmd(command)
    assert(result.returncode==0)
