# Contributing

## Submitting an issue

If you experience any problem when using this software, feel free to directly submit an issue on this github repository. When doing so, we would be grateful if you could provide us with as much information as possible, so that we can identify and resolve the problem quickly. Mainly:
1. A clear and concise description of the bug you are experiencing
2. Specify the operating system you are working, as well as any additional information helpful to identify your environment (`uname -a`)
3. A list of steps required to reproduce the behaviour
4. A description of what you *actually* expected to happen when running these steps
5. Any relevant output log message which may help decipher the problem
6. When applicable, screenshots and/or input test files required to execute a minimally reproducible example are welcome.
7. Any additional context you deem useful to fully describe the issue.

Note that issue templates are provided when submitting an issue on github to help describe the problem.

## Development

### Running the extensive test-suite

Although the development environment of BADGER is supported by a set of integration tests, some of them will require a network access and can be quite expensive in terms of computing requirements. Thus, these additionnal tests are not executed by default when running the test suite, and must be explicitly requested. To execute the *'extensive'* test suite, developers may run the usual test command with the following additional arguments:

```bash
./badger/install.sh test --network --extensive
```
This command will run the following additional tests:

- `test_badger_cli.py::create-conda-envs` : Will attempt to create all of the conda environments required to run BADGER, using `badger setup`
- `test_badger_cli.py::fetch_data` : Will attempt to download all datasets required to run BADGER, using `badger setup`
- `test_badger_cli.py::run_module` : Will attempt to execute a single run of BADGER, using `badger run` module, and the provided test data
- `test_badger_cli.py::archive`: Will attempt to archive a single run of BADGER, using `badger archive`, and the provided test data (located in: `badger/tests/test-data/`)

Note that these integration tests may depend one-another and should be executed in a specific order. This feature is automatically handled by the program [pytest-dependency](https://pypi.org/project/pytest-dependency/)

### Pull-request process
1. Increase the version numbers in the `badger/envs/` conda environment definition files, `badger.install.sh` script. Please use the [SemVer](https://semver.org/) versioning scheme.
2. Check that the specified minimum required snakemake version is still relevant. If not, update the specified version at the beginning of `workflow/Snakefile`
2. Update the README.md and documentation with details reflecting any changes in the interface. Ensure any specified version numbers are updated accordingly.
4. Ensure the program can be installed properly and comply with basic integration tests using the installation script (`./badger/install.sh`)
5. We strongly recommend that you run the [Extensive test suite](#running-extensive-test-suite) before requesting a pull request. Note that pull requests will only be merged once all unit-tests integration-tests, and github action jobs have been successfully passed.

