
from .logger import Logger, LogLevel
from typing import Self
from copy import deepcopy
from packaging.version import Version
import sys, os, shutil
import snakemake
if Version(snakemake.__version__) >= Version("8.0.0"):
    import snakemake.cli as snakemake


DEFAULT_SMK_ARGS = [
    "--use-conda",
    "--conda-frontend", "mamba",
    "--rerun-incomplete",
]

class SmkHandler:
        def __init__(self, add_defaults = True):
            self.exe          = shutil.which("snakemake")
            self.logger       = Logger()
            self.targets      = []
            self.args         = []
            self.resources    = ["--resources"]
            self.add_defaults = add_defaults
            self.defaults     = DEFAULT_SMK_ARGS + [
                "--printshellcmds" if self.logger.level >= LogLevel.INFO  else None,
                "--verbose"        if self.logger.level >= LogLevel.DEBUG else None,
                "--quiet all"      if self.logger.level == LogLevel.QUIET else None,
            ]

        def add_arguments(self, args: [str]) -> Self:
            self.args.extend(args)
            return self

        def set_cores(self, cores: int | str) -> Self:
            if cores == "all":
                cores = int(os.cpu_count())
            if not str(cores).isnumeric() or int(cores) < 1:
                self.logger.error("Provided number of cores must either be a non-negative integer, or the string 'all'")
                sys.exit(1)

            self.logger.info(f"Setting {cores} cores.")
            self.args.extend(["--cores", f"{cores}"])
            self.resources.append(f"cores={cores}")
            return self

        def set_mem(self, mem_mb: int | None) -> Self:
            if mem_mb:
                self.logger.info(f"Setting Maximum Resident Set Size to {mem_mb} (MB).")
            else:
                self.logger.debug(f"Not setting any memory constraints, as --mem-mb is unset {mem_mb}")
            self.resources.append(f"mem_mb={mem_mb}")
            return self

        def wipe_defaults(self):
            self.defaults  = []
            self.resources = []

        def add_targets(self, targets) -> Self:
            if type(targets) is str:
                targets = [targets]
            self.targets.extend(targets)
            return self

        def add_resource(self, resources: [str]) -> Self:
            self.resources.extend(resources)
            return self

        def copy(self) -> Self:
            return deepcopy(self)

        def run(self, dry_run = False) -> None:
            if not self.add_defaults:
                self.wipe_defaults()

            smk_args = [self.exe] + self.targets + self.args + self.resources + self.defaults
            smk_args = [arg for arg in smk_args if arg]
            self.logger.info(f"Running snakemake with the following arguments: 'snakemake {' '.join(smk_args[1:])}'")
            try:
                if dry_run:
                    print('(dry-run)')
                    return
                sys.argv = smk_args
                snakemake.main()
            except SystemExit as status:
                self.logger.debug(f"Snakemake exit code {status.code}")
                if (status.code == 0):
                    self.logger.info("Done")
                else:
                    self.logger.error(f"Error encountered while runnine snakemake (exit code {status.code}) {status}")
                    sys.exit(status.code)
            except Exception as e:
                self.logger.error(f"Error encountered while runnine snakemake {e}")
                raise e