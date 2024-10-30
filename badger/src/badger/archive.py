from .logger import Logger
from .constants import DEFAULT_CORES, DEFAULT_MEM_MB
import snakemake


from .smk_handler import SmkHandler
from .module_components import BadgerCli

class ArchiveCli(BadgerCli):
    def help(self):
        return "Archive the results of BADGER"

def setup_cli(subparsers, **kwargs) -> None:
    parser = ArchiveCli(subparsers, **kwargs)
    parser.add_argument("-c", "--cores", action='store',
      help = "Set the maximum number of CPU cores for snakemake (Default: %(default)s)", 
      default = DEFAULT_CORES
    )
    parser.add_argument("-m", "--mem-mb", action='store',
      help = "Set the maximum amount of Preemptible memory (in MB) for snakemake (Default: %(default)s)",
      default = DEFAULT_MEM_MB
    )
    parser.add_argument("--no-defaults", action='store_true',
      help = f"Don't provide snakemake with default arguments."
    )


def run(args, smk_args, **kwargs):
    smk = SmkHandler(add_defaults = not args.no_defaults)
    smk = smk.set_cores(args.cores).add_targets('archive').set_mem(args.mem_mb).add_arguments(smk_args)
    smk = smk.add_arguments(["--rerun-triggers", "mtime"])
    smk.run()