from .logger import Logger, LogLevel
from .constants import DEFAULT_CORES, DEFAULT_MEM_MB
import snakemake

from .smk_handler import SmkHandler
from .module_components import BadgerCli, BadgerSmkCli

class RunCli(BadgerCli, BadgerSmkCli):
    def help(self): return "Run a single run of BADGER"
    def description(self): return "Run a single run of BADGER in the current directory"

def setup_cli(subparsers, **kwargs):
    parser = RunCli(subparsers, **kwargs)
    parser.add_smk_args()


def run(args, smk_args, **kwargs):
    smk = SmkHandler(add_defaults = not args.no_defaults)
    smk.set_cores(args.cores).add_targets('all').set_mem(args.mem_mb).add_arguments(smk_args)
    smk.run()

if __name__ == '__main__':
    sys.exit(run())
