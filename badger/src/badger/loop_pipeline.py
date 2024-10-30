from .module_components import BadgerCli, BadgerSmkCli
from .logger import Logger
from .constants import DEFAULT_CORES, DEFAULT_MEM_MB, RESULTS_DIR

import shutil, sys
from os import path

class ModuleCli(BadgerCli, BadgerSmkCli):
    def help(self): return "Run n iterations of BADGER"
    def description(self): return "Run and archive n iterations of BADGER using a set configuration"


def setup_cli(subparsers, **kwargs):
    parser = ModuleCli(subparsers, **kwargs)
    parser.add_smk_args()
    parser.add_argument("-i", "--iterations", required=True,
        help = "Specify the number of iterations"
    )
    parser.add_argument('--no-remove',
        help = "Do not remove the 'results' directory on the first run. This can be useful a previous attempt resulted in a snakemake failure.",
        action = "store_true",
        default = False
    )



from .run import run as badger_run
from .archive import run as badger_archive

def run(args, smk_args, **kwargs):


    logger = Logger()
    N = int(args.iterations)
    for n in range(N):
        prefix_i = lambda msg: f"(Iteration {n+1}/{N}) {msg}"
        logger.info(prefix_i("Starting loop."))
        logger.hr()
        if n == 1 and args.no_remove:
            logger.info(prefix_("Skipping removal of results for the first run, as --no-remove was requested"))
        else:
            logger.info(prefix_i("Removing previous results"))
            shutil.rmtree(RESULTS_DIR, ignore_errors = True)
            if path.isdir(RESULTS_DIR):
                logger.error(prefix_("Failed to remove directory of previous run..."))
                sys.exit(1)
        
        try:
            logger.info(prefix_i("Running snakemake pipeline..."))
            badger_run(args, smk_args, **kwargs)
            logger.hr()
            logger.info(prefix_i("Archiving results..."))
            badger_archive(args, smk_args, **kwargs)
        except Exception as e:
            logger.error(prefix_i("Failed run (See above for details...)."))
            sys.exit(1)


if __name__ == '__main__':
    sys.exit(run())
