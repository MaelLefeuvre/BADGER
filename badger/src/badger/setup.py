from .logger import Logger, LogLevel
from .constants import DEFAULT_CORES, DEFAULT_MEM_MB
import snakemake

from .smk_handler import SmkHandler
from .module_components import BadgerCli, BadgerSmkCli


class RunCli(BadgerCli, BadgerSmkCli):
    def help(self): return "Pre-configure all conda environments and download all datasets required to run BADGER"


def setup_cli(subparsers, **kwargs):
    parser = RunCli(subparsers, **kwargs)
    parser.add_smk_args()
    parser.add_argument("--no-create-envs", action='store_true',
      help = "Do not attempt to pre-configure conda environments", 
    )
    parser.add_argument("--no-fetch-data", action='store_true',
      help = "Do not attempt to download datasets", 
    )

def run(args, smk_args, **kwargs):
    logger = Logger()
    smk = SmkHandler(add_defaults = not args.no_defaults)
    smk.set_cores(args.cores).set_mem(args.mem_mb).add_arguments(smk_args)
 
    if args.no_fetch_data is False:
        logger.info("Fetching datasets...")
        logger.hr()
        try:
            smk_fetch_data = smk.copy().add_targets(['fetch_data'])
            smk_fetch_data.run()
        except Exception as e:
            logger.error("Failed to download the required datasets form BADGER.")
            raise e
    print("hey")
    if args.no_create_envs is False:
        logger.info("Creating all conda environments...")
        logger.hr()

        try:
            smk_create_conda_envs = smk.copy().add_targets(['all', 'archive']).add_arguments(['--conda-create-envs-only'])
            smk_create_conda_envs.run()
        except Exception as e:
            logger.error("Failed to pre-configure the required conda environments")
            raise e


if __name__ == '__main__':
    sys.exit(run())
