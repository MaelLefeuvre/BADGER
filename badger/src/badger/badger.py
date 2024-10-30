import sys
from .logger import Logger, LogLevel
from . import run, rerun, loop_pipeline, unpack, archive, setup


BADGER_MODULES = {
    "setup"        : setup,
    "run"          : run,
    "archive"      : archive,
    "loop-pipeline": loop_pipeline,
    "rerun"        : rerun,
    "unpack"       : unpack
}

"""
Separate badger args from snakemake args into two separate lists.
"""
def _split_args(args = sys.argv, split_flag='--'):
    try:
        switch = args.index(split_flag)
    except ValueError:
        switch = len(args)

    return {'prog':args[0], 'badger': args[1:switch], 'snakemake': args[switch+1:]}

# ----------------------------
import argparse
def _badger_cli():
    # ---- main badger argument parser.
    parser     = argparse.ArgumentParser(prog="badger")
    subparsers = parser.add_subparsers(dest="cmd")

    # ---- Parent parser with default arguments for all modules.
    template_parser = argparse.ArgumentParser(add_help=False)
    template_parser.set_defaults(loglevel=LogLevel.WARN)
    template_parser.add_argument(
      "-q", "--quiet", dest="loglevel", action="store_const", const = LogLevel.QUIET,
      help = "Do not output any logging information"
    )
    template_parser.add_argument(
      "-v", "--verbose", dest="loglevel", action="store_const", const = LogLevel.INFO,
      help = "Output basic logging information"
    )
    template_parser.add_argument(
        "--debug", dest="loglevel", action="store_const", const = LogLevel.DEBUG,
      help = "Output additional logging information."
    )

    # ---- Basic run command
    for module in BADGER_MODULES.values():
        module.setup_cli(subparsers, parents = [template_parser])

    return parser



def main(args = sys.argv):
    parser = _badger_cli()
    args   = _split_args(args)

    if len(args['badger']) < 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    badger_args = parser.parse_args(args=args['badger'])
    
    # ---- setup Logger
    logger = Logger(level=badger_args.loglevel)

    logger.debug(f"Running command: {badger_args.cmd}")
    BADGER_MODULES[badger_args.cmd].run(badger_args, smk_args=args['snakemake'])



if __name__ == '__main__':
    sys.exit(main())
