import os

from .module_components import BadgerCli
from .constants import DEFAULT_CORES

from .module_components import BadgerCli, BadgerSmkCli


class RerunCli(BadgerCli, BadgerSmkCli):
    def help(self): return "Re-apply a BADGER run from an archive"

def setup_cli(subparsers, **kwargs):
    parser = RerunCli(subparsers, **kwargs)
    parser.add_smk_args()

def run(args, smk_args):
    print("Hello from rerun")
    print(args)

if __name__ == '__main__':
    sys.exit(run())
