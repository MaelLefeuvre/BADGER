import inspect
from abc import ABC, abstractmethod, abstractproperty
from os.path import basename, splitext
from typing import Self

from .constants import DEFAULT_CORES, DEFAULT_MEM_MB

class BadgerCli(ABC):
    _parser = None

    def __init__(self, subparsers, **kwargs):
        self._parser = subparsers.add_parser(
            name        = self.name(),
            help        = self.help(),
            description = self.description(),
            usage       = self.usage(),
            **kwargs
        )

    def description(self) -> str:
        return self.help()

    def usage(self) -> str:
        return "%(prog)s [options] -- [snakemake options]"

    def name(self) -> str:
        return splitext(basename(inspect.getfile(self.__class__)))[0].replace("_", "-")

    def add_argument(self, *args, **kwargs) -> Self :
        self._parser.add_argument(*args, **kwargs)
        return self

    @abstractmethod
    def help(self) -> str:
        pass


class BadgerSmkCli(ABC):
    def add_smk_args(self):
        self._parser.add_argument("-c", "--cores", action='store',
            help = "Set the maximum number of CPU cores for snakemake (Default: %(default)s)", 
            default = DEFAULT_CORES
        )
        self._parser.add_argument("-m", "--mem-mb", action='store',
            help = "Set the maximum amount of Preemptible memory (in MB) for snakemake (Default: %(default)s)",
            default = DEFAULT_MEM_MB
        )
        self._parser.add_argument("--no-defaults", action='store_true',
            help = f"Don't provide snakemake with default arguments."
        )
