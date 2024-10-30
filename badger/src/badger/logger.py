import sys, os
from dataclasses import dataclass, field
from colorama import init as colorama_init
from colorama import Fore, Style
from enum import Enum
import inspect

def _singleton(cls):
    instances = {}
    def wrapper(*args, **kwargs):
        if cls not in instances:
          instances[cls] = cls(*args, **kwargs)
        return instances[cls]
    return wrapper

@dataclass
class LogLevelData:
    priority: int
    ansi    : str
    tail    : bool = field(repr=False, default=True)

class LogLevel(LogLevelData, Enum):
    QUIET = 0, Fore.BLACK  + Style.DIM
    ERROR = 1, Fore.RED    + Style.DIM
    WARN  = 2, Fore.YELLOW + Style.BRIGHT,
    INFO  = 3, Fore.GREEN  + Style.BRIGHT,
    DEBUG = 4, Fore.BLUE   + Style.BRIGHT,
    def __eq__(self, other):
        return self.priority == other.priority
    
    def __lt__(self, other):
        return self.priority < other.priority

    def __le__(self, other):
        return self.__lt__(other) or self.__eq__(other) 


@_singleton
class Logger:
    def __init__(self, level: LogLevel = LogLevel.INFO):
        self.level = level
        self.base_depth = len(inspect.stack())
        colorama_init()
        
    def hr(self, level: LogLevel = LogLevel.INFO):
        if level > self.level:
            return
        width: int = os.get_terminal_size().columns if sys.stdout.isatty() else 80
        print("\u2500" * width)

    def log(self, msg, file=sys.stdout, level: LogLevel =LogLevel.INFO, indent = False):
        if level > self.level:
            return
        color=level.ansi
        depth = len(inspect.stack()) - self.base_depth - 1 if indent else 0
        print(f"{'  ' * depth}[{color}{level.name}{Style.RESET_ALL}] {msg}", file=file)
        
    def error(self, msg, file=sys.stderr):
        self.log(msg, file=file, level=LogLevel.ERROR)

    def warn(self, msg, file=sys.stderr):
        self.log(msg, file=file, level=LogLevel.WARN)

    def info(self, msg, file=sys.stdout):
        self.log(msg, file=file, level=LogLevel.INFO)

    def debug(self, msg, file=sys.stdout):
        self.log(msg, file=file, level=LogLevel.DEBUG)

    def set_level(self, level: LogLevel):
        self.level = level
