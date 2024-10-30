import os, psutil

DEFAULT_CORES  = int(os.cpu_count() / 2)
DEFAULT_MEM_MB = int(psutil.virtual_memory()[0] / 1_000_000 / 2)

RESULTS_DIR = "results"