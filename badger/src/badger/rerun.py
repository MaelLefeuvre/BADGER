import os, sys, re, pathlib, glob
from pathlib import Path
from time import sleep, mktime
from datetime import datetime
from types import SimpleNamespace

from .logger import Logger
from .module_components import BadgerCli
from .constants import DEFAULT_CORES, DEFAULT_CONFIG_PATH, RESULTS_DIR

from .module_components import BadgerCli, BadgerSmkCli

from .unpack import run as badger_unpack
from .run import run as badger_run
from .archive import run as badger_archive

TOUCHEABLE_FILE_PATTERNS=[
    "results/meta/pipeline-metadata.yml",
    "results/00-ped-sim/pedigree.def",
    "results/00-ped-sim/???-pedigrees.*",
    "results/00-ped-sim/???-pedigrees-everyone.fam",
    "results/meta/pipeline-metadata.yml",
    "results/00-ped-sim/???-pedigrees.*",
    "results/00-ped-sim/???-pedigrees-everyone.fam",
    "results/00-ped-sim/???-pedigrees-doppleganged-merged.vcf.gz",
    "results/00-ped-sim/sample_names.tsv",
    "results/01-gargammel/contaminants/contaminants.tsv",
    "results/02-preprocess/05-dedup/**/*.bam",
    "results/02-preprocess/06-mapdamage/**/*"
]

class RerunCli(BadgerCli, BadgerSmkCli):
    def help(self): return "Re-apply a BADGER run from an archive"

def setup_cli(subparsers, **kwargs):
    parser = RerunCli(subparsers, **kwargs)
    parser.add_smk_args()
    parser.add_argument('-i', "--input-archive-dir", action='store',
      required = True,
      help = "Input archive directory that you wish to rerun.",
    )
    parser.add_argument('-o', "--output-archive-dir", action='store',
      help = "Output directory where results should be archived once the rerun is completed.",
    )
    parser.add_argument('-n', "--no-rerun-triggers", action='store_true',
      help = "Force snakemake to never apply rule rerun trigger.",
    )

def tamper_timestamps(file: []):
    timestamp = time.mktime(datetime.strptime("2024-01-01", "%Y-%m-%d").timetuple())
    os.utime(file, (timestamp, timestamp))

def find_archive_dir(config_file, keyword = r"archive-dir:"):
    config = pathlib.Path(config_file)
    if not config.exists():
        Logger().error(f"Cannot find config file {config_file}")
        sys.exit(1)

    matched_line = None
    regex  = re.compile(keyword)
    with open(config, 'r') as file:
        for line in file:
            if re.search(regex, line):
                matched_line = line
                break

    if matched_line is None:
        Logger().error(f"Cannot find any line matching keyword {keyword} in {config_file}")
        sys.exit(1)

    return matched_line.strip().split(':')[1].split('#')[0].replace(' ', '').strip('"')

def find_archived_runs(archive_directory):
    return glob.glob(f'{archive_directory}/run-*')

def run(args, smk_args):

    logger = Logger()
    if not args.output_archive_dir:
        logger.info(
            "No output archive directory provided through the command line. "
            f"Attempting find one specified within {DEFAULT_CONFIG_PATH}."
        )
        args.output_archive_dir = find_archive_dir(DEFAULT_CONFIG_PATH)

    try:
        Path(args.output_archive_dir).mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Failed to icreate {args.output_archive_dir} ({e})")
        sys.exit(1)

    logger.info(f"Looking through input archive directory: {args.input_archive_dir}")
    if not Path(args.input_archive_dir).exists():
        Logger().error(f"Directory {args.input_archive_dir} does not appear to exist.")
        sys.exit(1)

    archived_runs = find_archived_runs(args.input_archive_dir)
    archived_runs.sort()
    logger.info(f"found {len(archived_runs)} within the archive directory.")

    logger.info(f"Starting reruns...")
    logger.hr()

    unpack_args={
        'filetype': ['metadata', 'contaminants', 'pedigree', 'uncram', 'misincorporation'],
        'cores' : args.cores, 'output_dir': RESULTS_DIR, 'force': False
    }


    for run in archived_runs:
        run_id = os.path.basename(run)
        output_run = Path(os.path.join(args.output_archive_dir, run_id))
        print(f"path: {output_run}")
        if output_run.exists():
            logger.info(f"{run_id} already found in '{args.output_archive_dir}'. Skipping...")
            sleep(1)
            continue

        # ---- Unpack previously archived results...    
        logger.info(f"Unpacking {run_id}...")
        unpack_args['archive_dir'] = run
        badger_unpack(args=SimpleNamespace(**unpack_args))

        # ---- Modify the timestamps if user requested it.
        if args.no_rerun_triggers:
            logger.info("Changing timestamps of unpacked files...")
            tamper_timestamps(
                [match for match in glob.glob(pattern) for pattern in TOUCHEABLE_FILE_PATTERNS]
            )

        # ---- Run the archive from this checkpoint
        badger_run(args, smk_args)

        # ---- Re-archive these results in the specified archive directory.
        badger_archive(args, smk_args)






if __name__ == '__main__':
    sys.exit(run())
