import os, sys, re, pathlib, glob, shutil
from pathlib import Path
from time import time, sleep, mktime
from datetime import datetime
from types import SimpleNamespace

from .logger import Logger
from .module_components import BadgerCli
from .constants import DEFAULT_CORES, DEFAULT_CONFIG_PATH, RESULTS_DIR

from .module_components import BadgerCli, BadgerSmkCli
from .smk_handler import SmkHandler

from .unpack import run as badger_unpack
from .run import run as badger_run
from .archive import run as badger_archive

TOUCHEABLE_FILE_PATTERNS=[
    "results/meta/pipeline-metadata.yml",
    "results/00-ped-sim/pedigree.def",
    "results/00-ped-sim/???-pedigrees.*",
    "results/00-ped-sim/???-pedigrees-everyone.fam",
    "results/00-ped-sim/???-pedigrees-doppleganged-merged.vcf.gz",
    "results/00-ped-sim/???-pedigrees-doppleganged-merged.vcf.gz.tbi",
    "results/00-ped-sim/???-pedigree-simulations.pdf",
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
    parser.add_argument("-f", "--force", action = "store_true",
      help = "Allow the program to overwrite existing write-protected files.",
      default = False
    )

def tamper_timestamps(files: []):
    #timestamp = mktime(datetime.strptime("1990-01-01", "%Y-%m-%d").timetuple())
    timestamp = time()
    if type(files) is not list:
        files = [files]
    for file in files:
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
    logger.info(f"Searching output archive directory in {DEFAULT_CONFIG_PATH}.")
    output_archive_dir = find_archive_dir(DEFAULT_CONFIG_PATH)
    logger.info(f"Found archive dir: {output_archive_dir}.")

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
        'cores' : args.cores, 'output_dir': "./", 'force': args.force
    }

    # ---- Start a snakemake handler
    smk = SmkHandler(add_defaults = not args.no_defaults)
    smk = smk.set_cores(args.cores).set_mem(args.mem_mb).add_arguments(smk_args)
    smk = smk.set_cores(args.cores).add_targets(["correctKin", "KIN", "GRUPS", "READ", "READv2", "TKGWV2"])
    if args.no_rerun_triggers:
        smk = smk.add_arguments(["--rerun-triggers", "mtime"])

    for run in archived_runs:
        run_id = os.path.basename(run)
        output_run = Path(os.path.join(output_archive_dir, run_id))
        if output_run.exists():
            logger.info(f"{run_id} already found in '{output_archive_dir}'. Skipping...")
            sleep(1)
            continue

        # ---- Unpack previously archived results...    
        logger.info(f"Unpacking {run_id}...")
        unpack_args['archive_dir'] = run
        badger_unpack(args=SimpleNamespace(**unpack_args))

        # ---- Modify the timestamps if user requested it.
        if args.no_rerun_triggers:
            #smk_args.extend(["--rerun-triggers", "mtime"])
            smk = smk.add_arguments(["--rerun-triggers", "mtime"])
            logger.info("Changing timestamps of unpacked files...")
            for pattern in TOUCHEABLE_FILE_PATTERNS:
                matches = glob.glob(pattern, recursive = True)
                logger.debug(f"  - {pattern}: {matches}")
                tamper_timestamps(matches)
            #for smk_dir in [".snakemake/metadata/", ".snakemake/incomplete/", ".snakemake/locks/"]:
            #    shutil.rmtree(smk_dir, ignore_errors=True)

        # ---- Run the archive from this checkpoint
        smk.copy().add_arguments(['--touch', '--keep-going', '--quiet', 'all']).run()
        smk.run()

        # ---- Re-archive these results in the specified archive directory.
        badger_archive(args, smk_args)






if __name__ == '__main__':
    sys.exit(run())
