from .logger import Logger
from .constants import DEFAULT_CORES
from os import path
from typing import Iterator
import sys, os, re, glob, tarfile, lzma, shutil
import pysam

from .module_components import BadgerCli

class UnpackCli(BadgerCli):
    def help(self): return "Unpack an archive of BADGER results"
    def description(self): return "Unpack an archive of BADGER results in a specified directory"
    def usage(self): return "%(prog)s [filetype, ...] [options]"


_ALL_FILETYPES_KW  = "all"
_ALLOWED_FILETYPES = [
    _ALL_FILETYPES_KW, "uncram", "metadata", "contaminants",
    "pedigree", "variants", "misincorporation", "correctKin", 
    "GRUPS", "KIN", "READ", "READv2", "TKGWV2"
]

def setup_cli(subparsers, **kwargs) -> None:
    parser = UnpackCli(subparsers, **kwargs)
    parser.add_argument("-c", "--cores", action='store',
      help = "Set the maximum number of CPU cores for snakemake (Default: %(default)s)", 
      default = DEFAULT_CORES
    )
    parser.add_argument("-a", "--archive-dir", action = "store",
      help = "Specify the location of the archive directory to unpack",
      required = True
    )
    parser.add_argument("-f", "--force", action = "store_true",
      help = "Allow the program to overwrite existing write-protected files.",
      default = False
    )
    parser.add_argument("-o", "--output-dir", action = "store",
      help = "Specify the output directory where the archive should be unpacked.",
      default = "."
    )
    parser.add_argument('filetype', nargs='+',
        choices=_ALLOWED_FILETYPES,
        default=[_ALL_FILETYPES_KW], 
        action = "store",
        help='Specify which type of file(s) within the archive should be extracted.'
    )


def _file_search(target_dir: str, pattern: str = "**", recursive: bool = False) -> Iterator[str]:
    Logger().debug(f"Searching files matching pattern '{pattern}' in target_directory '{target_dir}' (recursive={recursive})")
    regex = re.compile(pattern)
    return (f for f in glob.iglob(f"{target_dir}/**", recursive=recursive) if re.search(regex, f))


def _build_output_path(file: str, archive_dir: str, output_dir: str) -> None:
    return path.join(
        output_dir,
        path.relpath(file, archive_dir)
    )


def _try_mkdirs(dirs: str, **kwargs) -> None:
    dirs = dirs if isinstance(dirs, list) else [dirs]
    for dir in dirs:
        try:
            Logger().debug(f"Creating target directory '{dir}'")
            os.makedirs(dir, exist_ok=True)
        except PermissionError as e:
            Logger().error(f"Cannot create target directory '{dir}' ({type(e).__name__})")
            sys.exit(1)
    return

def _try_chdir(dir: str, **kwargs) -> None:
    try:
        Logger().debug(f"Changing into directory {dir}")
        os.chdir(dir)
    except (PermissionError, FileNotFoundError) as e:
        Logger().error(f"Cannot change into {dir} ({type(e).__name__})")
        sys.exit(1)

def uncram(cramfile: str, target_dir: str, cores: int = 1) -> None:
    logger = Logger()
    # ---- Extract the pedigree number from the name of the cramfile
    pattern = re.compile('ped([0-9]+)')
    rep = re.search(pattern, path.basename(cramfile)).group(1)
    logger.debug(f"Matched pedigree n°{rep}")

    # ---- Change CWD to target directory
    main_wd=os.getcwd()
    _try_chdir(target_dir)

    # ---- Decompress and split .cram file
    logger.info(f"Decompressing samples from {cramfile}...") 
    try:
        split_args = ['-@', str(cores), '-f', '%!.bam', '--output-fmt', 'BAM', cramfile]
        logger.debug(f"Running 'samtools split {' '.join(split_args)}'")
        pysam.split(*split_args)
    except pysam.SamtoolsError as e:
        logger.error(f"Failed to run samtools split ({type(e).__name__})")
        pysam.sort.get_messages()
        raise e

    # ---- Put each unpacked bam in its own subdirectory
    for bam in glob.iglob("./**.bam", recursive=False):
        bam_dir = path.splitext(bam)[0]
        _try_mkdirs(bam_dir)
        try:
            logger.debug(f"Moving {bam} in {target_dir}/{bam_dir}")
            os.chmod(bam, 0o444)
            os.rename(bam, f"{bam_dir}/{bam_dir}.srt.rmdup.bam")
        except PermissionError as e:
            logger.error(f"Failed to move {bam} in its subdirectory (PermissionError)")
    _try_chdir(main_wd)

def uncram_all(archive_dir: str, output_dir: str, cores: int = 1, **kwargs):
    logger = Logger()
    for cram in _file_search(archive_dir, r".cram$", recursive=True):
        target_dir = path.dirname(_build_output_path(cram, archive_dir, output_dir))
        _try_mkdirs(target_dir) 
        try:
            logger.debug(f"Attempting to uncram '{cram}'")
            uncram(cram, target_dir, cores=cores)
        except pysam.SamtoolsError as e:
            logger.error(f"Failed to uncram {cram}")
            sys.exit(1)

def copy_one_file(archive_dir: str, output_dir: str, pattern: str, force: bool = False, **kwargs):
    logger = Logger()
    file = list(_file_search(archive_dir, pattern=pattern, recursive=True))
    if len(file) > 1:
        logger.error(
            f"Found multiple candidate files matching pattern '{pattern}' in {archive_dir}. "
            f"Will only use first occurence: {file[0]}"
        )
        sys.exit(1)
    try:
        file = file[0]
    except IndexError:
        logger.error(f"Failed to find any file matching pattern '{pattern}' in {archive_dir}")
        sys.exit(1)

    target_dir=path.dirname(_build_output_path(file, archive_dir, output_dir))
    _try_mkdirs(target_dir)
    logger.info(f"Copying {file} in {output_dir}")
    output_file = path.join(target_dir, path.basename(file))
    try:
        if path.exists(output_file) and force: 
            logger.debug(f"{output_file} already exists, while overwrite was requested. Attempting to change file permissions.")
            os.chmod(output_file, 0o644)
        _output = shutil.copy2(file, output_file)
        shutil.copystat(file, output_file)
    except PermissionError as e :
        logger.error(f"Failed to write {file} in {target_dir} ({type(e).__name__})")
        if path.exists(output_file):
            logger.error(
              "Some output files appear to already exist in target directory. "
              "Use --force argument to force overwrite."
            )

def copy_files(archive_dir: str, output_dir: str, pattern: str, force: bool = False, **kwargs) -> None:
    for file in _file_search(archive_dir, pattern=pattern, recursive = True):
        Logger().info(f"Copying {file} in {output_dir}")
        file_pattern=_build_output_path(file, archive_dir,"") + "$"
        copy_one_file(archive_dir, output_dir, file_pattern, force=force)


def untar_one_archive(archive_dir: str, output_dir: str, pattern: str, **kwargs) -> None:
    logger = Logger()
    file = list(_file_search(archive_dir, pattern, recursive=True))
    if len(file) > 1:
        logger.error(
            f"Found multiple candidate tar archives matching pattern '{pattern}' in {archive_dir}. "
            f"Will only use first occurence: {file[0]}"
        )
        sys.exit(1)
    try:
        file = file[0]
    except IndexError:
        logger.error(f"Failed to find any archive matching pattern '{pattern}' in {archive_dir}")
        sys.exit(1)

    _try_mkdirs(output_dir)
    logger.info(f"Decompressing {file} in {output_dir}")

    tar = tarfile.open(file, mode='r')
    tar.extractall(output_dir, filter = "data")


def untar_archive(archive_dir: str, output_dir: str, pattern:str , **kwargs) -> None:
    for tarball in _file_search(archive_dir, pattern=pattern, recursive = True):
        Logger().info(f"Decompressing {tarball} in {output_dir}")
        tar_pattern=_build_output_path(tarball, archive_dir,"")
        untar_one_archive(archive_dir, output_dir, tar_pattern)
    return

def unpack_lzma(archive_dir: str, output_dir: str , pattern: str, **kwargs) -> None:
    logger = Logger()
    xz = list(_file_search(archive_dir, pattern, recursive=True))
    if len(xz) > 1:
        logger.error(
            f"Found multiple candidate LZMA-compressed files matching pattern '{pattern}' in {archive_dir}. "
            f"Will only use first occurence: {xz[0]}"
        )
        sys.exit(1)
    try:
        xz = xz[0]
    except IndexError:
        logger.error(f"Failed to find any LZMA-compressed file matching pattern '{pattern}' in {archive_dir}")
        sys.exit(1)

    target_dir=path.dirname(_build_output_path(xz, archive_dir, output_dir))
    _try_mkdirs(target_dir)
    logger.info(f"Copying {xz} in {output_dir}")
    output_file = path.join(target_dir, path.basename(xz))
    try:
        with lzma.open(xz, "rb") as stream:
            with open(output_file, "wb") as dest:
                shutil.copyfileobj(stream, dest)
    except Exception as e:
        logger.error(
            "Failed to decompress {xz} contents in {output_file}.\n"
            "Error message: {e}"
        )
        sys.exit(1)

def run(args, **kwargs) -> None:
    logger = Logger()
    if _ALL_FILETYPES_KW in args.filetype:
        args.filetype = [filetype for filetype  in _ALLOWED_FILETYPES if filetype != _ALL_FILETYPES_KW] 

    for request in args.filetype:
        logger.hr()
        logger.info(f"Unpacking {request}")
        match request:
            case "uncram":
                uncram_all(**vars(args))
            case "metadata" | "contaminants":
                pattern = {
                  "metadata": r"meta/pipeline-metadata.yml$",
                  "contaminants": r"contaminants.tsv$"
                }[request]
                #copy_one_file(args.archive_dir, args.output_dir, pattern, args.force, logger)
                copy_one_file(**vars(args), pattern=pattern)

            case "pedigree":
                pattern = r"results/00-ped-sim.tar.xz"
                untar_one_archive(**vars(args), pattern=pattern)

            case "variants":
                pattern = r"variants-intersect-[A-Z]{3}_maf[.0-9]+.ucscbed.xz$"
                unpack_lzma(**vars(args), pattern=pattern)

            case "misincorporation" | "TKGWV2":
                logger.info(f"Unpacking {request}")
                pattern = {
                  "misincorporation": r"misincorporation.txt$",
                  "TKGWV2": r"ped[0-9]+-TKGWV2_Results.txt$"
                }[request]
                copy_files(**vars(args), pattern=pattern)

            case "correctKin" | "GRUPS" | "KIN" | "READ" | "READv2":
                pattern = {
                  "correctKin": r"correctKin/ped[0-9]+[.]tar[.]xz$",
                  "GRUPS":      r"GRUPS/ped[0-9]+[.]tar[.]xz$",
                  "KIN":        r"KIN-results.tar.xz$",
                  "READ":       r"READ/ped[0-9]+[.]tar[.]xz$",
                  "READv2":     r"READv2/ped[0-9]+[.]tar[.]xz$",
                }[request]
                untar_archive(**vars(args), pattern=pattern)

                




