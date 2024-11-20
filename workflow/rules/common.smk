from os import path
from requests.packages import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

configfile: "config/config.yml"
configfile: "config/netrules.yml"

class ReferenceGenome:
    _name = config["reference-genome"].lower()
    _dict = config["netrules"]["reference-genomes"]
    @staticmethod
    def get_name():
        return config["reference-genome"].lower()

    @classmethod
    def get_url(cls):
        return cls._dict[ReferenceGenome.get_name()]["url"]

    @classmethod
    def get_path(cls):
        return cls._dict[ReferenceGenome.get_name()]["path"]

    @classmethod
    def list_available_references(cls):
        for key in cls._dict.keys():
            print(f" - {key}")

# ------------------------------------------------------------------------------------------------------------------- #
# ----- Default files when not supplied by the user

def get_snp_targets(ext=".snp"):
    user_defined = config["kinship"]["targets"]
    default      = config["netrules"]["aadr-1240k"]["default-targets"]
    return path.splitext(user_defined or default)[0] + ext

# ------------------------------------------------------------------------------------------------------------------- #
# ---- Generic / Utility rules and functions.

def get_generations():
    n = config['ped-sim']['replicates']
    return [f"ped{i}" for i in range(1, n+1)]

def get_samples_ids(wildcards):
    """
    Get each pedigree sample's id based on each unique pedigree sample id + the
    number of pedigree replicates. This list is sorted to match UNIX's 
    """
    # Run through the initial samples files and extract pedigree ids 
    with checkpoints.get_samples.get().output[0].open() as f:
        samples = str.split(f.readline().replace('\n', ''), '\t')
        ids     = set([sample.split('_')[1] for sample in samples])
        return sorted(expand("{generation}_{ids}",
            ids=ids, generation="{generation}"), key=str.casefold
        )

def get_all_samples_ids(wildcards):
    n = config['ped-sim']['replicates']
    return expand(get_samples_ids(wildcards), generation=[f'ped{i}' for i in range(1, n+1)])


"""Simple decorator to exclude ped-sim user-selected samples from the workflow"""
def exclude_samples(function):
    def _decorated(wildcards):
        excluded_samples = config['kinship']['exclude-samples']
        x = list(filter(lambda x: not any(id in x for id in excluded_samples), function(wildcards)))
        return x
    return _decorated

@exclude_samples
def get_samples_ids_filtered(wildcards):
    return get_samples_ids(wildcards)

@exclude_samples
def get_all_samples_ids_filtered(wildcards):
    return get_all_samples_ids(wildcards)

