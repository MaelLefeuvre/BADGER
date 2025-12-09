from os import path
from requests.packages import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

configfile: "config/config.yml"
configfile: "config/netrules.yml"

class AADRDataset:
    _version                = config["aadr-dataset"].lower()
    _server_url             = config["netrules"]["aadr-dataset"]["server-url"]
    _persistent_id          = config["netrules"]["aadr-dataset"]["persistent-id"]
    _dataverse_version_dict = config["netrules"]["aadr-dataset"]["dataverse-version-dict"]

    @staticmethod
    def get_version():
        return config["aadr-dataset"].lower()

    @classmethod
    def get_url(cls):
        try:
            dataverse_version = cls._dataverse_version_dict[AADRDataset.get_version()]
        except KeyError:
            error_msg = "Invalid version number for 'aadr-dataset' parameter. Available versions are:\n"
            for key, value in cls._dataverse_version_dict.items():
                error_msg += f" - {key}  (Dataverse version {value})\n"
            raise RuntimeError(error_msg)

        return (
            f"{cls._server_url}" +
            "/api/access/dataset/:persistentId/versions/" +
            f"{dataverse_version}" +
            f"?persistentId={cls._persistent_id}"
        )

    @classmethod
    def get_path(cls):
        return f"data/Reich-dataset/{cls._version}"

    @classmethod
    def get_output(cls):
        return multiext(f"data/Reich-dataset/{cls._version}/aadr_{cls._version}_1240K_public", ".snp", ".anno", ".geno", ".ind")

    @classmethod
    def get_default_target(cls):
        return f"data/Reich-dataset/{cls._version}/aadr_{cls._version}_1240K_public.snp"

    @classmethod
    def list_available_versions(cls):
        for key in cls._dataverse_version_dict.keys():
            print(f" - {key}")


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
    default      = AADRDataset.get_default_target()
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

