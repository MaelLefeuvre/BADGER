def get_generations():
    n = config['ped-sim']['replicates']
    return [f"ped{i}" for i in range(1, n+1)]

# ------------------------------------------------------------------------------------------------------------------- #
# ---- Generic / Utility rules and functions.

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

