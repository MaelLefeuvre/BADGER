#!/usr/bin/env python

import os, sys, time, yaml, requests
import pytest
from abc import ABC, abstractmethod

NETRULES_YAML = os.path.join(os.path.dirname(__file__), '../../config/netrules.yml')

def load_netrules(config_path="config/netrules.yml"):
    with open(config_path) as stream:
        try:
            netrules =yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc, file=sys.stderr)
    return(netrules['netrules'])

class BadgerURL(ABC):
    def __init__(self):
        self.config_yaml = load_netrules(NETRULES_YAML)
        self.config_keys = self._CONFIG_KEYS

    @abstractmethod
    def parse_urls(self) -> [str]:
        ...

    def urls_are_accessible(self) -> bool:
        cls_name = f"{self.__class__.__name__}"
        for url in self.parse_urls():
            if not (url.startswith('https://') or url.startswith('http://')):
                url = 'http://' + url
            resp = requests.head(url, allow_redirects = False)
            time.sleep(0.1)
            resp.raise_for_status()
            print(f"OK ({cls_name}): {url}")
        return True

class G1KPhase3(BadgerURL):
    _CONFIG_KEYS = ["1000g-phase3-url"]
    def parse_urls(self):
        return [self.config_yaml[self.config_keys[0]]]

class HapMapII(BadgerURL):
    _CONFIG_KEYS = ["hapmapII"]
    def parse_urls(self):
        return [self.config_yaml[self.config_keys[0]]]

class AADR(BadgerURL):
    _CONFIG_KEYS = ["aadr-dataset", "server-url", "persistent-id", "dataverse-version-dict"]
    def parse_urls(self):
        config = self.config_yaml[self.config_keys[0]]
        server = config[self.config_keys[1]]
        persistent_id = config[self.config_keys[2]]

        output = []
        for version in config[self.config_keys[3]].values():
            url = f"{server}/api/access/dataset/:persistentId/versions/{version}?persistentId={persistent_id}"
            output.append(url)
        return output

class RefGen(BadgerURL):
    _CONFIG_KEYS=["reference-genomes"]
    def parse_urls(self):
        config = self.config_yaml["reference-genomes"]
        return [ values['url'] for values in config.values()]

class RefinedGeneticMap(BadgerURL):
    _CONFIG_KEYS=["ped-sim", "refined-genetic-map-url"]
    def parse_urls(self):
        return [self.config_yaml[self.config_keys[0]][self.config_keys[1]]]

class InterferenceMap(BadgerURL):
    _CONFIG_KEYS=["ped-sim", "interference-map-url"]
    def parse_urls(self):
        return [self.config_yaml[self.config_keys[0]][self.config_keys[1]]]

class TKGWV2(BadgerURL):
    _CONFIG_KEYS=["TKGWV2", "support-files-url"]
    def parse_urls(self):
        return [self.config_yaml[self.config_keys[0]][self.config_keys[1]]]

@pytest.mark.urls
def test_badger_urls():
    config_yaml = load_netrules()
    for obj in [G1KPhase3, HapMapII, AADR, RefGen, RefinedGeneticMap, InterferenceMap, TKGWV2]:
        obj_instance = obj()
        _return = obj_instance.urls_are_accessible()

if __name__ == '__main__':
    test_badger_urls()

