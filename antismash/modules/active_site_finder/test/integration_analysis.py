# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

import antismash
from antismash.config import build_config
from antismash.common.record_processing import parse_input_sequence
from antismash.common.test import helpers
from antismash.modules import active_site_finder

class TestAnalyses(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--asf", "--clusterhmmer", "--minimal"],
                            isolated=True,
                            modules=antismash.get_all_modules())

    def test_everything(self):
        datafile = helpers.get_path_to_balhymicin_genbank()
        results = helpers.run_and_regenerate_results_for_module(datafile, active_site_finder, self.options)