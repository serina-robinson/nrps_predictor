from unittest import TestCase

from minimock import mock, restore

from antismash.common.fasta import write_fasta, read_fasta
from antismash.common import path, subprocessing
from antismash.common.secmet.features import AntismashDomain, PFAMDomain, FeatureLocation
from antismash.modules import nrps_pks

import extract_sig

class TestAngstromGeneration(TestCase):
    def setUp(self):
        self.aligns = read_fasta(path.get_full_path(nrps_pks.__file__, "test", "data", "nrpspred_aligns.fasta"))
        mock("subprocessing.run_muscle_single", returns=self.aligns)

    def tearDown(self):
        restore()

    def test_angstrom(self):
        domain = AntismashDomain(FeatureLocation(1, 2), "test")
        domain.domain_id = "query"
        domain.translation = self.aligns[domain.domain_id].replace("-", "")

        sig = extract_sig.get_34_aa_signature(domain)
        assert sig == "L--SFDASLFEMYLLTGGDRNMYGPTEATMCATW"