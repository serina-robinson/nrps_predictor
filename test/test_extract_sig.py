import unittest 

from minimock import mock, restore

from antismash.common.fasta import write_fasta, read_fasta
from antismash.common import path, subprocessing
from antismash.common.secmet.features import AntismashDomain, PFAMDomain, FeatureLocation
from antismash.modules import nrps_pks
from bin import extract_sig

class TestConversion(unittest.TestCase):
    def test_conversion(self):
        domain = AntismashDomain(FeatureLocation(1, 3, 1), tool="test")
        domain.domain_subtype = "subtest"
        domain.specificity = ["a", "c", "f"]
        domain.asf.add("first")
        domain.asf.add("second")
        assert domain.tool == "test"
        assert domain.created_by_antismash

        bio = domain.to_biopython()
        assert len(bio) == 1
        assert bio[0].qualifiers["aSTool"] == ["test"]
        assert bio[0].qualifiers["tool"] == ["antismash"]
        new_domain = AntismashDomain.from_biopython(bio[0])
        assert new_domain.domain_subtype == domain.domain_subtype == "subtest"
        assert new_domain.specificity == domain.specificity == ["a", "c", "f"]
        assert new_domain.asf.hits == domain.asf.hits
        assert new_domain.asf.hits == ["first", "second"]
        assert new_domain.tool == domain.tool == "test"
        assert new_domain.created_by_antismash

class TestAngstromGeneration(unittest.TestCase):
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

