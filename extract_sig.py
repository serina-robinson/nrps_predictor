from collections import defaultdict
import os
import sys
from typing import Any, Dict, List, Set

from antismash.common import path, subprocessing
from antismash.common.secmet import AntismashDomain
from antismash.config import ConfigType
from antismash.modules import nrps_pks
from antismash.modules.nrps_pks.nrps_predictor import build_position_list, read_positions, extract, verify_good_sequence

REF_SEQUENCE = "P0C062_A1"
A34_POSITIONS_FILENAME = path.get_full_path(nrps_pks.__file__, "external", "NRPSPredictor2", "A34positions.txt")
APOSITION_FILENAME = path.get_full_path(nrps_pks.__file__, "external", "NRPSPredictor2", "Apositions.txt")
KNOWN_CODES = path.get_full_path(nrps_pks.__file__, "knowncodes.fasta")
ILLEGAL_CHARS = "!@#$%^&*(){}:\"<>?/.,';][`~1234567890*-+-=_\\|"
ADOMAINS_FILENAME = path.get_full_path(nrps_pks.__file__, "external", "NRPSPredictor2", "A_domains_muscle.fasta")
START_POSITION = 66

def get_34_aa_signature(domain: AntismashDomain) -> str:
    """ Extract 10 / 34 AA NRPS signatures from A domains """
    assert " " not in domain.get_name()
    assert verify_good_sequence(domain.translation)

    # Run muscle and collect sequence positions from file
    alignments = subprocessing.run_muscle_single(domain.get_name(), domain.translation, ADOMAINS_FILENAME)

    domain_alignment = alignments[domain.get_name()]
    reference_alignment = alignments[REF_SEQUENCE]

    positions = read_positions(APOSITION_FILENAME, START_POSITION)
    # Count residues in ref sequence and put positions in list
    poslist = build_position_list(positions, reference_alignment)

    # Extract positions from query sequence
    query_sig_seq = extract(domain_alignment, poslist)
    # Add fixed lysine 517
    query_sig_seq += "K"

    # repeat with 34 AA codes
    angpositions = read_positions(A34_POSITIONS_FILENAME, START_POSITION)
    poslist = build_position_list(angpositions, reference_alignment)

    return extract(domain_alignment, poslist)

if __name__ == '__main__':

	# signatures = [get_34_aa_signature(a_domain) for a_domain in a_domains]
