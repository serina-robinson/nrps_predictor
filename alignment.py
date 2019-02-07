import sys
from typing import Dict, Iterable, List
import warnings
import re
from tempfile import NamedTemporaryFile

from antismash.common.fasta import write_fasta, read_fasta
from antismash.common.subprocessing import execute

def run_mafft(toalign: Dict[str, str], gapopen: float = 3.4) -> Dict[str, str]:
    """ Runs mafft for initial alignment for stachelhaus code extraction
        Arguments:
            toalign: fasta Dictionary of seq names (str) to seqs (str)
            gapopen: gapopen penalty
        Returns:
            aligned fasta Dictionary of all seq names (str) to seqs (str)
    """
    if len(toalign) < 2:
        raise ValueError("Requires at least two sequences")

    if is_invalid_sequence(toalign.values()):
        raise ValueError("Invalid sequence")   

    fnames, fseqs = [], []
    for k, v in toalign.items():
        fnames.append(k)
        fseqs.append(v)
    with NamedTemporaryFile(mode="w+") as query_unaligned:
        write_fasta(fnames, fseqs, query_unaligned.name)
        with NamedTemporaryFile(mode="w+") as all_aligned:
            mafft_result = execute(["mafft", "--quiet",
                                     "--op", str(gapopen),
                                     "--namelength", "60",
                                     query_unaligned.name],
                                   stdout=all_aligned)
            if not mafft_result.successful():
                raise RuntimeError("mafft returned %d: %r while performing ASM seed alignment" % (
                    mafft_result.return_code, mafft_result.stderr.replace("\n", "") ))
            return read_fasta(all_aligned.name)


def verify_good_sequence(sequence: str) -> bool:
    """ Ensures a sequence is valid """
    for char in ILLEGAL_CHARS:
        if char in sequence:
            return False
    return True


