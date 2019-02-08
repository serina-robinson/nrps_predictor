import sys
from typing import Dict, Iterable, List
import warnings
import re
from tempfile import NamedTemporaryFile

from antismash.common.fasta import write_fasta, read_fasta
from antismash.common.subprocessing import execute
import alignment

def run_asm(queryfa: Dict[str, str], stachfa: Dict[str, str], seedfa: Dict[str, str]) -> [str, int]:
    """ Active site motif (ASM) substrate prediction
    Arguments:
        queryfa: seq id to seq dictionary
        stachfa: seq name to seq for stachelhaus codes
        seedfa: seq name to seq for seed alignment for stachelhaus code extraction
    Returns:                                                                                                                            substrate specificity prediction (str)
        number of identical matches (int)
    """ 
    ## ASM settings
    gapopen = 3.4
    properlen = 117 ## true length
    grsAcode = {4:1,5:1,8:1,47:1,68:1,70:1,91:1,99:1,100:1} ## positions in grsA for code
    ## Alignment
    toalign = {**queryfa, **seedfa}
    aligned2seed = subprocessing.mafft_sandpuma_asm(toalign, gapopen)
    ## Loop through alignment to find new positions for code
    qname = next(iter(queryfa))
    pos = 0
    newcode = []
    for p, val in enumerate(aligned2seed['phe_grsA']):
        if(val=='-'):
            continue
        else:
            pos += 1
            if(pos in grsAcode):
                newcode.append(p)
    ## Extract codes
    extractedcode = {}
    for seqname in aligned2seed:
        code = ''
        for n in newcode:
            code = code + aligned2seed[seqname][n]
            extractedcode[seqname] = code
    ## Error checking
    truth = {'phe_grsA':'DAWTIAAIC', 'asp_stfA-B2':'DLTKVGHIG','orn_grsB3':'DVGEIGSID','val_cssA9':'DAWMFAAVL'}
    for seqname in extractedcode:
        if seqname == qname:
            continue
        else:
            if extractedcode[seqname] != truth[seqname]:
                #print("\t".join([seqname, extractedcode[seqname], truth[seqname]]))
                return('no_call','0') ## Issue with the alignment
    ## Score each
    scores = {}
    for sname in stachfa:
        match = 0
        split_id = re.split("_+", sname)
        if re.match(r"\|", split_id[-1]) is not None: 
            spec = re.split("|", split_id[-1])
        else:
            spec = [split_id[-1]]
        for p, val in enumerate(stachfa[sname]):
            if val == extractedcode[qname][p]:
                match += 1
        if str(match) in scores:
            for s in spec:
                if s in scores[str(match)]:
                    scores[str(match)][s] += 1
                else:
                    scores[str(match)][s] = 1
        else:
            scores[str(match)] = {}
            for s in spec:
                scores[str(match)][s] = 1
    ## Dereplicate and return spec predictions
    for i in range(0,10):
        m = str(9-i)
        if m in scores:
            seen = {}
            for s in scores[m]:
                if s.count('|') > 0:
                    for ss in s.split('|'):
                        seen[ss] = 1
                else:
                    seen[s] = 1
            return('|'.join(sorted(seen)), m )
    return('no_call','0')


def main(queryfa, stachfa, seedfa):
    run_asm(queryfa, stachfa, seedfa)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Not enough arguments")

    if len(sys.argv) == 3:
        main(read_fasta(sys.argv[1]), 'data/fullset0_smiles.stach.faa', read_fasta(sys.argv[3]))



