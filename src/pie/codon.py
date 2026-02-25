"""Nei-Gojobori codon lookup tables.

Precomputed at import time:
- CODON_TO_INDEX / INDEX_TO_CODON: codon <-> integer (ACGT lexicographic)
- AMINO_ACID: amino acid for each codon index (stop = "*")
- N_SITES / S_SITES: fractional nonsynonymous/synonymous site counts (64, 3)
- N_DIFFS / S_DIFFS: pairwise N/S differences between all codon pairs (64, 64)
"""

from itertools import permutations

import numpy as np

# ---------------------------------------------------------------------------
# 1. Codon indexing — ACGT lexicographic order
# ---------------------------------------------------------------------------
_BASES = "ACGT"

CODON_TO_INDEX: dict[str, int] = {}
INDEX_TO_CODON: dict[int, str] = {}

_idx = 0
for b1 in _BASES:
    for b2 in _BASES:
        for b3 in _BASES:
            codon = b1 + b2 + b3
            CODON_TO_INDEX[codon] = _idx
            INDEX_TO_CODON[_idx] = codon
            _idx += 1


def codon_to_index(codon: str) -> int:
    """Return the integer index (0-63) for a 3-letter codon string."""
    return CODON_TO_INDEX[codon.upper()]


# ---------------------------------------------------------------------------
# 2. Standard genetic code
# ---------------------------------------------------------------------------
_GENETIC_CODE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
}

AMINO_ACID: np.ndarray = np.array(
    [_GENETIC_CODE[INDEX_TO_CODON[i]] for i in range(64)], dtype="U1"
)


def is_stop_codon(index: int) -> bool:
    """Return True if the codon at *index* is a stop codon."""
    return AMINO_ACID[index] == "*"


# ---------------------------------------------------------------------------
# 3. N_SITES / S_SITES — fractional site counts per position
# ---------------------------------------------------------------------------
N_SITES: np.ndarray = np.zeros((64, 3), dtype=np.float64)
S_SITES: np.ndarray = np.zeros((64, 3), dtype=np.float64)

_base_index = {b: i for i, b in enumerate(_BASES)}

for _i in range(64):
    _codon = INDEX_TO_CODON[_i]
    _aa = AMINO_ACID[_i]
    if _aa == "*":
        # Stop codons get 0 at all positions
        continue
    for _pos in range(3):
        _n_count = 0
        _s_count = 0
        for _alt in _BASES:
            if _alt == _codon[_pos]:
                continue
            _mutant = _codon[:_pos] + _alt + _codon[_pos + 1 :]
            _mut_aa = _GENETIC_CODE[_mutant]
            if _mut_aa == _aa:
                _s_count += 1
            else:
                _n_count += 1
        N_SITES[_i, _pos] = _n_count / 3.0
        S_SITES[_i, _pos] = _s_count / 3.0

# ---------------------------------------------------------------------------
# 4. N_DIFFS / S_DIFFS — pairwise differences
# ---------------------------------------------------------------------------
N_DIFFS: np.ndarray = np.zeros((64, 64), dtype=np.float64)
S_DIFFS: np.ndarray = np.zeros((64, 64), dtype=np.float64)


def _classify_pathway(
    codon_from: str, codon_to: str, positions: tuple[int, ...]
) -> tuple[float, float] | None:
    """Walk a single mutational pathway (ordered positions). Return (nd, sd)
    or None if the pathway passes through a stop codon."""
    nd = 0.0
    sd = 0.0
    current = list(codon_from)
    for pos in positions:
        prev_aa = _GENETIC_CODE["".join(current)]
        current[pos] = codon_to[pos]
        new_codon = "".join(current)
        new_aa = _GENETIC_CODE[new_codon]
        # If intermediate is a stop codon, pathway is invalid
        if new_aa == "*" and new_codon != codon_to:
            return None
        if new_aa != prev_aa:
            nd += 1.0
        else:
            sd += 1.0
    return nd, sd


for _i in range(64):
    for _j in range(_i + 1, 64):
        _ci = INDEX_TO_CODON[_i]
        _cj = INDEX_TO_CODON[_j]

        # Find differing positions
        _diff_pos = [p for p in range(3) if _ci[p] != _cj[p]]
        _ndiff = len(_diff_pos)

        if _ndiff == 0:
            continue

        if _ndiff == 1:
            # Single step — classify directly
            _aa_i = AMINO_ACID[_i]
            _aa_j = AMINO_ACID[_j]
            if _aa_i != _aa_j:
                N_DIFFS[_i, _j] = 1.0
            else:
                S_DIFFS[_i, _j] = 1.0
        else:
            # 2- or 3-step: average over all shortest pathways
            _total_nd = 0.0
            _total_sd = 0.0
            _valid = 0
            for _perm in permutations(_diff_pos):
                _result = _classify_pathway(_ci, _cj, _perm)
                if _result is not None:
                    _total_nd += _result[0]
                    _total_sd += _result[1]
                    _valid += 1
            if _valid > 0:
                N_DIFFS[_i, _j] = _total_nd / _valid
                S_DIFFS[_i, _j] = _total_sd / _valid

        # Symmetric
        N_DIFFS[_j, _i] = N_DIFFS[_i, _j]
        S_DIFFS[_j, _i] = S_DIFFS[_i, _j]

# Make arrays read-only
N_SITES.flags.writeable = False
S_SITES.flags.writeable = False
N_DIFFS.flags.writeable = False
S_DIFFS.flags.writeable = False
AMINO_ACID.flags.writeable = False
