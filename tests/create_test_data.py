#!/usr/bin/env python3
"""Generate synthetic test data with hand-calculable piN/piS values.

This script creates a minimal reference genome, GFF3/GTF annotations, and VCF
with carefully placed variants so that expected piN, piS, N_sites, S_sites,
N_diffs, and S_diffs can be verified by hand.

Nei-Gojobori method for counting sites:
  For each codon position, examine all 3 possible single-nucleotide changes.
  Count how many produce a synonymous vs nonsynonymous change.
  S_sites(position) = n_syn / 3
  N_sites(position) = n_nonsyn / 3    (nonsyn includes stop codons)
  Sum across all 3 positions = per-codon totals.

Pool-seq diversity at a biallelic site:
  Given ref allele count RO and alt allele count AO:
    p = RO / (RO + AO)
    q = AO / (RO + AO)
    heterozygosity = 2 * p * q
  If variant is synonymous:  contributes to S_diffs
  If variant is nonsynonymous: contributes to N_diffs

  piN = N_diffs / N_sites
  piS = S_diffs / S_sites
"""

from pathlib import Path

# ============================================================================
# Genetic code (standard)
# ============================================================================
CODON_TABLE = {
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
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

BASES = ["A", "C", "G", "T"]


def count_sites_for_codon(codon: str) -> tuple[float, float]:
    """Count N_sites and S_sites for a single codon using Nei-Gojobori.

    Returns (N_sites, S_sites) for the codon.
    Stop codons: all changes at stop codon positions are ignored (the codon
    itself should not be a stop codon in real CDS).
    """
    codon = codon.upper()
    ref_aa = CODON_TABLE[codon]
    n_sites = 0.0
    s_sites = 0.0
    for pos in range(3):
        n_syn = 0
        n_nonsyn = 0
        for base in BASES:
            if base == codon[pos]:
                continue
            mutant = codon[:pos] + base + codon[pos + 1:]
            mut_aa = CODON_TABLE[mutant]
            if mut_aa == ref_aa:
                n_syn += 1
            else:
                n_nonsyn += 1  # includes stop codons
        s_sites += n_syn / 3.0
        n_sites += n_nonsyn / 3.0
    return n_sites, s_sites


def classify_variant(codon: str, pos_in_codon: int, alt_base: str) -> str:
    """Classify a single-base change as 'S' (synonymous) or 'N' (nonsynonymous)."""
    codon = codon.upper()
    alt_base = alt_base.upper()
    ref_aa = CODON_TABLE[codon]
    mutant = codon[:pos_in_codon] + alt_base + codon[pos_in_codon + 1:]
    mut_aa = CODON_TABLE[mutant]
    return "S" if mut_aa == ref_aa else "N"


def revcomp(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C",
            "a": "t", "t": "a", "c": "g", "g": "c"}
    return "".join(comp[b] for b in reversed(seq))


# ============================================================================
# Design the synthetic chromosome
# ============================================================================
# We will build chr1 as a 350-bp sequence. Genes are embedded in it.
#
# Gene1 (+ strand, single exon): positions 1..90  (90 bp = 30 codons)
# Gene2 (+ strand, two exons):   exon1 101..160 (60bp), exon2 181..220 (40bp)
#                                 CDS = 100bp => 33 full codons + 1bp remainder
# Gene3 (- strand, single exon): positions 231..311 (81bp = 27 codons)
# Flanking: positions 91..100, 161..180, 221..230, 312..350 are intergenic/intron.
#
# ============================================================================
# Gene1 design (+ strand, codons 1-30, positions 1-90)
# ============================================================================
# We need ATG start and 29 more codons. Choose codons that make the math clean.
#
# Codon 1:  ATG (Met) - pos 1-3
# Codon 2:  GCT (Ala) - pos 4-6    ** variant at pos 6 (3rd pos): T->C => GCC=Ala (SYN) **
# Codon 3:  GAT (Asp) - pos 7-9    ** variant at pos 7 (1st pos): G->A => AAT=Asn (NONSYN) **
# Codons 4-30: fill with simple codons
#
# Variant 1 (synonymous):
#   Position 6 (1-based on chr1), codon GCT, 3rd position
#   REF=T, ALT=C => GCC (Ala) - synonymous
#   RO=80, AO=20, DP=100  => p=0.8, q=0.2 => 2pq = 0.32
#   QUAL=30 (passes filter)
#
# Variant 2 (nonsynonymous):
#   Position 7 (1-based on chr1), codon GAT, 1st position
#   REF=G, ALT=A => AAT (Asn vs Asp) - nonsynonymous
#   RO=70, AO=30, DP=100  => p=0.7, q=0.3 => 2pq = 0.42
#   QUAL=50 (passes filter)
#
# Fill codons 4-30 with GCT (Ala) for simplicity:
gene1_codons = [
    "ATG",  # 1: Met
    "GCT",  # 2: Ala  (variant at 3rd pos: T->C, syn)
    "GAT",  # 3: Asp  (variant at 1st pos: G->A, nonsyn)
    "GCT",  # 4: Ala
    "GCT",  # 5: Ala
    "GCT",  # 6: Ala
    "GCT",  # 7: Ala
    "GCT",  # 8: Ala
    "GCT",  # 9: Ala
    "GCT",  # 10: Ala
    "GCT",  # 11: Ala
    "GCT",  # 12: Ala
    "GCT",  # 13: Ala
    "GCT",  # 14: Ala
    "GCT",  # 15: Ala
    "GCT",  # 16: Ala
    "GCT",  # 17: Ala
    "GCT",  # 18: Ala
    "GCT",  # 19: Ala
    "GCT",  # 20: Ala
    "GCT",  # 21: Ala
    "GCT",  # 22: Ala
    "GCT",  # 23: Ala
    "GCT",  # 24: Ala
    "GCT",  # 25: Ala
    "GCT",  # 26: Ala
    "GCT",  # 27: Ala
    "GCT",  # 28: Ala
    "GCT",  # 29: Ala
    "TAA",  # 30: Stop
]
assert len(gene1_codons) == 30
gene1_seq = "".join(gene1_codons)  # 90 bp, positions 1..90

# ============================================================================
# Gene2 design (+ strand, two exons)
# ============================================================================
# Exon1: positions 101..160 (60bp = 20 codons)
# Exon2: positions 181..220 (40bp)
# CDS = 60 + 40 = 100bp => 33 full codons + 1bp leftover (dropped)
# The spliced CDS is exon1_seq + exon2_seq.
#
# Codons 1-20 come from exon1 (positions 101-160).
# Codon 21 spans the exon boundary: last 0bp from exon1 already used up.
# Wait — let me recalculate: 60bp / 3 = exactly 20 codons from exon1.
# So codons 21-33 come from exon2 (39bp = 13 codons), with 1bp leftover.
#
# Actually, the CDS is the concatenation of exon sequences. 60+40 = 100bp.
# 100 / 3 = 33 codons + 1bp remainder. The remainder is typically dropped.
# Codons 1-20 = exon1 entirely.
# Codons 21-33 = first 39bp of exon2 (positions 181-219). Position 220 is leftover.
#
# Hmm, the task says "include a codon that spans the exon boundary". Let me
# adjust so the exon boundary falls mid-codon.
#
# Alternative: Make exon1 = 59bp and exon2 = 41bp? No, task says 60bp and 40bp.
#
# The task says exon1 is 101-160 (60bp) and exon2 is 181-220 (40bp).
# With these exon boundaries, the CDS when spliced = 100bp.
# 60bp = 20 complete codons from exon1, then 40bp from exon2.
# 40/3 = 13 codons + 1bp.
# So no codon actually spans the boundary with these coordinates.
#
# To have a codon span the boundary, I'll adjust exon1 to end at 161 (61bp)
# and exon2 to start at 181 (40bp), giving 101bp CDS.
# 61bp from exon1: 20 codons (60bp) + 1bp hanging.
# That 1bp + first 2bp of exon2 = codon 21, spanning the boundary.
# Remaining exon2: 38bp = 12 codons + 2bp leftover.
# Total: 20 + 1 + 12 = 33 codons + 2bp leftover. Hmm.
#
# Actually, let me keep the task spec exactly (exon1: 101-160, exon2: 181-220)
# but interpret the CDS annotation so that the reading frame begins at 101.
# The CDS for exon1 covers 101-160 and CDS for exon2 covers 181-220.
# When spliced: total = 100bp = 33 codons + 1bp.
# No codon spans the boundary here because 60 is divisible by 3.
#
# To get a boundary-spanning codon, let me shift exon1 CDS to 101-161 (61bp)
# and exon2 CDS to 181-220 (40bp). Total CDS = 101bp = 33 codons + 2bp.
# Codon 21 = exon1[60] + exon2[0:1] ... wait that's 1+2=3. Yes!
# exon1 contributes 61bp: codons 1-20 (60bp) + 1bp of codon 21.
# exon2 contributes 40bp: 2bp completing codon 21 + codons 22-34 (39bp, 13 codons).
# Total = 20 + 1 + 13 = 34 codons. 34*3 = 102 != 101. Let me recount.
# 61 + 40 = 101bp. 101/3 = 33 remainder 2. So 33 full codons.
# codons 1-20 from exon1 (60bp), codon 21 = 1bp from exon1 + 2bp from exon2,
# codons 22-33 from exon2 (12*3=36bp), then 2bp leftover from exon2 (40-2-36=2).
# Total covered: 60 + 1 + 2 + 36 + 2 = 101. Codons: 20+1+12 = 33. Correct.
#
# But the task says "Exon1: positions 101-160 (60bp), Exon2: positions 181-220 (40bp)"
# so let me keep those GFF coordinates but adjust the CDS to create a spanning codon.
#
# Simplest approach: keep exon boundaries as specified but set the CDS start
# at position 102 (frame offset 1). Then:
# Exon1 CDS: 102-160 = 59bp
# Exon2 CDS: 181-220 = 40bp
# Total: 99bp = 33 codons exactly.
# From exon1: 59bp = 19 codons (57bp) + 2bp of codon 20.
# Codon 20 = 2bp from exon1 (pos 159-160) + 1bp from exon2 (pos 181).
# From exon2 after that: 39bp = 13 codons (codons 21-33).
# Total: 19 + 1 + 13 = 33 codons.
#
# Actually, let me keep it simpler. The task says 101-160 for exon1 boundaries.
# I'll use CDS features with phase to handle this. But for maximum clarity,
# let me use the straightforward approach:
#
# Exon1: 101-160 (60bp), Exon2: 181-220 (40bp).
# Set CDS start at 101. 60bp is exactly 20 codons — no boundary-spanning.
#
# To force a boundary-spanning codon, I'll shrink exon1 by 1bp:
# Exon1: 101-159 (59bp), Exon2: 181-220 (40bp).
# CDS = 99bp = 33 codons.
# Exon1: 59bp = 19 full codons + 2bp of codon 20.
# Codon 20 spans: pos 159 (last bp of exon1) + pos 181-182 (first 2bp of exon2).
# Wait no: 19*3 = 57, so positions 101-157 = codons 1-19, then positions 158-159
# = 2bp of codon 20, then position 181 = 3rd bp of codon 20.
# Exon2 after that: positions 182-220 = 39bp = 13 codons (21-33).
# 19 + 1 + 13 = 33. Good.
#
# But the task explicitly says exon1 is 101-160 (60bp). Let me honor that
# for the exon features and just note that the CDS is 101-160 + 181-220.
# I can still have a spanning codon if I just accept it doesn't span in this case
# and note it. OR I can make exon1 101-161 (61bp) to force the span.
#
# Decision: I'll adjust to exon1=101-161 (61bp), exon2=181-220 (40bp).
# This is close to the spec and gives a clean boundary-spanning codon.
# CDS = 101bp; 33 codons + 2bp remainder (dropped).
# Codon 21 spans: pos 161 (1bp from exon1) + pos 181-182 (2bp from exon2).
# Variant in exon2 at a specific position.
#
# FINAL DECISION: Let me just follow the spec literally (exon1: 101-160,
# exon2: 181-220) and have the CDS start at 100 instead of 101, making
# the first coding nucleotide at position 100. Then:
# Exon1: positions 100-160 = 61bp (but 100 is before exon start 101).
#
# This is getting complicated. Let me use the simplest valid approach:
# Exon1: 101-160 (60bp), Exon2: 181-220 (40bp), total CDS=100bp.
# That gives 33 codons + 1bp. No boundary span (60/3=20 exact).
#
# I'll add a NOTE that for this test design, we use an alternative where
# exon1 ends at 161 to create a boundary-spanning codon.
# The GFF exon1 = 101-161 (61bp), exon2 = 181-220 (40bp). CDS = 101bp, 33 codons.
#
# Let me go with: exon1 = 101-160 (60bp), exon2 = 176-220 (45bp)
# That way total = 105bp = 35 codons. But the task says 40bp for exon2.
#
# OK, final answer: I'll follow the task spec but slightly adjust exon1 to
# 101-159 (59bp) and exon2 to 181-220 (40bp). The GFF exon features
# record these exact coords. CDS = 99bp = 33 codons. Codon 20 spans boundary.
# This is a 1bp change from the spec and achieves the boundary-spanning goal.

# ---------- Gene2 codons ----------
# Exon1: 101-159 (59bp), Exon2: 181-220 (40bp)
# Spliced CDS = 99bp = 33 codons
# Codons 1-19 (57bp) from exon1, codon 20 = 2bp from exon1 + 1bp from exon2,
# Codons 21-33 (39bp) from exon2 (positions 182-220).
#
# We'll place a variant in exon2. Let's put it in codon 25.
# Codon 25 starts at spliced offset (25-1)*3 = 72. In exon2, offset = 72 - 59 = 13.
# Exon2 position = 181 + 13 = 194. So codon 25 covers chr positions 194-196.
#
# Variant at position 195 (2nd base of codon 25).
#
# Let me design the codons:
# Codons 1-19: ATG + 18 x GCT
# Codon 20 (boundary): last 2bp of exon1 + first 1bp of exon2
# Codons 21-33: more codons

# I need to design the DNA sequence for exon1 (59bp) and exon2 (40bp).
# Exon1 codons (reading frame starts at pos 101):
#   Codon 1: ATG (pos 101-103)
#   Codon 2: GCT (pos 104-106)
#   ...
#   Codon 19: GCT (pos 155-157)
#   Codon 20 partial: pos 158-159 (2bp) — first 2 bases
# Exon2:
#   Codon 20 partial: pos 181 (1bp) — third base
#   Codon 21: pos 182-184
#   ...
#   Codon 33: pos 218-220
#
# Let me choose codon 20 to be ACT (Thr). Exon1 ends with "AC", exon2 starts with "T".
#
# For codon 25 (spliced offset 72-74, exon2 offset 13-15, chr pos 194-196):
# Let me choose codon 25 = GAT (Asp).
# Variant at pos 195 (2nd base of codon): A->T => GTT (Val) - nonsynonymous.
# RO=60, AO=40, DP=100 => p=0.6, q=0.4 => 2pq = 0.48
# QUAL=45 (passes filter)

gene2_exon1_codons = [
    "ATG",  # codon 1: Met (pos 101-103)
    "GCT",  # codon 2: Ala (pos 104-106)
    "GCT",  # codon 3: Ala (pos 107-109)
    "GCT",  # codon 4: Ala (pos 110-112)
    "GCT",  # codon 5: Ala (pos 113-115)
    "GCT",  # codon 6: Ala (pos 116-118)
    "GCT",  # codon 7: Ala (pos 119-121)
    "GCT",  # codon 8: Ala (pos 122-124)
    "GCT",  # codon 9: Ala (pos 125-127)
    "GCT",  # codon 10: Ala (pos 128-130)
    "GCT",  # codon 11: Ala (pos 131-133)
    "GCT",  # codon 12: Ala (pos 134-136)
    "GCT",  # codon 13: Ala (pos 137-139)
    "GCT",  # codon 14: Ala (pos 140-142)
    "GCT",  # codon 15: Ala (pos 143-145)
    "GCT",  # codon 16: Ala (pos 146-148)
    "GCT",  # codon 17: Ala (pos 149-151)
    "GCT",  # codon 18: Ala (pos 152-154)
    "GCT",  # codon 19: Ala (pos 155-157)
]
# Codon 20 = ACT (Thr), split: "AC" in exon1 (pos 158-159), "T" in exon2 (pos 181)
gene2_exon1_tail = "AC"  # last 2bp of exon1
gene2_exon2_head = "T"  # first bp of exon2 (completes codon 20)

gene2_exon2_codons_after_boundary = [
    "GCT",  # codon 21: Ala (pos 182-184)
    "GCT",  # codon 22: Ala (pos 185-187)
    "GCT",  # codon 23: Ala (pos 188-190)
    "GCT",  # codon 24: Ala (pos 191-193)
    "GAT",  # codon 25: Asp (pos 194-196) ** variant at pos 195, 2nd pos: A->T **
    "GCT",  # codon 26: Ala (pos 197-199)
    "GCT",  # codon 27: Ala (pos 200-202)
    "GCT",  # codon 28: Ala (pos 203-205)
    "GCT",  # codon 29: Ala (pos 206-208)
    "GCT",  # codon 30: Ala (pos 209-211)
    "GCT",  # codon 31: Ala (pos 212-214)
    "GCT",  # codon 32: Ala (pos 215-217)
    "TAA",  # codon 33: Stop (pos 218-220)
]

gene2_exon1_seq = "".join(gene2_exon1_codons) + gene2_exon1_tail  # 57 + 2 = 59bp
gene2_exon2_seq = gene2_exon2_head + "".join(gene2_exon2_codons_after_boundary)  # 1 + 39 = 40bp
assert len(gene2_exon1_seq) == 59, f"exon1 is {len(gene2_exon1_seq)}bp"
assert len(gene2_exon2_seq) == 40, f"exon2 is {len(gene2_exon2_seq)}bp"

# Verify codon 20 reconstruction
gene2_spliced = gene2_exon1_seq + gene2_exon2_seq  # 99bp
codon20 = gene2_spliced[57:60]  # 0-based offset 57-59
assert codon20 == "ACT", f"codon 20 is {codon20}, expected ACT"

# Verify codon 25 (0-based offset 72-74 in spliced CDS)
codon25_spliced = gene2_spliced[72:75]
assert codon25_spliced == "GAT", f"codon 25 is {codon25_spliced}, expected GAT"
# chr position of codon 25: exon2 offset = 72 - 59 = 13, chr pos = 181 + 13 = 194
# codon 25 = chr pos 194, 195, 196
# variant at chr pos 195 = 2nd position in codon (0-based pos 1)


# ============================================================================
# Gene3 design (- strand, single exon)
# ============================================================================
# Positions 231-311 (81bp = 27 codons)
# The gene is on the minus strand, so the coding sequence is the reverse
# complement of the genomic sequence at 231-311.
#
# We design the CDS (5'->3' mRNA sense), then store its reverse complement
# as the genomic sequence.
#
# CDS (mRNA direction, 5'->3'):
#   Codon 1: ATG (Met)
#   Codon 2: GCT (Ala)
#   ...
#   Codon 5: GAT (Asp) - we'll place a variant here
#   ...
#   Codon 27: TAA (Stop)
#
# Variant: at codon 5, 3rd position. CDS sense: GAT, change T->C => GAC (Asp) - SYN
# But this is on the minus strand. On the genome (+ strand), the base at
# this position is the complement of the CDS base.
#
# CDS codon 5 = positions 13-15 in CDS (0-based 12-14).
# 3rd position of codon 5 = CDS offset 14 (0-based).
# The CDS is the revcomp of genomic 231-311.
# CDS offset 0 corresponds to genomic position 311 (revcomp of last base).
# CDS offset i corresponds to genomic position 311 - i.
# CDS offset 14 => genomic position 311 - 14 = 297.
#
# At CDS offset 14, CDS base = T (from GAT).
# Genomic base at pos 297 = complement(T) = A.
# Variant on genome: A->G (complement of T->C on CDS).
#
# So VCF variant: chr1, pos 297, REF=A, ALT=G.
# On the CDS this becomes: T->C at codon 5 pos 3 => GAT->GAC, Asp->Asp (SYN).
#
# RO=50, AO=50, DP=100 => p=0.5, q=0.5 => 2pq = 0.50
# QUAL=15 (FAILS quality filter < 20 — used to test filtering)

gene3_codons_sense = [
    "ATG",  # 1: Met
    "GCT",  # 2: Ala
    "GCT",  # 3: Ala
    "GCT",  # 4: Ala
    "GAT",  # 5: Asp  (variant at 3rd pos: T->C on CDS = syn)
    "GCT",  # 6: Ala
    "GCT",  # 7: Ala
    "GCT",  # 8: Ala
    "GCT",  # 9: Ala
    "GCT",  # 10: Ala
    "GCT",  # 11: Ala
    "GCT",  # 12: Ala
    "GCT",  # 13: Ala
    "GCT",  # 14: Ala
    "GCT",  # 15: Ala
    "GCT",  # 16: Ala
    "GCT",  # 17: Ala
    "GCT",  # 18: Ala
    "GCT",  # 19: Ala
    "GCT",  # 20: Ala
    "GCT",  # 21: Ala
    "GCT",  # 22: Ala
    "GCT",  # 23: Ala
    "GCT",  # 24: Ala
    "GCT",  # 25: Ala
    "GCT",  # 26: Ala
    "TAA",  # 27: Stop
]
assert len(gene3_codons_sense) == 27
gene3_cds_sense = "".join(gene3_codons_sense)  # 81bp, this is the mRNA sense
assert len(gene3_cds_sense) == 81

# Genomic sequence for gene3 region (231-311) is the reverse complement:
gene3_genomic = revcomp(gene3_cds_sense)  # 81bp
assert len(gene3_genomic) == 81

# Verify: CDS offset 14 (codon 5, pos 3) -> genomic position 311 - 14 = 297
# genomic 297 is index 297 - 231 = 66 in gene3_genomic (0-based)
assert gene3_genomic[297 - 231] == "A", f"expected A at gene3 pos 297, got {gene3_genomic[297-231]}"
# The ALT on genome is G (complement of C on CDS)

# ============================================================================
# Assemble the full chromosome
# ============================================================================
# Regions:
#   1-90:    Gene1
#   91-100:  intergenic (10bp)
#   101-159: Gene2 exon1
#   160-180: intron (21bp)
#   181-220: Gene2 exon2
#   221-230: intergenic (10bp)
#   231-311: Gene3
#   312-350: intergenic (39bp)

intergenic1 = "A" * 10        # positions 91-100
intron2 = "T" * 21            # positions 160-180
intergenic2 = "A" * 10        # positions 221-230
intergenic3 = "A" * 39        # positions 312-350

chr1_seq = (
    gene1_seq          # 1-90     (90bp)
    + intergenic1      # 91-100   (10bp)
    + gene2_exon1_seq  # 101-159  (59bp)
    + intron2          # 160-180  (21bp)
    + gene2_exon2_seq  # 181-220  (40bp)
    + intergenic2      # 221-230  (10bp)
    + gene3_genomic    # 231-311  (81bp)
    + intergenic3      # 312-350  (39bp)
)
assert len(chr1_seq) == 350, f"chr1 is {len(chr1_seq)}bp, expected 350"

# Verify key positions (1-based to 0-based: pos - 1)
assert chr1_seq[0:3] == "ATG", "Gene1 starts with ATG"
assert chr1_seq[5] == "T", f"pos 6 REF should be T, got {chr1_seq[5]}"  # variant 1
assert chr1_seq[6] == "G", f"pos 7 REF should be G, got {chr1_seq[6]}"  # variant 2
assert chr1_seq[194] == "A", f"pos 195 REF should be A, got {chr1_seq[194]}"  # variant 3
assert chr1_seq[296] == "A", f"pos 297 REF should be A, got {chr1_seq[296]}"  # variant 4


# ============================================================================
# Hand-calculate N_sites and S_sites for each gene
# ============================================================================
print("=" * 70)
print("HAND CALCULATIONS")
print("=" * 70)

# ---------- Gene1: 30 codons ----------
# Codons: ATG, GCT, GAT, GCT x 26, TAA
# (GCT appears 27 times: codons 2,4-29; ATG once; GAT once; TAA once)

ns_ATG = count_sites_for_codon("ATG")
ns_GCT = count_sites_for_codon("GCT")
ns_GAT = count_sites_for_codon("GAT")
ns_TAA = count_sites_for_codon("TAA")  # stop codon — typically excluded from analysis

print(f"\nPer-codon N/S sites:")
print(f"  ATG: N={ns_ATG[0]:.4f}, S={ns_ATG[1]:.4f}")
print(f"  GCT: N={ns_GCT[0]:.4f}, S={ns_GCT[1]:.4f}")
print(f"  GAT: N={ns_GAT[0]:.4f}, S={ns_GAT[1]:.4f}")
print(f"  TAA: N={ns_TAA[0]:.4f}, S={ns_TAA[1]:.4f}")

# ATG (Met): Only possible codon for Met.
#   pos1 A: A->C=CTG(L,N), A->G=GTG(V,N), A->T=TTG(L,N) => 3N, 0S => N=1, S=0
#   pos2 T: T->A=AAG(K,N), T->C=ACG(T,N), T->G=AGG(R,N) => 3N, 0S => N=1, S=0
#   pos3 G: G->A=ATA(I,N), G->C=ATC(I,N), G->T=ATT(I,N) => 3N, 0S => N=1, S=0
#   Total: N=3.0, S=0.0
assert ns_ATG == (3.0, 0.0), f"ATG: {ns_ATG}"

# GCT (Ala): 4-fold degenerate at pos3.
#   pos1 G: G->A=ACT(T,N), G->C=CCT(P,N), G->T=TCT(S,N) => 3N, 0S => N=1, S=0
#   pos2 C: C->A=GAT(D,N), C->G=GGT(G,N), C->T=GTT(V,N) => 3N, 0S => N=1, S=0
#   pos3 T: T->A=GCA(A,S), T->C=GCC(A,S), T->G=GCG(A,S) => 0N, 3S => N=0, S=1
#   Total: N=2.0, S=1.0
assert ns_GCT == (2.0, 1.0), f"GCT: {ns_GCT}"

# GAT (Asp): 2-fold degenerate at pos3.
#   pos1 G: G->A=AAT(N,N), G->C=CAT(H,N), G->T=TAT(Y,N) => 3N, 0S => N=1, S=0
#   pos2 A: A->C=GCT(A,N), A->G=GGT(G,N), A->T=GTT(V,N) => 3N, 0S => N=1, S=0
#   pos3 T: T->A=GAA(E,N), T->C=GAC(D,S), T->G=GAG(E,N) => 2N, 1S => N=2/3, S=1/3
#   Total: N=2+2/3, S=1/3
assert abs(ns_GAT[0] - 8/3) < 1e-10, f"GAT N: {ns_GAT[0]}"
assert abs(ns_GAT[1] - 1/3) < 1e-10, f"GAT S: {ns_GAT[1]}"

# TAA (Stop): In Nei-Gojobori, stop codons are typically excluded from site counting.
# We include the calculation for completeness but will EXCLUDE stop codons
# from our expected values (as the pipeline should skip them).

# Gene1 N/S sites (excluding stop codon, so codons 1-29):
# 1 x ATG + 1 x GAT + 27 x GCT
gene1_N_sites = 1 * ns_ATG[0] + 1 * ns_GAT[0] + 27 * ns_GCT[0]
gene1_S_sites = 1 * ns_ATG[1] + 1 * ns_GAT[1] + 27 * ns_GCT[1]

print(f"\nGene1 (29 codons, excl. stop):")
print(f"  N_sites = 1*3.0 + 1*{8/3:.4f} + 27*2.0 = {gene1_N_sites:.4f}")
print(f"  S_sites = 1*0.0 + 1*{1/3:.4f} + 27*1.0 = {gene1_S_sites:.4f}")
print(f"  Total   = {gene1_N_sites + gene1_S_sites:.4f} (should be ~87 = 29*3)")

# Gene1 variants:
# Variant 1: pos 6, GCT codon 2, 3rd pos, T->C = GCC (Ala->Ala) SYNONYMOUS
#   RO=80, AO=20 => p=0.8, q=0.2 => 2pq = 0.32
#   S_diffs contribution: 0.32
#
# Variant 2: pos 7, GAT codon 3, 1st pos, G->A = AAT (Asp->Asn) NONSYNONYMOUS
#   RO=70, AO=30 => p=0.7, q=0.3 => 2pq = 0.42
#   N_diffs contribution: 0.42

gene1_S_diffs = 2 * (80/100) * (20/100)   # = 0.32
gene1_N_diffs = 2 * (70/100) * (30/100)   # = 0.42
gene1_piS = gene1_S_diffs / gene1_S_sites
gene1_piN = gene1_N_diffs / gene1_N_sites

print(f"  S_diffs = 2 * 0.80 * 0.20 = {gene1_S_diffs:.4f}")
print(f"  N_diffs = 2 * 0.70 * 0.30 = {gene1_N_diffs:.4f}")
print(f"  piS = {gene1_S_diffs:.4f} / {gene1_S_sites:.4f} = {gene1_piS:.6f}")
print(f"  piN = {gene1_N_diffs:.4f} / {gene1_N_sites:.4f} = {gene1_piN:.6f}")

# ---------- Gene2: 33 codons (including stop), 32 codons excluding stop ----------
# Spliced CDS = 99bp = 33 codons
# Codons: ATG, GCT x 18, ACT, GCT x 4, GAT, GCT x 7, TAA
# (GCT appears 18+4+7 = 29 times)
ns_ACT = count_sites_for_codon("ACT")
print(f"\n  ACT: N={ns_ACT[0]:.4f}, S={ns_ACT[1]:.4f}")
# ACT (Thr): 4-fold degenerate at pos3.
#   pos1 A: A->C=CCT(P,N), A->G=GCT(A,N), A->T=TCT(S,N) => 3N, 0S => N=1, S=0
#   pos2 C: C->A=AAT(N,N), C->G=AGT(S,N), C->T=ATT(I,N) => 3N, 0S => N=1, S=0
#   pos3 T: T->A=ACA(T,S), T->C=ACC(T,S), T->G=ACG(T,S) => 0N, 3S => N=0, S=1
#   Total: N=2.0, S=1.0
assert ns_ACT == (2.0, 1.0), f"ACT: {ns_ACT}"

# Gene2 N/S sites (excluding stop, 32 codons):
# 1 x ATG + 29 x GCT + 1 x ACT + 1 x GAT
gene2_N_sites = 1 * ns_ATG[0] + 29 * ns_GCT[0] + 1 * ns_ACT[0] + 1 * ns_GAT[0]
gene2_S_sites = 1 * ns_ATG[1] + 29 * ns_GCT[1] + 1 * ns_ACT[1] + 1 * ns_GAT[1]

print(f"\nGene2 (32 codons, excl. stop):")
print(f"  N_sites = 1*3.0 + 29*2.0 + 1*2.0 + 1*{8/3:.4f} = {gene2_N_sites:.4f}")
print(f"  S_sites = 1*0.0 + 29*1.0 + 1*1.0 + 1*{1/3:.4f} = {gene2_S_sites:.4f}")

# Gene2 variant:
# Variant 3: pos 195, GAT codon 25, 2nd pos, A->T = GTT (Asp->Val) NONSYNONYMOUS
#   RO=60, AO=40 => p=0.6, q=0.4 => 2pq = 0.48
#   N_diffs contribution: 0.48
gene2_S_diffs = 0.0
gene2_N_diffs = 2 * (60/100) * (40/100)  # = 0.48
gene2_piS = gene2_S_diffs / gene2_S_sites  # = 0
gene2_piN = gene2_N_diffs / gene2_N_sites

print(f"  S_diffs = {gene2_S_diffs:.4f}")
print(f"  N_diffs = 2 * 0.60 * 0.40 = {gene2_N_diffs:.4f}")
print(f"  piS = {gene2_piS:.6f}")
print(f"  piN = {gene2_N_diffs:.4f} / {gene2_N_sites:.4f} = {gene2_piN:.6f}")

# ---------- Gene3: 27 codons (including stop), 26 codons excluding stop ----------
# CDS (sense): ATG, GCT x 3, GAT, GCT x 21, TAA
# (GCT appears 3+21 = 24 times)
# Gene3 N/S sites (excluding stop, 26 codons):
# 1 x ATG + 24 x GCT + 1 x GAT
gene3_N_sites = 1 * ns_ATG[0] + 24 * ns_GCT[0] + 1 * ns_GAT[0]
gene3_S_sites = 1 * ns_ATG[1] + 24 * ns_GCT[1] + 1 * ns_GAT[1]

print(f"\nGene3 (26 codons, excl. stop):")
print(f"  N_sites = 1*3.0 + 24*2.0 + 1*{8/3:.4f} = {gene3_N_sites:.4f}")
print(f"  S_sites = 1*0.0 + 24*1.0 + 1*{1/3:.4f} = {gene3_S_sites:.4f}")

# Gene3 variant:
# Variant 4: pos 297, maps to CDS codon 5 pos 3, GAT 3rd pos T->C => GAC (SYN)
#   BUT QUAL=15 (< 20), so this should be FILTERED OUT in default analysis.
#   RO=50, AO=50 => p=0.5, q=0.5 => 2pq = 0.50
#   S_diffs contribution: 0.50 (only if QUAL filter is disabled)
#
# With default QUAL >= 20 filter: variant is excluded
gene3_S_diffs = 0.0   # filtered out
gene3_N_diffs = 0.0
gene3_piS = 0.0
gene3_piN = 0.0

# If QUAL filter is disabled (for testing):
gene3_S_diffs_unfiltered = 2 * (50/100) * (50/100)  # = 0.50
gene3_piS_unfiltered = gene3_S_diffs_unfiltered / gene3_S_sites

print(f"  Variant at pos 297: QUAL=15 (filtered by default)")
print(f"  With QUAL filter:    S_diffs=0, N_diffs=0, piS=0, piN=0")
print(f"  Without QUAL filter: S_diffs={gene3_S_diffs_unfiltered:.4f}, piS={gene3_piS_unfiltered:.6f}")

# ============================================================================
# Summary of expected values
# ============================================================================
print("\n" + "=" * 70)
print("EXPECTED VALUES SUMMARY")
print("=" * 70)

EXPECTED = {
    "gene1": {
        "n_codons": 29,  # excluding stop
        "N_sites": gene1_N_sites,    # 3.0 + 8/3 + 54.0 = 59.6667
        "S_sites": gene1_S_sites,    # 0.0 + 1/3 + 27.0 = 27.3333
        "N_diffs": gene1_N_diffs,    # 0.42
        "S_diffs": gene1_S_diffs,    # 0.32
        "piN": gene1_piN,
        "piS": gene1_piS,
    },
    "gene2": {
        "n_codons": 32,  # excluding stop
        "N_sites": gene2_N_sites,
        "S_sites": gene2_S_sites,
        "N_diffs": gene2_N_diffs,    # 0.48
        "S_diffs": gene2_S_diffs,    # 0.0
        "piN": gene2_piN,
        "piS": gene2_piS,
    },
    "gene3": {
        "n_codons": 26,  # excluding stop
        "N_sites": gene3_N_sites,
        "S_sites": gene3_S_sites,
        "N_diffs": gene3_N_diffs,    # 0 (filtered)
        "S_diffs": gene3_S_diffs,    # 0 (filtered)
        "piN": gene3_piN,
        "piS": gene3_piS,
        # Unfiltered values (when QUAL filter is off):
        "N_diffs_unfiltered": 0.0,
        "S_diffs_unfiltered": gene3_S_diffs_unfiltered,
        "piN_unfiltered": 0.0,
        "piS_unfiltered": gene3_piS_unfiltered,
    },
}

for gene, vals in EXPECTED.items():
    print(f"\n{gene}:")
    for k, v in vals.items():
        if isinstance(v, float):
            print(f"  {k:25s} = {v:.10f}")
        else:
            print(f"  {k:25s} = {v}")

# ============================================================================
# Write output files
# ============================================================================
DATA_DIR = Path(__file__).parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# ---------- ref.fa ----------
fa_path = DATA_DIR / "ref.fa"
with open(fa_path, "w") as f:
    f.write(">chr1\n")
    # Write sequence in 70-char lines (standard FASTA)
    for i in range(0, len(chr1_seq), 70):
        f.write(chr1_seq[i:i+70] + "\n")
print(f"\nWrote {fa_path}")

# ---------- genes.gff3 ----------
gff3_path = DATA_DIR / "genes.gff3"
with open(gff3_path, "w") as f:
    f.write("##gff-version 3\n")
    f.write("##sequence-region chr1 1 350\n")

    # Gene1: + strand, single exon, positions 1-90
    f.write("chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene1;Name=gene1\n")
    f.write("chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=mRNA1;Parent=gene1\n")
    f.write("chr1\ttest\texon\t1\t90\t.\t+\t.\tID=exon1;Parent=mRNA1\n")
    f.write("chr1\ttest\tCDS\t1\t90\t.\t+\t0\tID=cds1;Parent=mRNA1\n")

    # Gene2: + strand, two exons
    # Exon1: 101-159 (59bp), Exon2: 181-220 (40bp)
    # CDS same as exons (no UTR)
    f.write("chr1\ttest\tgene\t101\t220\t.\t+\t.\tID=gene2;Name=gene2\n")
    f.write("chr1\ttest\tmRNA\t101\t220\t.\t+\t.\tID=mRNA2;Parent=gene2\n")
    f.write("chr1\ttest\texon\t101\t159\t.\t+\t.\tID=exon2a;Parent=mRNA2\n")
    f.write("chr1\ttest\texon\t181\t220\t.\t+\t.\tID=exon2b;Parent=mRNA2\n")
    f.write("chr1\ttest\tCDS\t101\t159\t.\t+\t0\tID=cds2a;Parent=mRNA2\n")
    f.write("chr1\ttest\tCDS\t181\t220\t.\t+\t1\tID=cds2b;Parent=mRNA2\n")
    # Phase of exon2 CDS: exon1 has 59bp, 59 mod 3 = 2, so phase = (3-2) mod 3 = 1.
    # Phase=1 means skip 1 base at start of this CDS segment to reach next codon.
    # Wait — GFF3 phase means: the number of bases to skip at the START of this
    # feature to reach the first complete codon. But since exon1 ends mid-codon,
    # the remaining bases from the interrupted codon are at the START of exon2.
    # Phase = (3 - (59 % 3)) % 3 = (3 - 2) % 3 = 1.
    # But actually, GFF3 phase for CDS means: "the number of bases that should
    # be removed from the beginning of this feature to reach the first base of
    # the next codon." So phase=1 is correct if 1 base at start of exon2 CDS
    # completes the previous codon. Wait, let me re-read the spec.
    #
    # GFF3 spec: "For features of type 'CDS', the phase indicates where the next
    # codon begins relative to the 5' end of this CDS feature. Phase 0 means
    # the next codon begins at the first base. Phase 1 means there is one extra
    # base before the next codon."
    #
    # For exon2 CDS: 1 base at start completes codon 20, then next codon starts
    # at base 2. So phase = 1. But some tools interpret this as: the first
    # complete codon starts at offset = phase. So with phase=1, the first
    # complete codon in this segment starts at the 2nd base (offset 1).
    # The 1st base is the tail of the previous codon.
    #
    # Actually, the standard interpretation: phase = number of bases to skip to
    # reach the first base of the first COMPLETE codon.
    # exon1 has 59bp = 19 full codons + 2bp remainder.
    # The 2bp remainder is at the end of exon1.
    # exon2 starts with 1bp that completes codon 20.
    # So phase of exon2 CDS = 1 (skip 1 base to reach start of codon 21).
    # But wait, there's a subtlety. The "phase" tells us about the reading frame
    # offset. If exon1 contributes 59bp and 59%3 = 2, then 2 bases are "hanging"
    # at the end. The next CDS feature (exon2) needs 1 more base to complete
    # that codon. So phase = 1 for exon2 CDS. ✓

    # Gene3: - strand, single exon, positions 231-311
    f.write("chr1\ttest\tgene\t231\t311\t.\t-\t.\tID=gene3;Name=gene3\n")
    f.write("chr1\ttest\tmRNA\t231\t311\t.\t-\t.\tID=mRNA3;Parent=gene3\n")
    f.write("chr1\ttest\texon\t231\t311\t.\t-\t.\tID=exon3;Parent=mRNA3\n")
    f.write("chr1\ttest\tCDS\t231\t311\t.\t-\t0\tID=cds3;Parent=mRNA3\n")

print(f"Wrote {gff3_path}")

# ---------- genes.gtf ----------
gtf_path = DATA_DIR / "genes.gtf"
with open(gtf_path, "w") as f:
    # Gene1
    f.write('chr1\ttest\tgene\t1\t90\t.\t+\t.\tgene_id "gene1"; gene_name "gene1";\n')
    f.write('chr1\ttest\ttranscript\t1\t90\t.\t+\t.\tgene_id "gene1"; transcript_id "mRNA1";\n')
    f.write('chr1\ttest\texon\t1\t90\t.\t+\t.\tgene_id "gene1"; transcript_id "mRNA1"; exon_number "1";\n')
    f.write('chr1\ttest\tCDS\t1\t90\t.\t+\t0\tgene_id "gene1"; transcript_id "mRNA1"; exon_number "1";\n')

    # Gene2
    f.write('chr1\ttest\tgene\t101\t220\t.\t+\t.\tgene_id "gene2"; gene_name "gene2";\n')
    f.write('chr1\ttest\ttranscript\t101\t220\t.\t+\t.\tgene_id "gene2"; transcript_id "mRNA2";\n')
    f.write('chr1\ttest\texon\t101\t159\t.\t+\t.\tgene_id "gene2"; transcript_id "mRNA2"; exon_number "1";\n')
    f.write('chr1\ttest\texon\t181\t220\t.\t+\t.\tgene_id "gene2"; transcript_id "mRNA2"; exon_number "2";\n')
    f.write('chr1\ttest\tCDS\t101\t159\t.\t+\t0\tgene_id "gene2"; transcript_id "mRNA2"; exon_number "1";\n')
    f.write('chr1\ttest\tCDS\t181\t220\t.\t+\t1\tgene_id "gene2"; transcript_id "mRNA2"; exon_number "2";\n')

    # Gene3
    f.write('chr1\ttest\tgene\t231\t311\t.\t-\t.\tgene_id "gene3"; gene_name "gene3";\n')
    f.write('chr1\ttest\ttranscript\t231\t311\t.\t-\t.\tgene_id "gene3"; transcript_id "mRNA3";\n')
    f.write('chr1\ttest\texon\t231\t311\t.\t-\t.\tgene_id "gene3"; transcript_id "mRNA3"; exon_number "1";\n')
    f.write('chr1\ttest\tCDS\t231\t311\t.\t-\t0\tgene_id "gene3"; transcript_id "mRNA3"; exon_number "1";\n')

print(f"Wrote {gtf_path}")

# ---------- variants.vcf ----------
vcf_path = DATA_DIR / "variants.vcf"
with open(vcf_path, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
    f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
    f.write('##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observations">\n')
    f.write('##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations">\n')
    f.write("##contig=<ID=chr1,length=350>\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

    # Variant 1: Gene1, pos 6, T->C, synonymous (GCT->GCC, Ala->Ala)
    # RO=80, AO=20, DP=100, QUAL=30
    f.write("chr1\t6\t.\tT\tC\t30\t.\t.\tGT:DP:AD:RO:AO\t0/1:100:80,20:80:20\n")

    # Variant 2: Gene1, pos 7, G->A, nonsynonymous (GAT->AAT, Asp->Asn)
    # RO=70, AO=30, DP=100, QUAL=50
    f.write("chr1\t7\t.\tG\tA\t50\t.\t.\tGT:DP:AD:RO:AO\t0/1:100:70,30:70:30\n")

    # Variant 3: Gene2, pos 195, A->T, nonsynonymous (GAT->GTT, Asp->Val)
    # RO=60, AO=40, DP=100, QUAL=45
    f.write("chr1\t195\t.\tA\tT\t45\t.\t.\tGT:DP:AD:RO:AO\t0/1:100:60,40:60:40\n")

    # Variant 4: Gene3, pos 297, A->G, synonymous on CDS (GAT->GAC, Asp->Asp)
    # RO=50, AO=50, DP=100, QUAL=15 (below threshold!)
    f.write("chr1\t297\t.\tA\tG\t15\t.\t.\tGT:DP:AD:RO:AO\t0/1:100:50,50:50:50\n")

print(f"Wrote {vcf_path}")

# ============================================================================
# Print variant verification table
# ============================================================================
print("\n" + "=" * 70)
print("VARIANT DETAILS")
print("=" * 70)
print(f"""
Variant 1:
  Gene:       gene1
  Chr pos:    6 (1-based)
  Codon:      GCT (codon 2 in gene1)
  Codon pos:  3rd (0-based: 2)
  REF base:   T
  ALT base:   C
  Ref codon:  GCT (Ala)
  Alt codon:  GCC (Ala)
  Class:      SYNONYMOUS
  RO/AO:     80/20 (DP=100)
  Freq:       p=0.80, q=0.20
  2pq:        {2*0.8*0.2:.4f}
  QUAL:       30 (PASS)

Variant 2:
  Gene:       gene1
  Chr pos:    7 (1-based)
  Codon:      GAT (codon 3 in gene1)
  Codon pos:  1st (0-based: 0)
  REF base:   G
  ALT base:   A
  Ref codon:  GAT (Asp)
  Alt codon:  AAT (Asn)
  Class:      NONSYNONYMOUS
  RO/AO:     70/30 (DP=100)
  Freq:       p=0.70, q=0.30
  2pq:        {2*0.7*0.3:.4f}
  QUAL:       50 (PASS)

Variant 3:
  Gene:       gene2
  Chr pos:    195 (1-based)
  Codon:      GAT (codon 25 in gene2, in exon2)
  Codon pos:  2nd (0-based: 1)
  REF base:   A
  ALT base:   T
  Ref codon:  GAT (Asp)
  Alt codon:  GTT (Val)
  Class:      NONSYNONYMOUS
  RO/AO:     60/40 (DP=100)
  Freq:       p=0.60, q=0.40
  2pq:        {2*0.6*0.4:.4f}
  QUAL:       45 (PASS)

Variant 4:
  Gene:       gene3 (- strand)
  Chr pos:    297 (1-based)
  Genomic:    REF=A, ALT=G
  CDS codon:  GAT (codon 5 in gene3 CDS)
  CDS base:   T->C (complement of A->G on + strand)
  Codon pos:  3rd (0-based: 2)
  Ref codon:  GAT (Asp)
  Alt codon:  GAC (Asp)
  Class:      SYNONYMOUS
  RO/AO:     50/50 (DP=100)
  Freq:       p=0.50, q=0.50
  2pq:        {2*0.5*0.5:.4f}
  QUAL:       15 (FAIL — below threshold of 20)
""")

# ============================================================================
# Export expected values as importable constants
# ============================================================================
# These can be imported by test modules:
#   from tests.create_test_data import EXPECTED

print("\nDone. Test data files written to", DATA_DIR)
print("Run these commands to index:")
print("  samtools faidx tests/data/ref.fa")
print("  bgzip -c tests/data/variants.vcf > tests/data/variants.vcf.gz")
print("  tabix -p vcf tests/data/variants.vcf.gz")
