"""Derive tRNA isotype from its anticodon via the genetic code.

A tRNA's amino-acid identity is, to first order, determined by its anticodon: the
anticodon base-pairs with the mRNA codon, and the codon specifies the amino acid.
So ``isotype = translate(reverse_complement(anticodon))``. This is far cheaper
than a covariance-model isotype scan and agrees with reference tRNAscan-SE for
the standard cases.

Special cases the anticodon alone does not resolve (the reference uses dedicated
models for these) are left to the optional CM-based classifier:

* **Initiator vs elongator Met** — both read ``CAT``; only the elongator is
  reported here as ``Met`` (initiator is ``iMet`` in the reference).
* **Selenocysteine** — ``TCA`` anticodon decodes a recoded ``UGA`` stop; reported
  as ``SeC`` by special case below.
* **Suppressor / pseudo / undetermined** — not inferable from the anticodon.
"""

from __future__ import annotations

from trnascan_py.fasta import reverse_complement

# Standard genetic code: DNA codon (sense strand, 5'->3') -> 3-letter amino acid.
# Stop codons map to "Stop".
_CODON_TABLE: dict[str, str] = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "Stop", "TAG": "Stop",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "Stop", "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}


# Anticodons whose isotype the genetic code alone cannot resolve, so a CM scan
# is needed: CAT decodes Met but also marks bacterial initiator fMet / Ile2
# (lysidine-modified) and eukaryotic initiator iMet; TCA is the recoded SeC/Sup.
AMBIGUOUS_ANTICODONS = frozenset({"CAT", "TCA"})


def isotype_from_anticodon(anticodon: str) -> str:
    """Return the 3-letter isotype for an anticodon (DNA letters, 5'->3').

    Returns ``"Undet"`` if the anticodon is malformed. ``TCA`` is reported as
    ``SeC`` (selenocysteine) rather than ``Sup`` since that is the biological
    isotype; other stop-decoding anticodons are reported as ``Sup`` (suppressor).
    """
    ac = anticodon.strip().upper().replace("U", "T")
    if len(ac) != 3 or any(b not in "ACGT" for b in ac):
        return "Undet"
    if ac == "TCA":
        return "SeC"
    codon = reverse_complement(ac)
    aa = _CODON_TABLE.get(codon, "Undet")
    if aa == "Stop":
        return "Sup"
    return aa


__all__ = ["isotype_from_anticodon", "AMBIGUOUS_ANTICODONS"]
