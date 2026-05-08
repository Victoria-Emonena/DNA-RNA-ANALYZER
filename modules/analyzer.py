"""
modules/analyzer.py — Core DNA/RNA analysis logic
"""
import re

# ── Codon table (DNA codons → amino acid) ──────────────────────────────────
CODON_TABLE = {
    # Phenylalanine
    "UUU": ("Phenylalanine", "Phe", "F"),
    "UUC": ("Phenylalanine", "Phe", "F"),
    # Leucine
    "UUA": ("Leucine", "Leu", "L"),
    "UUG": ("Leucine", "Leu", "L"),
    "CUU": ("Leucine", "Leu", "L"),
    "CUC": ("Leucine", "Leu", "L"),
    "CUA": ("Leucine", "Leu", "L"),
    "CUG": ("Leucine", "Leu", "L"),
    # Isoleucine
    "AUU": ("Isoleucine", "Ile", "I"),
    "AUC": ("Isoleucine", "Ile", "I"),
    "AUA": ("Isoleucine", "Ile", "I"),
    # Methionine (Start)
    "AUG": ("Methionine (Start)", "Met", "M"),
    # Valine
    "GUU": ("Valine", "Val", "V"),
    "GUC": ("Valine", "Val", "V"),
    "GUA": ("Valine", "Val", "V"),
    "GUG": ("Valine", "Val", "V"),
    # Serine
    "UCU": ("Serine", "Ser", "S"),
    "UCC": ("Serine", "Ser", "S"),
    "UCA": ("Serine", "Ser", "S"),
    "UCG": ("Serine", "Ser", "S"),
    "AGU": ("Serine", "Ser", "S"),
    "AGC": ("Serine", "Ser", "S"),
    # Proline
    "CCU": ("Proline", "Pro", "P"),
    "CCC": ("Proline", "Pro", "P"),
    "CCA": ("Proline", "Pro", "P"),
    "CCG": ("Proline", "Pro", "P"),
    # Threonine
    "ACU": ("Threonine", "Thr", "T"),
    "ACC": ("Threonine", "Thr", "T"),
    "ACA": ("Threonine", "Thr", "T"),
    "ACG": ("Threonine", "Thr", "T"),
    # Alanine
    "GCU": ("Alanine", "Ala", "A"),
    "GCC": ("Alanine", "Ala", "A"),
    "GCA": ("Alanine", "Ala", "A"),
    "GCG": ("Alanine", "Ala", "A"),
    # Tyrosine
    "UAU": ("Tyrosine", "Tyr", "Y"),
    "UAC": ("Tyrosine", "Tyr", "Y"),
    # Stop codons
    "UAA": ("Stop (Ochre)", "Stop", "*"),
    "UAG": ("Stop (Amber)", "Stop", "*"),
    "UGA": ("Stop (Opal)", "Stop", "*"),
    # Histidine
    "CAU": ("Histidine", "His", "H"),
    "CAC": ("Histidine", "His", "H"),
    # Glutamine
    "CAA": ("Glutamine", "Gln", "Q"),
    "CAG": ("Glutamine", "Gln", "Q"),
    # Asparagine
    "AAU": ("Asparagine", "Asn", "N"),
    "AAC": ("Asparagine", "Asn", "N"),
    # Lysine
    "AAA": ("Lysine", "Lys", "K"),
    "AAG": ("Lysine", "Lys", "K"),
    # Aspartic acid
    "GAU": ("Aspartic acid", "Asp", "D"),
    "GAC": ("Aspartic acid", "Asp", "D"),
    # Glutamic acid
    "GAA": ("Glutamic acid", "Glu", "E"),
    "GAG": ("Glutamic acid", "Glu", "E"),
    # Cysteine
    "UGU": ("Cysteine", "Cys", "C"),
    "UGC": ("Cysteine", "Cys", "C"),
    # Tryptophan
    "UGG": ("Tryptophan", "Trp", "W"),
    # Arginine
    "CGU": ("Arginine", "Arg", "R"),
    "CGC": ("Arginine", "Arg", "R"),
    "CGA": ("Arginine", "Arg", "R"),
    "CGG": ("Arginine", "Arg", "R"),
    "AGA": ("Arginine", "Arg", "R"),
    "AGG": ("Arginine", "Arg", "R"),
    # Glycine
    "GGU": ("Glycine", "Gly", "G"),
    "GGC": ("Glycine", "Gly", "G"),
    "GGA": ("Glycine", "Gly", "G"),
    "GGG": ("Glycine", "Gly", "G"),
}

DNA_COMPLEMENT = str.maketrans("ATCGatcg", "TAGCtagc")
RNA_COMPLEMENT = str.maketrans("AUCGaucg", "UAGCuagc")


def clean_sequence(raw: str) -> str:
    """Strip whitespace, FASTA headers, and newlines; uppercase."""
    lines = raw.strip().splitlines()
    # Remove FASTA header lines starting with >
    seq_lines = [l for l in lines if not l.startswith(">")]
    return "".join(seq_lines).upper().replace(" ", "").replace("\t", "")


def detect_sequence_type(seq: str) -> tuple[str, str]:
    """
    Returns (type, explanation).
    type is one of: 'DNA', 'RNA', 'INVALID'
    """
    dna_bases = set("ATCG")
    rna_bases = set("AUCG")

    bases = set(seq)

    if not bases:
        return "INVALID", "Sequence is empty."

    invalid = bases - (dna_bases | rna_bases | {"N"})
    if invalid:
        return "INVALID", f"Sequence contains invalid characters: {', '.join(sorted(invalid))}"

    has_t = "T" in bases
    has_u = "U" in bases

    if has_t and has_u:
        return "INVALID", "Sequence contains both T (thymine) and U (uracil), which is not valid for a single strand."
    if has_u:
        return "RNA", "Sequence contains uracil (U) — this is an **RNA** strand."
    return "DNA", "Sequence contains thymine (T) and no uracil (U) — this is a **DNA** strand."


def get_complement_dna(seq: str) -> str:
    return seq.translate(DNA_COMPLEMENT)


def dna_to_mrna(seq: str, strand_type: str) -> tuple[str, str]:
    """
    Convert a DNA strand to mRNA.
    - Template strand: complement + replace T→U
    - Coding strand: directly replace T→U (same as mRNA)
    Returns (mRNA sequence, explanation)
    """
    if strand_type == "Template":
        # mRNA is synthesised complementary to the template strand
        complement = get_complement_dna(seq)
        mrna = complement.replace("T", "U")
        explanation = (
            "The **template strand** is used directly by RNA polymerase. "
            "The mRNA is the complement of the template strand, with T replaced by U."
        )
    else:  # Coding
        mrna = seq.replace("T", "U")
        explanation = (
            "The **coding strand** has the same sequence as the mRNA (except T→U). "
            "The mRNA is produced by simply replacing thymine (T) with uracil (U)."
        )
    return mrna, explanation


def translate_mrna(mrna: str) -> list[dict]:
    """
    Translate mRNA into a list of codon dicts.
    Each dict: {codon, amino_acid_full, amino_acid_3, amino_acid_1, is_stop}
    """
    results = []
    for i in range(0, len(mrna) - 2, 3):
        codon = mrna[i:i+3]
        if len(codon) < 3:
            break
        info = CODON_TABLE.get(codon, ("Unknown", "???", "?"))
        is_stop = info[1] == "Stop"
        results.append({
            "codon": codon,
            "amino_acid_full": info[0],
            "amino_acid_3": info[1],
            "amino_acid_1": info[2],
            "is_stop": is_stop,
        })
        if is_stop:
            break
    return results


def rna_to_protein_chain(codons: list[dict]) -> str:
    """Build a one-letter protein string."""
    return "".join(c["amino_acid_1"] for c in codons if not c["is_stop"])


def compute_stats(seq: str, seq_type: str) -> dict:
    """Compute base composition and other stats."""
    length = len(seq)
    counts = {b: seq.count(b) for b in (set(seq) - {"N"})}

    if seq_type == "DNA":
        gc = counts.get("G", 0) + counts.get("C", 0)
    else:
        gc = counts.get("G", 0) + counts.get("C", 0)

    gc_content = (gc / length * 100) if length else 0.0

    return {
        "length": length,
        "counts": counts,
        "gc_content": round(gc_content, 2),
    }


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA format; returns list of (header, sequence) tuples."""
    results = []
    current_header = ""
    current_seq = []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                results.append((current_header, "".join(current_seq).upper()))
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_seq:
        results.append((current_header, "".join(current_seq).upper()))
    return results
