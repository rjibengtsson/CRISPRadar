from __future__ import annotations

from dataclasses import dataclass

IUPAC_BASES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"C", "G"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

COMPLEMENTS = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
)


@dataclass(frozen=True)
class CandidateGuide:
    guide: str
    pam: str
    guide_start: int
    guide_end: int
    pam_start: int
    pam_end: int
    strand: str


def reverse_complement(sequence: str) -> str:
    normalized = sequence.upper()
    _validate_bases(normalized, "sequence")
    return normalized.translate(COMPLEMENTS)[::-1]


def scan_sequence(sequence: str, guide_length: int, pam: str) -> list[CandidateGuide]:
    if not isinstance(guide_length, int) or isinstance(guide_length, bool):
        raise ValueError("guide_length must be an integer")

    if guide_length <= 0:
        raise ValueError("guide_length must be a positive integer")

    normalized_sequence = sequence.upper()
    normalized_pam = pam.upper()

    if not normalized_pam:
        raise ValueError("pam must not be empty")

    _validate_bases(normalized_sequence, "sequence")
    _validate_bases(normalized_pam, "pam")

    matches: list[CandidateGuide] = []
    matches.extend(
        _scan_strand(normalized_sequence, guide_length, normalized_pam, strand="+")
    )

    reverse_sequence = reverse_complement(normalized_sequence)
    sequence_length = len(normalized_sequence)
    for candidate in _scan_strand(reverse_sequence, guide_length, normalized_pam, strand="-"):
        reverse_guide_start = sequence_length - candidate.guide_end
        reverse_guide_end = sequence_length - candidate.guide_start
        reverse_pam_start = sequence_length - candidate.pam_end
        reverse_pam_end = sequence_length - candidate.pam_start
        matches.append(
            CandidateGuide(
                guide=candidate.guide,
                pam=candidate.pam,
                guide_start=reverse_guide_start,
                guide_end=reverse_guide_end,
                pam_start=reverse_pam_start,
                pam_end=reverse_pam_end,
                strand="-",
            )
        )

    return sorted(matches, key=lambda match: (match.guide_start, match.pam_start, match.strand))


def _scan_strand(
    sequence: str,
    guide_length: int,
    pam: str,
    strand: str,
) -> list[CandidateGuide]:
    pam_length = len(pam)
    window_length = guide_length + pam_length
    matches: list[CandidateGuide] = []

    for start in range(0, len(sequence) - window_length + 1):
        guide = sequence[start : start + guide_length]
        pam_sequence = sequence[start + guide_length : start + window_length]
        if _pam_matches(pam_sequence, pam):
            matches.append(
                CandidateGuide(
                    guide=guide,
                    pam=pam_sequence,
                    guide_start=start,
                    guide_end=start + guide_length,
                    pam_start=start + guide_length,
                    pam_end=start + window_length,
                    strand=strand,
                )
            )

    return matches


def _pam_matches(sequence_fragment: str, pam: str) -> bool:
    return all(
        sequence_base in IUPAC_BASES[pam_base]
        for sequence_base, pam_base in zip(sequence_fragment, pam)
    )


def _validate_bases(sequence: str, field_name: str) -> None:
    invalid_bases = sorted({base for base in sequence if base not in IUPAC_BASES})
    if invalid_bases:
        raise ValueError(
            f"{field_name} contains unsupported bases: {', '.join(invalid_bases)}"
        )
