from __future__ import annotations

from dataclasses import dataclass, replace

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
    target_seq: str
    pam: str
    guide_start: int
    guide_end: int
    pam_start: int
    pam_end: int
    strand: str


def reverse_complement(sequence: str) -> str:
    if not isinstance(sequence, str):
        raise ValueError("sequence must be a string")

    normalized = sequence.upper()
    _validate_bases(normalized, "sequence")
    return _reverse_complement_validated(normalized)


def scan_sequence(
        sequence: str, 
        guide_length: int, 
        pam: str,
        pam_position: str = "3prime") -> list[CandidateGuide]:

    if pam_position not in {"5prime", "3prime"}:
        raise ValueError("pam_position must be '5prime' or '3prime'")

    if not isinstance(sequence, str):
        raise ValueError("sequence must be a string")

    if not isinstance(pam, str):
        raise ValueError("pam must be a string")

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

    if len(normalized_sequence) < guide_length + len(normalized_pam):
        return []

    matches: list[CandidateGuide] = []
    matches.extend(
        _scan_strand(
            normalized_sequence, 
            guide_length, 
            normalized_pam, 
            strand="+",
            pam_position=pam_position
        )
    )

    reverse_sequence = _reverse_complement_validated(normalized_sequence)
    sequence_length = len(normalized_sequence)
    for candidate in _scan_strand(reverse_sequence, 
                                  guide_length, normalized_pam, 
                                  strand="-",
                                  pam_position=pam_position):
        reverse_guide_start = sequence_length - candidate.guide_end
        reverse_guide_end = sequence_length - candidate.guide_start
        reverse_pam_start = sequence_length - candidate.pam_end
        reverse_pam_end = sequence_length - candidate.pam_start
        matches.append(
            CandidateGuide(
                target_seq=candidate.target_seq,
                pam=candidate.pam,
                guide_start=reverse_guide_start,
                guide_end=reverse_guide_end,
                pam_start=reverse_pam_start,
                pam_end=reverse_pam_end,
                strand="-",
            )
        )

    one_based_matches = [_to_one_based(match) for match in matches]
    return sorted(
        one_based_matches,
        key=lambda match: (match.guide_start, match.pam_start, match.strand),
    )


def _to_one_based(candidate: CandidateGuide) -> CandidateGuide:
    # Scanning uses 0-based half-open coordinates; report 1-based inclusive.
    return replace(
        candidate,
        guide_start=candidate.guide_start + 1,
        pam_start=candidate.pam_start + 1,
    )


def _scan_strand(
    sequence: str,
    guide_length: int,
    pam: str,
    strand: str,
    pam_position: str = "3prime"
) -> list[CandidateGuide]:
    pam_length = len(pam)
    window_length = guide_length + pam_length
    pam_first = pam_position == "5prime"
    matches: list[CandidateGuide] = []


    for start in range(0, len(sequence) - window_length + 1):
        # Only the offsets within the window change; the window itself is the same.
        if pam_first:
            pam_start = start
            guide_start = start + pam_length
        else:
            guide_start = start
            pam_start = start + guide_length

        guide_end = guide_start + guide_length
        pam_end = pam_start + pam_length

        pam_sequence = sequence[pam_start:pam_end]
        if _pam_matches(pam_sequence, pam):
            matches.append(
                CandidateGuide(
                    target_seq=sequence[guide_start:guide_end],
                    pam=pam_sequence,
                    guide_start=guide_start,
                    guide_end=guide_end,
                    pam_start=pam_start,
                    pam_end=pam_end,
                    strand=strand,
                )
            )

    return matches


def _pam_matches(sequence_fragment: str, pam: str) -> bool:
    return all(
        sequence_base in IUPAC_BASES[pam_base]
        for sequence_base, pam_base in zip(sequence_fragment, pam)
    )


def _reverse_complement_validated(sequence: str) -> str:
    return sequence.translate(COMPLEMENTS)[::-1]


def _validate_bases(sequence: str, field_name: str) -> None:
    invalid_bases = sorted({base for base in sequence if base not in IUPAC_BASES})
    if invalid_bases:
        raise ValueError(
            f"{field_name} contains unsupported bases: {', '.join(invalid_bases)}"
        )
