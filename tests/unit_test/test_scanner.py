import unittest

from crispradar import CandidateGuide, reverse_complement, scan_sequence


class ScanSequenceTests(unittest.TestCase):
    def test_returns_all_forward_strand_matches_for_requested_guide_length(self):
        matches = scan_sequence("AAAAAGGTTTTAGG", guide_length=4, pam="AGG")

        self.assertEqual(
            [match for match in matches if match.strand == "+"],
            [
                CandidateGuide(
                    guide="AAAA",
                    pam="AGG",
                    guide_start=0,
                    guide_end=4,
                    pam_start=4,
                    pam_end=7,
                    strand="+",
                ),
                CandidateGuide(
                    guide="TTTT",
                    pam="AGG",
                    guide_start=7,
                    guide_end=11,
                    pam_start=11,
                    pam_end=14,
                    strand="+",
                ),
            ],
        )

    def test_matches_degenerate_pam_bases(self):
        matches = scan_sequence("AAAATGGCCCCAAA", guide_length=4, pam="TGN")

        self.assertEqual(
            [(match.guide, match.pam) for match in matches if match.strand == "+"],
            [("AAAA", "TGG")],
        )

    def test_ambiguous_sequence_bases_do_not_act_as_pam_wildcards(self):
        matches = scan_sequence("AAAANGG", guide_length=4, pam="AGG")

        self.assertEqual(matches, [])

    def test_scans_reverse_complement_strand(self):
        matches = scan_sequence("CCTAAAA", guide_length=4, pam="AGG")

        self.assertEqual(
            matches,
            [
                CandidateGuide(
                    guide="TTTT",
                    pam="AGG",
                    guide_start=3,
                    guide_end=7,
                    pam_start=0,
                    pam_end=3,
                    strand="-",
                )
            ],
        )

    def test_scans_reverse_complement_strand_with_degenerate_pam(self):
        matches = scan_sequence("CCAAAAA", guide_length=4, pam="TGN")

        self.assertEqual(
            matches,
            [
                CandidateGuide(
                    guide="TTTT",
                    pam="TGG",
                    guide_start=3,
                    guide_end=7,
                    pam_start=0,
                    pam_end=3,
                    strand="-",
                )
            ],
        )

    def test_returns_no_matches_when_sequence_is_too_short(self):
        self.assertEqual(scan_sequence("AAAA", guide_length=4, pam="NGG"), [])

    def test_reverse_complement_supports_degenerate_bases(self):
        self.assertEqual(reverse_complement("ATGRYN"), "NRYCAT")

    def test_reverse_complement_rejects_non_string_sequences(self):
        with self.assertRaisesRegex(ValueError, "sequence must be a string"):
            reverse_complement(None)

    def test_rejects_invalid_bases(self):
        with self.assertRaisesRegex(ValueError, "unsupported bases"):
            scan_sequence("ABCZ", guide_length=2, pam="NGG")

    def test_rejects_non_integer_guide_lengths(self):
        for invalid_guide_length in (4.5, True):
            with self.subTest(invalid_guide_length=invalid_guide_length):
                with self.assertRaisesRegex(ValueError, "guide_length must be an integer"):
                    scan_sequence("AAAATGG", guide_length=invalid_guide_length, pam="NGG")

    def test_rejects_non_positive_guide_lengths(self):
        for invalid_guide_length in (0, -1):
            with self.subTest(invalid_guide_length=invalid_guide_length):
                with self.assertRaisesRegex(ValueError, "guide_length must be a positive integer"):
                    scan_sequence("AAAATGG", guide_length=invalid_guide_length, pam="NGG")

    def test_rejects_invalid_pam_bases(self):
        with self.assertRaisesRegex(ValueError, "unsupported bases"):
            scan_sequence("AAAATGG", guide_length=4, pam="NGZ")

    def test_rejects_empty_pam(self):
        with self.assertRaisesRegex(ValueError, "pam must not be empty"):
            scan_sequence("AAAATGG", guide_length=4, pam="")

    def test_rejects_non_string_sequence_or_pam(self):
        with self.assertRaisesRegex(ValueError, "sequence must be a string"):
            scan_sequence(None, guide_length=4, pam="NGG")

        with self.assertRaisesRegex(ValueError, "pam must be a string"):
            scan_sequence("AAAATGG", guide_length=4, pam=None)


if __name__ == "__main__":
    unittest.main()
