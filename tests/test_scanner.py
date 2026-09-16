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

    def test_reverse_complement_supports_degenerate_bases(self):
        self.assertEqual(reverse_complement("ATGRYN"), "NRYCAT")

    def test_rejects_invalid_bases(self):
        with self.assertRaisesRegex(ValueError, "unsupported bases"):
            scan_sequence("ABCZ", guide_length=2, pam="NGG")


if __name__ == "__main__":
    unittest.main()
