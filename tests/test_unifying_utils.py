#!/usr/bin/env python
# coding: utf-8
"""
Tests for unifying_utils.py

Tests the modular unifying proof generation functions.
"""

import unittest
from unifying_utils import (
    generate_one_case_unifying_proof,
    generate_unifying_proof,
    generate_two_case_unifying_proof,
    generate_three_case_unifying_proof,
    generate_unifying_proof_from_legacy_params,
)


class TestOneCaseUnifyingProof(unittest.TestCase):
    """Test the one-case proof generation."""

    def test_basic_one_case(self):
        """Test basic one-case proof generation."""
        lemma_name = "test_lemma"
        result = generate_one_case_unifying_proof(lemma_name)
        
        self.assertIn("full_domain_soundness_lemma", result)
        self.assertIn(lemma_name, result)
        self.assertIn("SKEEP", result)
        self.assertIn("SKOLETIN*", result)
        self.assertIn("EXPAND \"g_1\"", result)
        self.assertIn("QED", result)


class TestTwoCaseUnifyingProof(unittest.TestCase):
    """Test the two-case proof generation (backward compatibility)."""

    def test_basic_two_case(self):
        """Test basic two-case proof generation."""
        domain_1 = "left_open(0)"
        domain_2 = "right_open(0)"
        trajectory_1 = "LAMBDA (x: real): COND x > 0 -> g(0), ELSE -> g(x) ENDCOND"
        trajectory_2 = "LAMBDA (x: real): COND x >= 0 -> g(x), ELSE -> g(0) ENDCOND"
        full_traj = "LAMBDA (x: real): g(x)"
        piecewise_split_bools = [False]
        domain_splits = ["0"]
        
        result = generate_two_case_unifying_proof(
            domain_1, domain_2, trajectory_1, trajectory_2,
            full_traj, piecewise_split_bools, domain_splits
        )
        
        self.assertIn("full_domain_soundness_lemma", result)
        self.assertIn("left_open", result)
        self.assertIn("right_open", result)
        self.assertIn("SKOLETIN", result)  # Can be "SKOLETIN*" or "SKOLETIN *"
        self.assertIn("PROPAX", result)
        self.assertIn("QED", result)
        # Should have 10 PROPAX calls for 2 cases
        self.assertGreaterEqual(result.count("(PROPAX)"), 10)  # At least 10

    def test_two_case_with_piecewise_splits(self):
        """Test two-case proof with piecewise splits (should still work)."""
        domain_1 = "left_open(0)"
        domain_2 = "right_open(0)"
        trajectory_1 = "LAMBDA (x: real): COND x > 0 -> g(0), ELSE -> g(x) ENDCOND"
        trajectory_2 = "LAMBDA (x: real): COND x >= 0 -> g(x), ELSE -> g(0) ENDCOND"
        full_traj = "LAMBDA (x: real): g(x)"
        piecewise_split_bools = [True]  # Mark as piecewise split
        domain_splits = ["0"]
        
        result = generate_two_case_unifying_proof(
            domain_1, domain_2, trajectory_1, trajectory_2,
            full_traj, piecewise_split_bools, domain_splits
        )
        
        # Should still generate valid proof
        self.assertIn("full_domain_soundness_lemma", result)
        self.assertIn("QED", result)


class TestThreeCaseUnifyingProof(unittest.TestCase):
    """Test the three-case proof generation (backward compatibility)."""

    def test_basic_three_case(self):
        """Test basic three-case proof generation."""
        domain_1 = "left_open(0)"
        domain_2 = "ci(0, 4)"
        domain_3 = "right_open(4)"
        trajectory_1 = "LAMBDA (x: real): COND x > 0 -> g(0), ELSE -> g(x) ENDCOND"
        trajectory_2 = "LAMBDA (x: real): COND x >= 0 -> g(x), ELSE -> g(0) ENDCOND"
        trajectory_3 = "LAMBDA (x: real): COND x >= 4 -> g(x), ELSE -> g(4) ENDCOND"
        full_traj = "LAMBDA (x: real): g(x)"
        piecewise_split_bools = [False, True]  # Second split is piecewise
        domain_splits = ["0", "4"]
        
        result = generate_three_case_unifying_proof(
            domain_1, domain_2, domain_3,
            trajectory_1, trajectory_2, trajectory_3,
            full_traj, piecewise_split_bools, domain_splits
        )
        
        self.assertIn("full_domain_soundness_lemma", result)
        self.assertIn("left_open", result)
        self.assertIn("ci", result)
        self.assertIn("right_open", result)
        self.assertIn("SKOLETIN 1", result)  # Should use SKOLETIN 1 for 3+ cases
        self.assertIn("CASE", result)  # Should use CASE statements
        self.assertIn("QED", result)
        # Should have 16 PROPAX calls at start + 6 at end = 22 total
        # (But may have more due to PROPAX in branches)
        self.assertGreaterEqual(result.count("(PROPAX)"), 22)

    def test_three_case_no_piecewise_splits(self):
        """Test three-case proof without piecewise splits (should skip)."""
        domain_1 = "left_open(0)"
        domain_2 = "ci(0, 4)"
        domain_3 = "right_open(4)"
        trajectory_1 = "LAMBDA (x: real): g(x)"
        trajectory_2 = "LAMBDA (x: real): g(x)"
        trajectory_3 = "LAMBDA (x: real): g(x)"
        full_traj = "LAMBDA (x: real): g(x)"
        piecewise_split_bools = [False, False]  # No piecewise splits
        domain_splits = ["0", "4"]
        
        result = generate_three_case_unifying_proof(
            domain_1, domain_2, domain_3,
            trajectory_1, trajectory_2, trajectory_3,
            full_traj, piecewise_split_bools, domain_splits
        )
        
        # Should return skip message
        self.assertIn("skipping", result.lower())


class TestModularUnifyingProof(unittest.TestCase):
    """Test the modular generate_unifying_proof function."""

    def test_two_case_modular(self):
        """Test modular function with 2 cases."""
        cases = [
            {
                "domain": "left_open(0)",
                "trajectory": "LAMBDA (x: real): COND x > 0 -> g(0), ELSE -> g(x) ENDCOND",
                "domain_type": "left_open",
                "domain_split": "0"
            },
            {
                "domain": "right_open(0)",
                "trajectory": "LAMBDA (x: real): COND x >= 0 -> g(x), ELSE -> g(0) ENDCOND",
                "domain_type": "right_open",
                "domain_split": ""
            }
        ]
        full_traj = "LAMBDA (x: real): g(x)"
        
        result = generate_unifying_proof(cases, full_traj)
        
        self.assertIn("full_domain_soundness_lemma", result)
        self.assertIn("left_open", result)
        self.assertIn("right_open", result)
        self.assertIn("SKOLETIN", result)  # Can be "SKOLETIN*" or "SKOLETIN *"
        self.assertIn("QED", result)

    def test_three_case_modular(self):
        """Test modular function with 3 cases."""
        cases = [
            {
                "domain": "left_open(0)",
                "trajectory": "LAMBDA (x: real): COND x > 0 -> g(0), ELSE -> g(x) ENDCOND",
                "domain_type": "left_open",
                "domain_split": "0"
            },
            {
                "domain": "ci(0, 4)",
                "trajectory": "LAMBDA (x: real): COND x >= 0 -> g(x), ELSE -> g(0) ENDCOND",
                "domain_type": "ci",
                "domain_split": ""
            },
            {
                "domain": "right_open(4)",
                "trajectory": "LAMBDA (x: real): COND x >= 4 -> g(x), ELSE -> g(4) ENDCOND",
                "domain_type": "right_open",
                "domain_split": ""
            }
        ]
        full_traj = "LAMBDA (x: real): g(x)"
        piecewise_splits_map = {2: ["4"]}  # Last domain (index 2) has piecewise split at 4
        
        result = generate_unifying_proof(cases, full_traj, piecewise_splits_map)
        
        self.assertIn("full_domain_soundness_lemma", result)
        self.assertIn("left_open", result)
        self.assertIn("ci", result)
        self.assertIn("right_open", result)
        self.assertIn("SKOLETIN 1", result)
        self.assertIn("CASE", result)
        self.assertIn("x!1=4", result)  # Should have piecewise split case
        self.assertIn("QED", result)

    def test_four_case_modular(self):
        """Test modular function with 4 cases (extending beyond original)."""
        cases = [
            {
                "domain": "left_open(0)",
                "trajectory": "LAMBDA (x: real): g(x)",
                "domain_type": "left_open",
                "domain_split": "0"
            },
            {
                "domain": "ci(0, 2)",
                "trajectory": "LAMBDA (x: real): g(x)",
                "domain_type": "ci",
                "domain_split": ""
            },
            {
                "domain": "ci(2, 4)",
                "trajectory": "LAMBDA (x: real): g(x)",
                "domain_type": "ci",
                "domain_split": ""
            },
            {
                "domain": "right_open(4)",
                "trajectory": "LAMBDA (x: real): g(x)",
                "domain_type": "right_open",
                "domain_split": ""
            }
        ]
        full_traj = "LAMBDA (x: real): g(x)"
        
        result = generate_unifying_proof(cases, full_traj)
        
        self.assertIn("full_domain_soundness_lemma", result)
        self.assertIn("SKOLETIN 1", result)
        self.assertIn("CASE", result)
        self.assertIn("QED", result)
        # Should have PROPAX calls
        self.assertIn("PROPAX", result)

    def test_piecewise_splits_interface(self):
        """Test the piecewise splits mapping interface."""
        cases = [
            {
                "domain": "left_open(0)",
                "trajectory": "LAMBDA (x: real): g(x)",
                "domain_type": "left_open",
                "domain_split": "0"
            },
            {
                "domain": "ci(0, 4)",
                "trajectory": "LAMBDA (x: real): g(x)",
                "domain_type": "ci",
                "domain_split": ""
            },
            {
                "domain": "right_open(4)",
                "trajectory": "LAMBDA (x: real): g(x)",
                "domain_type": "right_open",
                "domain_split": ""
            }
        ]
        full_traj = "LAMBDA (x: real): g(x)"
        
        # Test with piecewise splits mapped to last domain (index 2)
        piecewise_splits_map = {2: ["4", "8"]}  # Domain 2 has splits at 4 and 8
        result = generate_unifying_proof(cases, full_traj, piecewise_splits_map)
        
        self.assertIn("x!1=4", result)  # Should use first piecewise split
        self.assertIn("QED", result)
        
        # Test with piecewise splits mapped to middle domain (index 1)
        # Note: For 3-case proofs, piecewise splits are only used for the last domain
        # So this will skip if no piecewise splits for last domain
        piecewise_splits_map = {1: ["2"]}  # Domain 1 has split at 2, but not last domain
        result = generate_unifying_proof(cases, full_traj, piecewise_splits_map)
        
        # For 3-case without piecewise splits at last domain, should skip
        # (This matches original behavior)
        self.assertIn("skipping", result.lower())

    def test_empty_cases(self):
        """Test with empty cases list."""
        result = generate_unifying_proof([], "LAMBDA (x: real): g(x)")
        self.assertIn("No cases provided", result)

    def test_single_case_modular(self):
        """Test modular function with 1 case (should delegate)."""
        cases = [
            {
                "domain": "left_open(0)",
                "trajectory": "LAMBDA (x: real): g(x)",
                "domain_type": "left_open",
                "domain_split": "0"
            }
        ]
        result = generate_unifying_proof(cases, "LAMBDA (x: real): g(x)")
        self.assertIn("Single case", result)


class TestLegacyParamsWrapper(unittest.TestCase):
    """Test the legacy parameters wrapper function."""

    def test_legacy_params_two_case(self):
        """Test legacy params wrapper with 2 cases."""
        domains = ["left_open(0)", "right_open(0)"]
        trajectories = [
            "LAMBDA (x: real): COND x > 0 -> g(0), ELSE -> g(x) ENDCOND",
            "LAMBDA (x: real): COND x >= 0 -> g(x), ELSE -> g(0) ENDCOND"
        ]
        full_traj = "LAMBDA (x: real): g(x)"
        piecewise_split_bools = [False]
        domain_splits = ["0"]
        
        result = generate_unifying_proof_from_legacy_params(
            domains, trajectories, full_traj,
            piecewise_split_bools, domain_splits
        )
        
        self.assertIn("full_domain_soundness_lemma", result)
        self.assertIn("QED", result)

    def test_legacy_params_three_case(self):
        """Test legacy params wrapper with 3 cases."""
        domains = ["left_open(0)", "ci(0, 4)", "right_open(4)"]
        trajectories = [
            "LAMBDA (x: real): g(x)",
            "LAMBDA (x: real): g(x)",
            "LAMBDA (x: real): g(x)"
        ]
        full_traj = "LAMBDA (x: real): g(x)"
        piecewise_split_bools = [False, True]  # Second split is piecewise
        domain_splits = ["0", "4"]
        
        result = generate_unifying_proof_from_legacy_params(
            domains, trajectories, full_traj,
            piecewise_split_bools, domain_splits
        )
        
        self.assertIn("full_domain_soundness_lemma", result)
        self.assertIn("QED", result)

    def test_legacy_params_four_case(self):
        """Test legacy params wrapper with 4 cases."""
        domains = ["left_open(0)", "ci(0, 2)", "ci(2, 4)", "right_open(4)"]
        trajectories = [
            "LAMBDA (x: real): g(x)",
            "LAMBDA (x: real): g(x)",
            "LAMBDA (x: real): g(x)",
            "LAMBDA (x: real): g(x)"
        ]
        full_traj = "LAMBDA (x: real): g(x)"
        piecewise_split_bools = [False, False, True]
        domain_splits = ["0", "2", "4"]
        
        result = generate_unifying_proof_from_legacy_params(
            domains, trajectories, full_traj,
            piecewise_split_bools, domain_splits
        )
        
        self.assertIn("full_domain_soundness_lemma", result)
        self.assertIn("QED", result)


class TestPiecewiseSplitsInterface(unittest.TestCase):
    """Test the piecewise splits interface explanation."""
    
    def test_piecewise_splits_explanation(self):
        """
        Demonstrate the piecewise splits interface.
        
        The piecewise_splits_map is a dictionary that maps domain indices
        to lists of piecewise boundary points for that domain.
        
        Example:
        - If you have 3 domains (indices 0, 1, 2)
        - And domain 2 (the last one) has piecewise boundaries at x=4 and x=8
        - You would pass: piecewise_splits_map = {2: ["4", "8"]}
        
        This allows the proof generator to handle piecewise trajectories
        that have discontinuities or changes at specific points within a domain.
        """
        cases = [
            {"domain": "left_open(0)", "trajectory": "LAMBDA (x: real): g(x)",
             "domain_type": "left_open", "domain_split": "0"},
            {"domain": "ci(0, 4)", "trajectory": "LAMBDA (x: real): g(x)",
             "domain_type": "ci", "domain_split": ""},
            {"domain": "right_open(4)", "trajectory": "LAMBDA (x: real): g(x)",
             "domain_type": "right_open", "domain_split": ""}
        ]
        
        # Example 1: Piecewise splits at the last domain
        piecewise_splits_map = {2: ["4"]}  # Domain index 2 has split at 4
        result1 = generate_unifying_proof(cases, "LAMBDA (x: real): g(x)", piecewise_splits_map)
        self.assertIn("x!1=4", result1)
        
        # Example 2: Multiple piecewise splits at a domain
        piecewise_splits_map = {2: ["4", "8", "12"]}  # Domain index 2 has splits at 4, 8, 12
        result2 = generate_unifying_proof(cases, "LAMBDA (x: real): g(x)", piecewise_splits_map)
        # Should use first split (4) in the proof
        self.assertIn("x!1=4", result2)
        
        # Example 3: No piecewise splits
        result3 = generate_unifying_proof(cases, "LAMBDA (x: real): g(x)", None)
        # For 3-case without piecewise splits, should skip
        self.assertIn("skipping", result3.lower())


if __name__ == "__main__":
    unittest.main()

