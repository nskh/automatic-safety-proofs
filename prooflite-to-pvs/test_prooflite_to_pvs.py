#!/usr/bin/env python3
"""Tests for prooflite_to_pvs converter."""

import unittest
from prooflite_to_pvs import (
    strip_prooflite_prefix,
    remove_then_spread,
    balance_parentheses,
    convert_prooflite_to_pvs,
)


class TestStripProoflitePrefix(unittest.TestCase):
    def test_simple_prefix(self):
        text = "%|- (SKEEP*)"
        result = strip_prooflite_prefix(text)
        self.assertEqual(result, "(SKEEP*)")

    def test_prefix_with_leading_space(self):
        text = "%|-  (SPREAD (CASE \"test\"))"
        result = strip_prooflite_prefix(text)
        self.assertEqual(result, " (SPREAD (CASE \"test\"))")

    def test_multiline(self):
        text = "%|- (THEN (SKEEP*)\n%|-  (FLATTEN))"
        result = strip_prooflite_prefix(text)
        self.assertEqual(result, "(THEN (SKEEP*)\n (FLATTEN))")


class TestRemoveThenSpread(unittest.TestCase):
    def test_remove_then(self):
        text = "(THEN (SKEEP*) (FLATTEN))"
        result = remove_then_spread(text)
        self.assertEqual(result, "(SKEEP*) (FLATTEN))")

    def test_remove_spread(self):
        text = "(SPREAD (CASE \"test\") ((ASSERT) (PROPAX)))"
        result = remove_then_spread(text)
        self.assertEqual(result, "(CASE \"test\") ((ASSERT) (PROPAX)))")

    def test_nested(self):
        text = "(THEN (FLATTEN) (SPREAD (SPLIT -1) ((ASSERT) (PROPAX))))"
        result = remove_then_spread(text)
        self.assertEqual(result, "(FLATTEN) (SPLIT -1) ((ASSERT) (PROPAX))))")


class TestBalanceParentheses(unittest.TestCase):
    def test_simple_extra_close(self):
        text = "(SKEEP*))"
        result = balance_parentheses(text)
        self.assertEqual(result, "(SKEEP*)")

    def test_multiple_extra_close(self):
        text = "(ASSERT))))"
        result = balance_parentheses(text)
        self.assertEqual(result, "(ASSERT)")

    def test_string_with_parens(self):
        text = '(CASE "x >= 0 AND (y < 1)"))'
        result = balance_parentheses(text)
        self.assertEqual(result, '(CASE "x >= 0 AND (y < 1)")')

    def test_balanced_unchanged(self):
        text = "(SKEEP*) (FLATTEN)"
        result = balance_parentheses(text)
        self.assertEqual(result, "(SKEEP*) (FLATTEN)")


class TestFullConversion(unittest.TestCase):
    def test_simple_conversion(self):
        prooflite = "%|- (THEN (SKEEP*) (FLATTEN))"
        result = convert_prooflite_to_pvs(prooflite)
        self.assertEqual(result, "(SKEEP*) (FLATTEN)")

    def test_spread_conversion(self):
        prooflite = '%|- (SPREAD (CASE "test") ((ASSERT) (PROPAX)))'
        result = convert_prooflite_to_pvs(prooflite)
        self.assertEqual(result, '(CASE "test") ((ASSERT) (PROPAX))')

    def test_nested_conversion(self):
        prooflite = """%|- (THEN (SKEEP*)
%|-  (SPREAD (CASE "x >= 0")
%|-   ((THEN (FLATTEN) (ASSERT))
%|-    (PROPAX))))"""
        result = convert_prooflite_to_pvs(prooflite)
        # After conversion, should have balanced parens and no THEN/SPREAD
        self.assertNotIn("(THEN ", result)
        self.assertNotIn("(SPREAD ", result)
        # Should have balanced parens
        self.assertEqual(result.count("("), result.count(")"))


class TestSpecExample(unittest.TestCase):
    """Test with the example from the spec file."""

    def test_spec_example_structure(self):
        """Test that the conversion produces valid structure."""
        prooflite_input = """%|- (THEN (SKEEP*) (SKOLETIN*) (FLATTEN) (SKEEP)
%|-  (SPREAD (CASE "xo - 2 >= 0 AND xo + 2 <= 4")
%|-   ((THEN (FLATTEN) (ASSERT))
%|-    (SPREAD (CASE "xo - 2 < 0 AND xo + 2 <= 4")
%|-     ((THEN (FLATTEN) (LEMMA "mvt_gen_ge_ci")
%|-       (SPREAD (INST -1 "f" "0" "4" "0" "xo + 2" "x")
%|-        ((SPREAD (SPLIT -1)
%|-          ((THEN (LEMMA "mvt_gen_ge_ci")
%|-            (SPREAD (INST -1 "f" "0" "4" "0" "x" "0")
%|-             ((SPREAD (SPLIT -1)
%|-               ((THEN (LEMMA "mvt_gen_le_ci")
%|-                 (SPREAD (INST -1 "f" "0" "4" "8" "xo + 2" "x")
%|-                  ((SPREAD (SPLIT -1)
%|-                    ((THEN (ASSERT) (LEMMA "mvt_gen_le_ci")
%|-                      (SPREAD (INST -1 "f" "0" "4" "8" "x" "0")
%|-                       ((SPREAD (SPLIT -1)
%|-                         ((THEN (EXPAND "f") (ASSERT :FLUSH? T))
%|-                          (ASSERT :FLUSH? T) (PROPAX) (ASSERT :FLUSH? T)
%|-                          (PROPAX)))
%|-                        (THEN (EXPAND "ci") (PROPAX))
%|-                        (THEN (EXPAND "ci") (PROPAX)))))
%|-                     (ASSERT) (PROPAX) (PROPAX) (PROPAX)))
%|-                   (THEN (EXPAND "ci") (ASSERT))
%|-                   (THEN (EXPAND "ci") (ASSERT)))))
%|-                (ASSERT) (PROPAX) (PROPAX) (PROPAX)))
%|-              (THEN (EXPAND "ci") (PROPAX))
%|-              (THEN (ASSERT) (EXPAND "ci") (PROPAX)))))
%|-           (ASSERT) (PROPAX) (PROPAX) (PROPAX)))
%|-         (THEN (EXPAND "ci") (ASSERT)) (THEN (EXPAND "ci") (ASSERT)))))
%|-      (SPREAD (CASE "xo - 2 >= 0 AND xo + 2 > 4")
%|-       ((THEN (FLATTEN) (LEMMA "mvt_gen_ge_ci")
%|-         (SPREAD (INST -1 "f" "0" "4" "0" "4" "x")
%|-          ((SPREAD (SPLIT -1)
%|-            ((THEN (LEMMA "mvt_gen_ge_ci")
%|-              (SPREAD (INST -1 "f" "0" "4" "0" "x" "xo - 2")
%|-               ((SPREAD (SPLIT -1)
%|-                 ((THEN (LEMMA "mvt_gen_le_ci")
%|-                   (SPREAD (INST -1 "f" "0" "4" "8" "4" "x")
%|-                    ((SPREAD (SPLIT -1)
%|-                      ((THEN (LEMMA "mvt_gen_le_ci")
%|-                        (SPREAD (INST -1 "f" "0" "4" "8" "x" "xo - 2")
%|-                         ((SPREAD (SPLIT -1)
%|-                           ((ASSERT) (ASSERT) (PROPAX) (ASSERT) (PROPAX)))
%|-                          (THEN (EXPAND "ci") (ASSERT))
%|-                          (THEN (EXPAND "ci") (ASSERT)))))
%|-                       (ASSERT) (PROPAX) (ASSERT) (PROPAX)))
%|-                     (THEN (EXPAND "ci") (ASSERT))
%|-                     (THEN (EXPAND "ci") (PROPAX)))))
%|-                  (ASSERT) (PROPAX) (ASSERT) (PROPAX)))
%|-                (THEN (EXPAND "ci") (ASSERT)) (THEN (EXPAND "ci") (ASSERT)))))
%|-             (ASSERT) (PROPAX) (ASSERT) (PROPAX)))
%|-           (THEN (EXPAND "ci") (ASSERT)) (THEN (EXPAND "ci") (PROPAX)))))
%|-        (SPREAD (CASE "xo - 2 < 0 AND xo + 2 > 4")
%|-         ((SPREAD (CASE "xo=x")
%|-           ((ASSERT)
%|-            (SPREAD (CASE "xo>x")
%|-             ((THEN (ASSERT) (LEMMA "mvt_gen_ge_ci")
%|-               (SPREAD (INST -1 "f" "0" "4" "0" "xo" "x")
%|-                ((SPREAD (SPLIT -1)
%|-                  ((THEN (LEMMA "mvt_gen_ge_ci")
%|-                    (SPREAD (INST -1 "f" "0" "4" "0" "x" "0")
%|-                     ((SPREAD (SPLIT -1)
%|-                       ((THEN (LEMMA "mvt_gen_le_ci")
%|-                         (SPREAD (INST -1 "f" "0" "4" "8" "xo" "x")
%|-                          ((SPREAD (SPLIT -1)
%|-                            ((ASSERT)
%|-                             (THEN (LEMMA "mvt_gen_le_ci")
%|-                              (SPREAD (INST -1 "f" "0" "4" "8" "x" "0")
%|-                               ((SPREAD (SPLIT -1)
%|-                                 ((THEN (EXPAND "f") (ASSERT :FLUSH? T))
%|-                                  (ASSERT) (PROPAX) (ASSERT) (PROPAX)))
%|-                                (THEN (EXPAND "ci" 1) (PROPAX))
%|-                                (THEN (EXPAND "ci" 1) (PROPAX)))))
%|-                             (PROPAX) (ASSERT) (PROPAX)))
%|-                           (THEN (EXPAND "ci" 1) (PROPAX))
%|-                           (THEN (EXPAND "ci" 1) (PROPAX)))))
%|-                        (ASSERT) (PROPAX) (ASSERT) (PROPAX)))
%|-                      (THEN (EXPAND "ci" 1) (PROPAX))
%|-                      (THEN (ASSERT) (EXPAND "ci" 1) (PROPAX)))))
%|-                   (ASSERT) (PROPAX) (ASSERT) (PROPAX)))
%|-                 (THEN (EXPAND "ci" 1) (PROPAX))
%|-                 (THEN (EXPAND "ci" 1) (PROPAX)))))
%|-              (SPREAD (CASE "xo < x")
%|-               ((THEN (LEMMA "mvt_gen_ge_ci")
%|-                 (SPREAD (INST -1 "f" "0" "4" "0" "x" "xo")
%|-                  ((SPREAD (SPLIT -1)
%|-                    ((THEN (LEMMA "mvt_gen_ge_ci")
%|-                      (SPREAD (INST -1 "f" "0" "4" "0" "xo" "0")
%|-                       ((SPREAD (SPLIT -1)
%|-                         ((THEN (LEMMA "mvt_gen_le_ci")
%|-                           (SPREAD (INST -1 "f" "0" "4" "8" "x" "xo")
%|-                            ((SPREAD (SPLIT -1)
%|-                              ((ASSERT)
%|-                               (THEN (LEMMA "mvt_gen_le_ci")
%|-                                (SPREAD
%|-                                 (INST -1 "f" "0" "4" "8" "xo" "0")
%|-                                 ((SPREAD (SPLIT -1)
%|-                                   ((THEN (EXPAND "f") (ASSERT :FLUSH? T))
%|-                                    (ASSERT) (PROPAX) (ASSERT) (PROPAX)))
%|-                                  (THEN (EXPAND "ci" 1) (PROPAX))
%|-                                  (THEN (EXPAND "ci" 1) (ASSERT)))))
%|-                               (PROPAX) (ASSERT) (PROPAX)))
%|-                             (THEN (EXPAND "ci" 1) (ASSERT))
%|-                             (THEN (EXPAND "ci" 1) (ASSERT)))))
%|-                          (ASSERT) (PROPAX) (ASSERT) (PROPAX)))
%|-                        (THEN (EXPAND "ci" 1) (PROPAX))
%|-                        (THEN (EXPAND "ci" 1) (ASSERT)))))
%|-                     (ASSERT) (PROPAX) (ASSERT) (PROPAX)))
%|-                   (THEN (EXPAND "ci" 1) (ASSERT))
%|-                   (THEN (EXPAND "ci" 1) (ASSERT)))))
%|-                (ASSERT)))))))
%|-          (THEN (HIDE-ALL-BUT (1 2 3 4)) (GRIND)))))))"""

        result = convert_prooflite_to_pvs(prooflite_input)

        # Verify no THEN or SPREAD remain
        self.assertNotIn("(THEN ", result)
        self.assertNotIn("(SPREAD ", result)

        # Verify no %|- prefix remains
        self.assertNotIn("%|-", result)

        # Verify parentheses are balanced
        self.assertEqual(result.count("("), result.count(")"))

        # Verify key commands are present
        self.assertIn("(SKEEP*)", result)
        self.assertIn("(SKOLETIN*)", result)
        self.assertIn("(FLATTEN)", result)
        self.assertIn('(CASE "xo - 2 >= 0 AND xo + 2 <= 4")', result)
        self.assertIn("(LEMMA \"mvt_gen_ge_ci\")", result)
        self.assertIn("(SPLIT -1)", result)
        self.assertIn("(EXPAND \"f\")", result)
        self.assertIn("(ASSERT :FLUSH? T)", result)
        self.assertIn("(PROPAX)", result)
        self.assertIn("(HIDE-ALL-BUT (1 2 3 4))", result)
        self.assertIn("(GRIND)", result)


if __name__ == "__main__":
    unittest.main()
