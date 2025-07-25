import unittest
from pvs_utils import ProofNode, ProofBuilder, ProofScript


class TestProofBuilder(unittest.TestCase):
    def test_spread_case_and_inst_where(self):
        pb = ProofBuilder()
        then = pb.create_then_sequence(
            ProofNode("SKEEP*"),
            ProofNode("SKOLETIN*"),
            ProofNode("FLATTEN"),
            ProofNode("SKEEP"),
            ProofNode("ASSERT"),
        )
        # Branch 1
        branch1 = pb.create_then_sequence(
            ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
            ProofNode("INST", ["-1", '"f"', '"xo + 2"', '"x"', '"0"']),
            ProofNode(
                "SPREAD",
                children=[
                    ProofNode("SPLIT", ["-1"]),
                    ProofNode(
                        "()",
                        children=[
                            pb.create_then_sequence(
                                ProofNode("ASSERT"),
                                ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
                                ProofNode(
                                    "INST", ["-1", '"f"', '"x"', '"xo-2"', '"0"']
                                ),
                                ProofNode("ASSERT"),
                                ProofNode("SKEEP"),
                                ProofNode("INST?"),
                                ProofNode("ASSERT"),
                            ),
                            ProofNode("PROPAX"),
                            ProofNode("ASSERT"),
                            pb.create_then_sequence(
                                ProofNode("SKEEP"),
                                ProofNode("INST?", [], {"WHERE": 1}),
                                ProofNode("ASSERT"),
                            ),
                        ],
                    ),
                ],
            ),
        )
        # Branch 2
        branch2 = pb.create_then_sequence(
            ProofNode("ASSERT"),
            ProofNode("EXPAND", ['"f"']),
            ProofNode("ASSERT"),
            ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
            ProofNode("INST", ["-1", '"f"', '"xo+2"', '"x"', '"0"']),
            ProofNode(
                "SPREAD",
                children=[
                    ProofNode("SPLIT", ["-1"]),
                    ProofNode(
                        "()",
                        children=[
                            pb.create_then_sequence(
                                ProofNode("ASSERT"),
                                ProofNode("EXPAND", ['"f"']),
                                ProofNode(
                                    "SPREAD",
                                    children=[
                                        ProofNode("CASE", ['"x=0"']),
                                        pb.create_then_sequence(
                                            ProofNode("ASSERT"),
                                            ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
                                            ProofNode(
                                                "INST",
                                                ["-1", '"f"', '"x"', '"0"', '"0"'],
                                            ),
                                            ProofNode("ASSERT"),
                                            ProofNode("SKEEP"),
                                            ProofNode("INST?"),
                                            ProofNode("ASSERT"),
                                        ),
                                    ],
                                ),
                            ),
                            ProofNode("EXPAND", ['"f"', "1"]),
                            ProofNode("PROPAX"),
                            ProofNode("ASSERT"),
                            pb.create_then_sequence(
                                ProofNode("SKEEP"),
                                ProofNode("INST?", [], {"WHERE": 1}),
                                ProofNode("ASSERT"),
                            ),
                        ],
                    ),
                ],
            ),
        )
        spread = pb.create_spread_case("xo-2 >=0", [branch1, branch2])
        then.add_child(spread)
        script = ProofScript("bound22_rect_function_bounded")
        script.root = then
        proof_text = script.generate()
        self.assertIn("SPREAD", proof_text)
        self.assertIn('(CASE "xo-2 >=0")', proof_text)
        self.assertIn("INST? :WHERE 1", proof_text)
        self.assertIn('EXPAND "f" 1', proof_text)
        self.assertIn("QED bound22_rect_function_bounded", proof_text)

    def test_bound22_rect_function_bounded_with_build_rect_proof(self):
        pb = ProofBuilder()
        # Branch 1: (THEN (LEMMA ...) (INST ...) (SPREAD (SPLIT -1) ...))
        branch1_inner1 = pb.create_then_sequence(
            ProofNode("ASSERT"),
            ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
            ProofNode("INST", ["-1", '"f"', '"x"', '"xo-2"', '"0"']),
            ProofNode("ASSERT"),
            ProofNode("SKEEP"),
            ProofNode("INST?"),
            ProofNode("ASSERT"),
        )
        branch1_inner2 = ProofNode("PROPAX")
        branch1_inner3 = ProofNode("ASSERT")
        branch1_inner4 = pb.create_then_sequence(
            ProofNode("SKEEP"),
            ProofNode("INST?", [], {"WHERE": 1}),
            ProofNode("ASSERT"),
        )
        branch1_spread = ProofNode(
            "SPREAD",
            children=[
                ProofNode("SPLIT", ["-1"]),
                ProofNode(
                    "()",
                    children=[
                        branch1_inner1,
                        branch1_inner2,
                        branch1_inner3,
                        branch1_inner4,
                    ],
                ),
            ],
        )
        branch1 = pb.create_then_sequence(
            ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
            ProofNode("INST", ["-1", '"f"', '"xo + 2"', '"x"', '"0"']),
            branch1_spread,
        )
        # Branch 2: (THEN (ASSERT) (EXPAND ...) ...)
        branch2_inner1 = pb.create_then_sequence(
            ProofNode("ASSERT"),
            ProofNode("EXPAND", ['"f"']),
            ProofNode(
                "SPREAD",
                children=[
                    ProofNode("CASE", ['"x=0"']),
                    pb.create_then_sequence(
                        ProofNode("ASSERT"),
                        ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
                        ProofNode("INST", ["-1", '"f"', '"x"', '"0"', '"0"']),
                        ProofNode("ASSERT"),
                        ProofNode("SKEEP"),
                        ProofNode("INST?"),
                        ProofNode("ASSERT"),
                    ),
                ],
            ),
        )
        branch2_inner2 = ProofNode("EXPAND", ['"f"', "1"])
        branch2_inner3 = ProofNode("PROPAX")
        branch2_inner4 = ProofNode("ASSERT")
        branch2_inner5 = pb.create_then_sequence(
            ProofNode("SKEEP"),
            ProofNode("INST?", [], {"WHERE": 1}),
            ProofNode("ASSERT"),
        )
        branch2_spread = ProofNode(
            "SPREAD",
            children=[
                ProofNode("SPLIT", ["-1"]),
                ProofNode(
                    "()",
                    children=[
                        branch2_inner1,
                        branch2_inner2,
                        branch2_inner3,
                        branch2_inner4,
                        branch2_inner5,
                    ],
                ),
            ],
        )
        branch2 = pb.create_then_sequence(
            ProofNode("ASSERT"),
            ProofNode("EXPAND", ['"f"']),
            ProofNode("ASSERT"),
            ProofNode("LEMMA", ['"mvt_gen_ge_bound"']),
            ProofNode("INST", ["-1", '"f"', '"xo+2"', '"x"', '"0"']),
            branch2_spread,
        )
        script = pb.build_rect_proof(
            lemma_name="bound22_rect_function_bounded",
            mvt_lemma="mvt_gen_ge_bound",
            inst_main=["f", "xo + 2", "x", "0"],
            spread_cases=[
                ("xo-2 >=0", [branch1, branch2]),
            ],
        )
        proof_text = script.generate()
        print(proof_text)
        self.assertIn("bound22_rect_function_bounded : PROOF", proof_text)
        self.assertIn("SPREAD", proof_text)
        self.assertIn('(CASE "xo-2 >=0")', proof_text)
        self.assertIn("INST? :WHERE 1", proof_text)
        self.assertIn("QED bound22_rect_function_bounded", proof_text)


if __name__ == "__main__":
    unittest.main()
