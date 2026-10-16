"""Stdlib-only command and receipt-isolation checks; no research DNA or model."""
import ast
import json
import tempfile
import types
import unittest
from pathlib import Path

SOURCE = Path(__file__).with_name("run_carbon_a_train.py")


class ContextChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="m28-carbon-context-check-")
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        scope = [("fixture", x, 100000) for x in (
            "NC_003070.9", "NC_003071.7", "NC_003076.8", "NC_089043.1", "NC_089044.1", "NC_089046.1")]
        namespace = dict(ROOT=root, OUT=root / "default", CONTEXT24_OUT=root / "short", ASSETS=root / "assets",
            SCOPE=scope, REVISION="fixture", json=json, N=types.SimpleNamespace(INPUT_MODE="fixture_raw_fasta",
                input_fasta=lambda seqid: root / "inputs" / seqid))
        tree = ast.parse(SOURCE.read_text())
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                     and node.name in ("output_root", "pipeline_commands", "completed_runs")]
        exec(compile(ast.Module(body=functions, type_ignores=[]), str(SOURCE), "exec"), namespace)
        self.ns, self.root = namespace, root
        for context, out in ((98304, namespace["OUT"]), (24576, namespace["CONTEXT24_OUT"])):
            for task, (species, seqid, length) in enumerate(scope):
                folder = out / seqid
                folder.mkdir(parents=True)
                (folder / "STATUS").write_text("COMPLETED\n")
                gff = folder / "prediction.gff"
                gff.touch()
                row = dict(completed=True, native_inference_executed=True, task=task, species=species, seqid=seqid,
                    length=length, fixed_revision="fixture", input_fasta=str(namespace["N"].input_fasta(seqid)),
                    input_mode="fixture_raw_fasta", optimizer_steps=0, reference_used=False,
                    labels_or_feature_cache_used=False, DEV_used=False, sealed_data_used=False, test_setaria_used=False,
                    output_check={"coordinate_cds_protein_roundtrip": True}, effective_attention="flash_attention_2",
                    viterbi_backend="gpu", raw_GFF=str(gff.relative_to(root)), context_bp=context,
                    context_tokens=context // 6, overlap_bp=6144, overlap_tokens=1024, window_policy_modified=context != 98304)
                (folder / "summary.json").write_text(json.dumps(row))
                if context == 98304 and seqid == "NC_003071.7":
                    recovered = root / "outputs/M28-CARBON-A-TRAIN-RAW-VALIDATION" / seqid
                    recovered.mkdir(parents=True)
                    row.update(source_inference_job="13503402", source_inference_job_state="FAILED",
                        literal_standard_codon_table_with_non_ACGT_to_X=True,
                        native_GPU_inference_reexecuted=False, raw_outputs_modified=False)
                    (recovered / "STATUS").write_text("COMPLETED\n")
                    (recovered / "summary.json").write_text(json.dumps(row))

    def test_native_commands_and_units(self):
        make = self.ns["pipeline_commands"]
        default = make(self.root / "input", self.root / "output", "NC_003071.7", 98304)
        short = make(self.root / "input", self.root / "output", "NC_003071.7", 24576)
        self.assertEqual(len(default), 1)
        self.assertEqual(short[1][:-1], default[0])
        self.assertEqual(short[1][-1], "--skip-inference")
        args = short[0]
        self.assertEqual(int(args[args.index("--context_length") + 1]) * 6, 24576)
        self.assertEqual(int(args[args.index("--overlap_length") + 1]) * 6, 6144)
        self.assertNotIn("--context_length", default[0])

    def test_original_recovery_and_all_short_receipts_are_separate(self):
        original = self.ns["completed_runs"]()
        short = self.ns["completed_runs"](24576)
        self.assertEqual(len(original), 6)
        self.assertEqual(len(short), 6)
        self.assertIn("TRAIN-RAW-VALIDATION", original["NC_003071.7"]["accepted_run_receipt"])
        for row in short.values():
            self.assertTrue(row["accepted_run_receipt"].startswith("short/"))
            self.assertEqual(row["context_bp"], 24576)

    def test_wrong_condition_receipt_is_rejected(self):
        path = self.ns["CONTEXT24_OUT"] / "NC_003071.7/summary.json"
        row = json.loads(path.read_text())
        row["context_bp"] = 98304
        path.write_text(json.dumps(row))
        with self.assertRaises(AssertionError):
            self.ns["completed_runs"](24576)


if __name__ == "__main__":
    unittest.main(verbosity=2)
