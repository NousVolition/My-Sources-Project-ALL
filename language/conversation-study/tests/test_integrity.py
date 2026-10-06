"""Evidence integrity checks and deliberately damaged fixtures."""
import copy
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('conversation_evidence_validator', ROOT / 'validate.py')
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


class EvidenceIntegrity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.saved = validator.read_bundle(ROOT)
        cls.raw = {s['id']: (ROOT / 'sources' / f'{s["id"]}.txt').read_bytes() for s in cls.saved['manifest']['sources']}

    def setUp(self):
        self.data = copy.deepcopy(self.saved)

    def test_saved_evidence_is_consistent(self):
        result = validator.validate_data(self.data, self.raw)
        self.assertGreater(result['exact_fragments'], 0)

    def test_edited_source_is_detected(self):
        changed = dict(self.raw)
        changed['source-01'] += b'changed'
        with self.assertRaisesRegex(ValueError, 'Source hash mismatch'):
            validator.validate_data(self.data, changed)

    def test_inexact_quote_is_detected(self):
        self.data['word-paths']['paths'][0]['nodes'][0]['quote'] = 'This quotation was not in the source.'
        with self.assertRaisesRegex(ValueError, 'Quote mismatch'):
            validator.validate_data(self.data, self.raw)

    def test_duplicate_episode_is_detected(self):
        self.data['evidence']['events'].append(self.data['evidence']['events'][0])
        with self.assertRaisesRegex(ValueError, 'Duplicate episode'):
            validator.validate_data(self.data, self.raw)

    def test_incorrect_overlap_is_detected(self):
        self.data['manifest']['overlap'][0]['new_lines'][0] = 2
        with self.assertRaisesRegex(ValueError, 'overlap is not exact'):
            validator.validate_data(self.data, self.raw)

    def test_different_request_cannot_claim_exact_match(self):
        self.data['confusion-tracks']['same_request']['refs'][1] = ['source-03', 2211, 2211]
        with self.assertRaisesRegex(ValueError, 'Repeated request mismatch'):
            validator.validate_data(self.data, self.raw)

    def test_unknown_category_is_detected(self):
        self.data['evidence']['events'][0]['tags'].append('invented-category')
        with self.assertRaisesRegex(ValueError, 'Unknown category'):
            validator.validate_data(self.data, self.raw)

    def test_out_of_range_source_is_detected(self):
        self.data['evidence']['events'][0]['refs'][0] = ['source-01', 99999, 99999]
        with self.assertRaisesRegex(ValueError, 'Invalid source span'):
            validator.validate_data(self.data, self.raw)

    def test_interactive_links_and_source_displays(self):
        self.assertGreater(validator.validate_pages(ROOT), 0)


if __name__ == '__main__':
    unittest.main()
