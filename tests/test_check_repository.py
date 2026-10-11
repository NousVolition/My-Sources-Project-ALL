"""Keep frozen historical plans distinct from current missing-link errors."""
from pathlib import Path
import tempfile
import unittest

from tools.check_repository import check


ARCHIVE = 'reports/clay-contact-stress/obstacle-only/MODEL-NOTES.md'


class RepositoryLinkTests(unittest.TestCase):
    def test_frozen_plan_does_not_hide_current_or_changed_document_errors(self):
        original = (Path(__file__).resolve().parents[1] / ARCHIVE).read_bytes()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            notes = root / ARCHIVE
            notes.parent.mkdir(parents=True)
            notes.write_bytes(original)
            current = root / 'README.md'
            current.write_text('[Current result](missing-current.html)\n',
                               encoding='utf-8')
            counts, errors = check(root)
            self.assertEqual(counts['preserved_historical_links'], 4)
            self.assertEqual(errors, [
                'README.md: missing link target missing-current.html'])

            current.unlink()
            notes.write_bytes(original + b'\n[New result](missing-new.html)\n')
            counts, errors = check(root)
            self.assertEqual(counts['preserved_historical_links'], 0)
            self.assertEqual(len(errors), 5)
            self.assertTrue(any('missing-new.html' in error for error in errors))


if __name__ == '__main__':
    unittest.main()
