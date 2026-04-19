"""Tests for src/utils.py — patient ID normalization and logger."""
import logging
import unittest


class TestNormalizePatientId(unittest.TestCase):
    """Tests for normalize_patient_id()."""

    def test_full_barcode_truncated_to_three_segments(self):
        """Full TCGA sample barcode is normalized to 12-char 3-segment patient ID."""
        from src.utils import normalize_patient_id
        result = normalize_patient_id("TCGA-BH-A18H-01A-11R-A115-07")
        self.assertEqual(result, "TCGA-BH-A18H")

    def test_already_three_segments_unchanged(self):
        """Barcode with exactly 3 segments is returned as-is."""
        from src.utils import normalize_patient_id
        result = normalize_patient_id("TCGA-E2-A14P")
        self.assertEqual(result, "TCGA-E2-A14P")

    def test_four_segment_barcode_truncated(self):
        """Four-segment barcode is truncated to first 3 segments."""
        from src.utils import normalize_patient_id
        result = normalize_patient_id("TCGA-AN-A04A-01A")
        self.assertEqual(result, "TCGA-AN-A04A")


class TestGetLogger(unittest.TestCase):
    """Tests for get_logger()."""

    def test_get_logger_returns_logger_with_correct_name(self):
        """get_logger returns a Logger with the specified name."""
        from src.utils import get_logger
        logger = get_logger("my_module")
        self.assertIsInstance(logger, logging.Logger)
        self.assertEqual(logger.name, "my_module")

    def test_get_logger_has_stream_handler(self):
        """Logger returned by get_logger has at least one StreamHandler."""
        from src.utils import get_logger
        logger = get_logger("test_handler_check")
        has_stream_handler = any(
            isinstance(h, logging.StreamHandler) for h in logger.handlers
        )
        self.assertTrue(has_stream_handler)


if __name__ == "__main__":
    unittest.main()
