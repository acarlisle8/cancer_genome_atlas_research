"""Tests for src/gdc_client.py — GDC manifest fetch and S3 download."""
import pathlib
import unittest
from unittest.mock import MagicMock, patch, mock_open


class TestFetchManifest(unittest.TestCase):
    """Tests for fetch_manifest()."""

    def _make_mock_response(self, hits):
        """Helper: build a mock requests.Response with given hits list."""
        mock_resp = MagicMock()
        mock_resp.raise_for_status.return_value = None
        mock_resp.json.return_value = {"data": {"hits": hits}}
        return mock_resp

    def test_fetch_manifest_returns_list_of_dicts(self):
        """fetch_manifest returns list of dicts with file_id, file_name, patient_id."""
        hits = [
            {
                "file_id": "744a6d3d-b666-49aa-8d26-47f34e3d1eb5",
                "file_name": "94027f46.rna_seq.augmented_star_gene_counts.tsv",
                "cases": [{"submitter_id": "TCGA-BH-A18H-01A-11R-A115-07"}],
            },
            {
                "file_id": "aabbccdd-1234-5678-abcd-ef1234567890",
                "file_name": "another_file.tsv",
                "cases": [{"submitter_id": "TCGA-E2-A14P-01A-22R-A115-07"}],
            },
        ]
        mock_resp = self._make_mock_response(hits)

        with patch("requests.post", return_value=mock_resp) as mock_post:
            from src.gdc_client import fetch_manifest
            result = fetch_manifest("TCGA-BRCA", "Transcriptome Profiling", "Gene Expression Quantification")

        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

        first = result[0]
        self.assertEqual(first["file_id"], "744a6d3d-b666-49aa-8d26-47f34e3d1eb5")
        self.assertEqual(first["file_name"], "94027f46.rna_seq.augmented_star_gene_counts.tsv")
        self.assertEqual(first["patient_id"], "TCGA-BH-A18H")

        second = result[1]
        self.assertEqual(second["patient_id"], "TCGA-E2-A14P")

    def test_fetch_manifest_empty_hits_returns_empty_list(self):
        """fetch_manifest returns [] when API response has no hits."""
        mock_resp = self._make_mock_response([])

        with patch("requests.post", return_value=mock_resp):
            from src.gdc_client import fetch_manifest
            result = fetch_manifest("TCGA-BRCA", "Transcriptome Profiling", "Gene Expression Quantification")

        self.assertEqual(result, [])

    def test_fetch_manifest_raises_on_http_error(self):
        """fetch_manifest raises requests.HTTPError on non-200 response."""
        import requests
        mock_resp = MagicMock()
        mock_resp.raise_for_status.side_effect = requests.HTTPError("404 Client Error")

        with patch("requests.post", return_value=mock_resp):
            from src.gdc_client import fetch_manifest
            with self.assertRaises(requests.HTTPError):
                fetch_manifest("INVALID-PROJECT", "Transcriptome Profiling", "Gene Expression Quantification")

    def test_fetch_manifest_posts_correct_payload(self):
        """fetch_manifest calls POST with correct GDC filter structure."""
        mock_resp = self._make_mock_response([])

        with patch("requests.post", return_value=mock_resp) as mock_post:
            from src.gdc_client import fetch_manifest
            fetch_manifest("TCGA-BRCA", "Transcriptome Profiling", "Gene Expression Quantification")

        mock_post.assert_called_once()
        call_kwargs = mock_post.call_args
        # Verify JSON body was passed
        json_body = call_kwargs[1].get("json") or call_kwargs[0][1] if len(call_kwargs[0]) > 1 else None
        if json_body is None and "json" in call_kwargs[1]:
            json_body = call_kwargs[1]["json"]
        self.assertIsNotNone(json_body)
        self.assertIn("filters", json_body)
        self.assertIn("fields", json_body)


class TestDownloadFile(unittest.TestCase):
    """Tests for download_file()."""

    def test_download_file_skips_existing_file(self):
        """download_file returns path immediately if file already exists (D-07 resume)."""
        with patch("s3fs.S3FileSystem") as mock_s3_cls:
            with patch("pathlib.Path.exists", return_value=True):
                from src.gdc_client import download_file
                dest_dir = pathlib.Path("/tmp/test_output")
                result = download_file(
                    "744a6d3d-b666-49aa-8d26-47f34e3d1eb5",
                    "some_file.tsv",
                    dest_dir,
                )

            # S3FileSystem should NOT have been instantiated
            mock_s3_cls.assert_not_called()

        expected_path = dest_dir / "some_file.tsv"
        self.assertEqual(result, expected_path)

    def test_download_file_downloads_from_s3_when_not_exists(self):
        """download_file uses s3fs with anon=True to download when file absent."""
        mock_s3_instance = MagicMock()
        mock_s3_file = MagicMock()
        mock_s3_file.__enter__ = MagicMock(return_value=mock_s3_file)
        mock_s3_file.__exit__ = MagicMock(return_value=False)
        mock_s3_file.read.return_value = b"file content"
        mock_s3_instance.open.return_value = mock_s3_file

        with patch("s3fs.S3FileSystem", return_value=mock_s3_instance) as mock_s3_cls:
            with patch("pathlib.Path.exists", return_value=False):
                with patch("builtins.open", mock_open()) as mock_file:
                    from src.gdc_client import download_file
                    dest_dir = pathlib.Path("/tmp/test_output")
                    result = download_file(
                        "744a6d3d-b666-49aa-8d26-47f34e3d1eb5",
                        "some_file.tsv",
                        dest_dir,
                    )

        # Verify anon=True was used
        mock_s3_cls.assert_called_once_with(anon=True)

        # Verify the correct S3 path was opened
        expected_s3_path = "s3://tcga-2-open/744a6d3d-b666-49aa-8d26-47f34e3d1eb5/some_file.tsv"
        mock_s3_instance.open.assert_called_once_with(expected_s3_path, "rb")

        expected_local_path = dest_dir / "some_file.tsv"
        self.assertEqual(result, expected_local_path)

    def test_download_file_returns_path_object(self):
        """download_file return value is a pathlib.Path."""
        with patch("pathlib.Path.exists", return_value=True):
            from src.gdc_client import download_file
            dest_dir = pathlib.Path("/tmp/test_output")
            result = download_file("abc", "file.tsv", dest_dir)

        self.assertIsInstance(result, pathlib.Path)


if __name__ == "__main__":
    unittest.main()
