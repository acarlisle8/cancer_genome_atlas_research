"""Shared helpers: patient ID normalization and logging setup."""
import logging


def normalize_patient_id(barcode: str) -> str:
    """
    Normalize a TCGA sample barcode to the 12-char patient ID.

    Takes the first 3 hyphen-delimited segments.

    Examples:
        "TCGA-BH-A18H-01A-11R-A115-07" -> "TCGA-BH-A18H"
        "TCGA-E2-A14P"                  -> "TCGA-E2-A14P"
        "TCGA-AN-A04A-01A"              -> "TCGA-AN-A04A"

    Args:
        barcode: Full TCGA sample barcode or already-normalized patient ID

    Returns:
        12-character TCGA patient ID (3 hyphen-delimited segments)
    """
    parts = barcode.split("-")
    return "-".join(parts[:3])


def get_logger(name: str) -> logging.Logger:
    """
    Return a logger with the given name. Adds a StreamHandler if none present.

    Format: "%(asctime)s %(levelname)s %(name)s: %(message)s"

    Args:
        name: Logger name (typically the module name)

    Returns:
        Configured logging.Logger instance
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
        )
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger
