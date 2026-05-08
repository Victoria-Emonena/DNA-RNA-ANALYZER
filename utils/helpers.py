"""
utils/helpers.py — Utility functions for DNA/RNA Analyzer
"""
import os
import io
import csv
from datetime import datetime


def save_uploaded_sequence(uploaded_file, upload_dir: str) -> tuple[str, str]:
    """Save an uploaded text/FASTA file and return (filename, content)."""
    os.makedirs(upload_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{timestamp}_{uploaded_file.name}"
    filepath = os.path.join(upload_dir, filename)
    content = uploaded_file.read().decode("utf-8", errors="ignore")
    with open(filepath, "w") as f:
        f.write(content)
    return filename, content


def records_to_csv(records: list[dict]) -> str:
    if not records:
        return ""
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=records[0].keys())
    writer.writeheader()
    writer.writerows(records)
    return output.getvalue()


def truncate_seq(seq: str, max_len: int = 60) -> str:
    if len(seq) <= max_len:
        return seq
    return seq[:max_len] + f"... (+{len(seq)-max_len} more)"
