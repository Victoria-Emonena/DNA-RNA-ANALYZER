"""
database/db.py — SQLite database for DNA & RNA Sequence Analyzer
"""
import sqlite3
import os
from datetime import datetime

DB_PATH = os.path.join(os.path.dirname(__file__), "dna_rna.db")


def get_connection():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def init_db():
    conn = get_connection()
    cur = conn.cursor()

    cur.execute("""
        CREATE TABLE IF NOT EXISTS sequences (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sequence_input TEXT NOT NULL,
            sequence_type TEXT NOT NULL,
            strand_type TEXT,
            mrna TEXT,
            protein TEXT,
            gc_content REAL,
            length INTEGER,
            source TEXT,
            timestamp TEXT NOT NULL
        )
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS protein_searches (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sequence_id INTEGER,
            query_protein TEXT,
            result_name TEXT,
            result_organism TEXT,
            result_function TEXT,
            accession TEXT,
            timestamp TEXT NOT NULL,
            FOREIGN KEY (sequence_id) REFERENCES sequences(id)
        )
    """)

    conn.commit()
    conn.close()


def save_sequence(
    sequence_input: str,
    sequence_type: str,
    strand_type: str,
    mrna: str,
    protein: str,
    gc_content: float,
    length: int,
    source: str,
) -> int:
    conn = get_connection()
    try:
        cur = conn.execute(
            """
            INSERT INTO sequences
              (sequence_input, sequence_type, strand_type, mrna, protein, gc_content, length, source, timestamp)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (sequence_input, sequence_type, strand_type, mrna, protein, gc_content, length, source, datetime.now().isoformat()),
        )
        conn.commit()
        return cur.lastrowid
    finally:
        conn.close()


def save_protein_search(sequence_id: int, query_protein: str, result_name: str,
                        result_organism: str, result_function: str, accession: str):
    conn = get_connection()
    try:
        conn.execute(
            """
            INSERT INTO protein_searches
              (sequence_id, query_protein, result_name, result_organism, result_function, accession, timestamp)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (sequence_id, query_protein, result_name, result_organism, result_function, accession, datetime.now().isoformat()),
        )
        conn.commit()
    finally:
        conn.close()


def get_all_sequences(search: str = ""):
    conn = get_connection()
    try:
        if search:
            rows = conn.execute(
                """SELECT * FROM sequences WHERE sequence_type LIKE ? OR protein LIKE ?
                   ORDER BY timestamp DESC""",
                (f"%{search}%", f"%{search}%"),
            ).fetchall()
        else:
            rows = conn.execute("SELECT * FROM sequences ORDER BY timestamp DESC").fetchall()
        return [dict(r) for r in rows]
    finally:
        conn.close()


def get_protein_searches():
    conn = get_connection()
    try:
        rows = conn.execute("SELECT * FROM protein_searches ORDER BY timestamp DESC").fetchall()
        return [dict(r) for r in rows]
    finally:
        conn.close()


def delete_sequence(seq_id: int):
    conn = get_connection()
    try:
        conn.execute("DELETE FROM protein_searches WHERE sequence_id = ?", (seq_id,))
        conn.execute("DELETE FROM sequences WHERE id = ?", (seq_id,))
        conn.commit()
    finally:
        conn.close()
