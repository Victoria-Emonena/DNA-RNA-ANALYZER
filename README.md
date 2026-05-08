# 🧬 DNA & RNA Sequence Analyzer

**CSC 442 — Computational Biology & Interdisciplinary Studies | Project 2**

A production-ready multi-page Streamlit web application for analysing DNA and RNA sequences, performing transcription & translation, searching the UniProt protein database, and visualising sequence statistics.

---

## 📌 Features

### 🔬 Sequence Analyzer
- Paste sequence or upload `.txt` / `.fasta` files (drag-and-drop supported)
- **Auto-detection** of DNA, RNA, or invalid sequences with clear explanation
- **Strand type selection** — Template or Coding strand
- **Transcription** — DNA → mRNA with step-by-step explanation
- **Translation** — mRNA → Codons → Amino Acids (full name, 3-letter, 1-letter)
- **Polypeptide display** with beginner-friendly explanation

### 🧫 Protein Search
- Integrates the **UniProt REST API** for real protein data
- Displays: protein name, organism, function, accession, link to UniProt page
- Results saved to SQLite

### 📊 Visualizations
- Base composition pie chart
- GC vs AT/AU content bar chart
- Amino acid frequency histogram
- Full sequence statistics table

### 📋 History
- All sequence analyses stored in SQLite
- Protein search history
- Search, browse, delete, and export CSV

---

