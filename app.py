"""
app.py — DNA & RNA Sequence Analyzer
A production-ready multi-page Streamlit application for CSC 442 Project 2.
"""

import os
import sys
import json
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, os.path.dirname(__file__))

from database.db import (
    init_db, save_sequence, save_protein_search,
    get_all_sequences, get_protein_searches, delete_sequence,
)
from modules.analyzer import (
    clean_sequence, detect_sequence_type, dna_to_mrna,
    translate_mrna, rna_to_protein_chain, compute_stats, parse_fasta,
)
from modules.protein_api import search_uniprot
from utils.helpers import save_uploaded_sequence, records_to_csv, truncate_seq

UPLOAD_DIR = os.path.join(os.path.dirname(__file__), "uploads")

# ── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="DNA & RNA Analyzer",
    page_icon="🧬",
    layout_size="wide",
    initial_sidebar_state="expanded",
)

# ── CSS ───────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
    .main-title {
        font-size: 2.1rem;
        font-weight: 700;
        background: linear-gradient(90deg, #1A2332, #00BCD4);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }
    .subtitle {color: #A5D6A7; font-size: 1rem; margin-bottom: 1.5rem;}
    .seq-display {
        font-family: 'Courier New', monospace;
        font-size: 1.1rem;
        background: #4CAF50;
        border: 1px solid #E0F7FA;
        border-radius: 5px;
        padding: 1rem;
        word-break: break-all;
        line-height: 1.7;
        color: #263850;
    }
    .badge-dna {background:#1565C0; color:white; padding:3px 10px; border-radius:20px; font-size:0.8rem;}
    .badge-rna {background:#558B2F; color:white; padding:3px 10px; border-radius:20px; font-size:0.8rem;}
    .badge-invalid {background:#B71C1C; color:white; padding:3px 10px; border-radius:20px; font-size:0.8rem;}
    .codon-card {
        display:inline-block;
        background:#1E2130;
        border:1px solid #37474F;
        border-radius:6px;
        padding:4px 8px;
        margin:3px;
        font-family:monospace;
        font-size:0.8rem;
        text-align:center;
        min-width:60px;
    }
    .codon-stop {border-color:#EF5350;}
    .codon-start {border-color:#66BB6A;}
    .protein-chain {
        font-family: monospace;
        font-size: 1.1rem;
        letter-spacing: 0.15rem;
        background: #1A2332;
        padding: 0.8rem;
        border-radius: 8px;
        word-break: break-all;
        color: #80DEEA;
    }
    .stat-box {
        background: #1E2130;
        border-radius: 10px;
        padding: 1rem;
        text-align: center;
        border: 1px solid #2D3250;
    }
    .result-section {
        background: #0F2027;
        border-left: 4px solid #4CAF50;
        border-radius: 6px;
        padding: 1rem 1.2rem;
        margin: 0.7rem 0;
    }
    .protein-result {
        background: #1A2332;
        border: 1px solid #26A69A;
        border-radius: 10px;
        padding: 1rem;
        margin: 0.5rem 0;
    }
</style>
""", unsafe_allow_html=True)

# ── Init DB & session ────────────────────────────────────────────────────────
init_db()

if "analysis" not in st.session_state:
    st.session_state.analysis = None
if "sequence_id" not in st.session_state:
    st.session_state.sequence_id = None

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## DNA/RNA Analyzer")
    st.markdown("---")
    page = st.radio(
        "Navigate",
        [" Sequence Analyzer", " Protein Search", "Visualizations", "History"],
        label_visibility="collapsed",
    )
    st.markdown("---")
    st.markdown("### Quick Reference")
    with st.expander("DNA Bases"):
        st.markdown("- **A** — Adenine\n- **T** — Thymine\n- **G** — Guanine\n- **C** — Cytosine")
    with st.expander("RNA Bases"):
        st.markdown("- **A** — Adenine\n- **U** — Uracil\n- **G** — Guanine\n- **C** — Cytosine")
    with st.expander("Central Dogma"):
        st.markdown("DNA → mRNA (Transcription) → Protein (Translation)")

# ═══════════════════════════════════════════════════════
# PAGE 1 — SEQUENCE ANALYZER
# ═══════════════════════════════════════════════════════
if page == "Sequence Analyzer":
    st.markdown('<div class="main-title">DNA & RNA Sequence Analyzer</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">CSC 442 — Computational Biology | Project 2</div>', unsafe_allow_html=True)

    tab_input, tab_results = st.tabs(["Input", " Results"])

    with tab_input:
        col_left, col_right = st.columns([1, 1], gap="large")

        with col_left:
            st.markdown("#### Enter Sequence")
            seq_text = st.text_area(
                "Paste a DNA or RNA sequence (plain or FASTA format)",
                height=180,
                placeholder=">Example\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
            )

            st.markdown("#### Or Upload a File")
            uploaded_seq = st.file_uploader("Upload .txt or .fasta file", type=["txt", "fasta", "fa"])
            if uploaded_seq:
                _, file_content = save_uploaded_sequence(uploaded_seq, UPLOAD_DIR)
                seq_text = file_content
                st.success(f"File loaded: {uploaded_seq.name}")

        with col_right:
            st.markdown("#### Strand Configuration")
            strand_type = st.radio(
                "If DNA, which strand are you entering?",
                ["Template", "Coding"],
                help="Template strand: RNA polymerase reads this directly.\nCoding strand: same sequence as the mRNA (with T instead of U).",
            )
            st.info(
                "**Template strand**: RNA polymerase reads it. mRNA = complement (T→U).\n\n"
                "**Coding strand**: Same as mRNA sequence (T→U substitution only)."
            )

        analyze_btn = st.button("Analyse Sequence", use_container_width=True, type="primary")

    with tab_results:
        if analyze_btn and seq_text.strip():
            # Parse FASTA if needed
            parsed = parse_fasta(seq_text)
            if parsed:
                _, raw_seq = parsed[0]
            else:
                raw_seq = clean_sequence(seq_text)

            seq = clean_sequence(raw_seq) if raw_seq else clean_sequence(seq_text)

            seq_type, type_explanation = detect_sequence_type(seq)
            stats = compute_stats(seq, seq_type)

            mrna = ""
            mrna_explanation = ""
            codons = []
            protein = ""

            if seq_type == "DNA":
                mrna, mrna_explanation = dna_to_mrna(seq, strand_type)
                codons = translate_mrna(mrna)
                protein = rna_to_protein_chain(codons)
            elif seq_type == "RNA":
                mrna = seq
                mrna_explanation = "Input is already RNA — used directly as mRNA for translation."
                codons = translate_mrna(mrna)
                protein = rna_to_protein_chain(codons)

            # Save to DB
            seq_id = save_sequence(
                sequence_input=seq[:2000],
                sequence_type=seq_type,
                strand_type=strand_type if seq_type == "DNA" else "N/A",
                mrna=mrna[:2000],
                protein=protein[:500],
                gc_content=stats["gc_content"],
                length=stats["length"],
                source="manual" if not uploaded_seq else uploaded_seq.name,
            )

            st.session_state.analysis = {
                "seq": seq,
                "seq_type": seq_type,
                "type_explanation": type_explanation,
                "strand_type": strand_type,
                "mrna": mrna,
                "mrna_explanation": mrna_explanation,
                "codons": codons,
                "protein": protein,
                "stats": stats,
            }
            st.session_state.sequence_id = seq_id

        if st.session_state.analysis:
            a = st.session_state.analysis

            # ── Detection ──
            badge = {"DNA": "badge-dna", "RNA": "badge-rna", "INVALID": "badge-invalid"}[a["seq_type"]]
            st.markdown(f'<span class="{badge}">{a["seq_type"]}</span>', unsafe_allow_html=True)
            st.markdown(f'<div class="result-section">{a["type_explanation"]}</div>', unsafe_allow_html=True)

            # Stats row
            s = a["stats"]
            c1, c2, c3, c4 = st.columns(4)
            with c1:
                st.markdown(f'<div class="stat-box"><div style="font-size:1.8rem;font-weight:800;color:#64B5F6">{s["length"]}</div><div>Length (bp)</div></div>', unsafe_allow_html=True)
            with c2:
                st.markdown(f'<div class="stat-box"><div style="font-size:1.8rem;font-weight:800;color:#81C784">{s["gc_content"]}%</div><div>GC Content</div></div>', unsafe_allow_html=True)
            with c3:
                n_codons = len(a["codons"])
                st.markdown(f'<div class="stat-box"><div style="font-size:1.8rem;font-weight:800;color:#FFB74D">{n_codons}</div><div>Codons</div></div>', unsafe_allow_html=True)
            with c4:
                aa_count = len(a["protein"])
                st.markdown(f'<div class="stat-box"><div style="font-size:1.8rem;font-weight:800;color:#CE93D8">{aa_count}</div><div>Amino Acids</div></div>', unsafe_allow_html=True)

            st.markdown("---")

            # ── Input sequence ──
            with st.expander("Input Sequence", expanded=False):
                st.markdown(f'<div class="seq-display">{truncate_seq(a["seq"], 500)}</div>', unsafe_allow_html=True)

            # ── Transcription ──
            if a["seq_type"] == "DNA" and a["mrna"]:
                with st.expander("Transcription (DNA → mRNA)", expanded=True):
                    st.markdown(f'<div class="result-section">{a["mrna_explanation"]}</div>', unsafe_allow_html=True)
                    st.markdown("**mRNA Sequence:**")
                    st.markdown(f'<div class="seq-display">{truncate_seq(a["mrna"], 500)}</div>', unsafe_allow_html=True)

            # ── Translation ──
            if a["codons"]:
                with st.expander(" Translation (mRNA → Protein)", expanded=True):
                    st.markdown("**Codons & Amino Acids:**")

                    # Codon cards
                    cards_html = ""
                    for c in a["codons"]:
                        css = "codon-card"
                        if c["is_stop"]:
                            css += " codon-stop"
                        elif c["codon"] == "AUG":
                            css += " codon-start"
                        cards_html += (
                            f'<div class="{css}">'
                            f'<div style="color:#80DEEA">{c["codon"]}</div>'
                            f'<div style="font-size:0.7rem;color:#B0BEC5">{c["amino_acid_3"]}</div>'
                            f'<div style="font-size:0.65rem;color:#78909C">{c["amino_acid_1"]}</div>'
                            f'</div>'
                        )
                    st.markdown(cards_html, unsafe_allow_html=True)

                    st.markdown("---")
                    st.markdown("**Full Codon Table:**")
                    rows = [
                        {"#": i+1, "Codon": c["codon"], "Amino Acid": c["amino_acid_full"],
                         "3-letter": c["amino_acid_3"], "1-letter": c["amino_acid_1"],
                         "Type": "Stop" if c["is_stop"] else ("Start" if c["codon"] == "AUG" else "—")}
                        for i, c in enumerate(a["codons"])
                    ]
                    st.dataframe(rows, use_container_width=True, hide_index=True)

            # ── Polypeptide ──
            if a["protein"]:
                with st.expander("🧬 Polypeptide Chain", expanded=True):
                    st.markdown("**Protein (1-letter code):**")
                    st.markdown(f'<div class="protein-chain">{a["protein"]}</div>', unsafe_allow_html=True)

                    st.markdown("""
                    <div class="result-section" style="margin-top:0.8rem">
                    <b>What does this mean?</b><br>
                    Each letter represents one <b>amino acid</b> in the protein chain.
                    They are joined by peptide bonds to form a polypeptide. This sequence
                    folds into a 3D shape that determines the protein's function in the cell.
                    </div>
                    """, unsafe_allow_html=True)
        elif not analyze_btn:
            st.info("Enter a sequence in the **Input** tab and click **Analyse Sequence**.")

# ═══════════════════════════════════════════════════════
# PAGE 2 — PROTEIN SEARCH
# ═══════════════════════════════════════════════════════
elif page == "Protein Search":
    st.markdown('<div class="main-title">Protein Database Search</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">Search the UniProt protein database</div>', unsafe_allow_html=True)

    col_q, col_btn = st.columns([4, 1])
    with col_q:
        protein_query = st.text_input("Search for a protein (name, function, or gene symbol)", placeholder="e.g. insulin, hemoglobin, BRCA1")
    with col_btn:
        st.markdown("<br>", unsafe_allow_html=True)
        search_btn = st.button(" Search", use_container_width=True, type="primary")

    if search_btn and protein_query.strip():
        with st.spinner("Querying UniProt..."):
            results = search_uniprot(protein_query.strip(), max_results=8)

        if not results or (len(results) == 1 and "error" in results[0]):
            err = results[0].get("error", "Unknown error") if results else "No results"
            st.warning(f"Could not retrieve results: {err}")
            st.info(" Tip: Check your internet connection. UniProt API is required for this feature.")
        else:
            st.success(f"Found {len(results)} result(s) from UniProt")
            for i, res in enumerate(results):
                with st.expander(f"#{i+1} — {res['name']} ({res['accession']})", expanded=(i == 0)):
                    c1, c2 = st.columns(2)
                    with c1:
                        st.markdown(f"**Protein Name:** {res['name']}")
                        st.markdown(f"**Organism:** *{res['organism']}*")
                        st.markdown(f"**Accession:** `{res['accession']}`")
                    with c2:
                        st.markdown(f"**Function:**")
                        st.markdown(res["function"])
                    st.markdown(f"[ View on UniProt]({res['url']})")

                    # Save to DB if we have a linked sequence
                    if st.session_state.sequence_id:
                        save_protein_search(
                            sequence_id=st.session_state.sequence_id,
                            query_protein=protein_query,
                            result_name=res["name"],
                            result_organism=res["organism"],
                            result_function=res["function"][:200],
                            accession=res["accession"],
                        )
    else:
        st.info("Enter a protein name or keyword and click **Search** to query UniProt.")

# ═══════════════════════════════════════════════════════
# PAGE 3 — VISUALIZATIONS
# ═══════════════════════════════════════════════════════
elif page == "Visualizations":
    st.markdown('<div class="main-title">Sequence Visualizations</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">Charts and statistics for the last analysed sequence</div>', unsafe_allow_html=True)

    if not st.session_state.analysis:
        st.info("No analysis yet. Go to **Sequence Analyzer** to analyse a sequence first.")
    else:
        a = st.session_state.analysis
        stats = a["stats"]

        col1, col2 = st.columns(2)

        # Base composition pie
        with col1:
            st.markdown("#### Base Composition")
            counts = stats["counts"]
            if counts:
                labels = list(counts.keys())
                values = [counts[b] for b in labels]
                colors = {"A": "#4CAF50", "T": "#2196F3", "G": "#FF9800", "C": "#E91E63",
                          "U": "#9C27B0", "N": "#607D8B"}
                pie_colors = [colors.get(b, "#78909C") for b in labels]

                fig, ax = plt.subplots(figsize=(5, 4))
                fig.patch.set_facecolor("#0E1117")
                ax.set_facecolor("#0E1117")
                wedges, texts, autotexts = ax.pie(
                    values, labels=labels, autopct="%1.1f%%",
                    colors=pie_colors, startangle=90,
                    textprops={"color": "white"},
                )
                for at in autotexts:
                    at.set_color("white")
                ax.set_title(f"Base Composition ({a['seq_type']})", color="white", pad=10)
                st.pyplot(fig)
                plt.close(fig)

        # GC vs AT bar
        with col2:
            st.markdown("#### GC vs AT/AU Content")
            counts = stats["counts"]
            gc = counts.get("G", 0) + counts.get("C", 0)
            at = counts.get("A", 0) + counts.get("T", 0) + counts.get("U", 0)
            total = gc + at

            fig, ax = plt.subplots(figsize=(5, 4))
            fig.patch.set_facecolor("#0E1117")
            ax.set_facecolor("#1A1F2E")
            bars = ax.bar(
                ["GC Content", "AT/AU Content"],
                [gc / total * 100 if total else 0, at / total * 100 if total else 0],
                color=["#26C6DA", "#AB47BC"],
                width=0.5,
            )
            ax.set_ylabel("Percentage (%)", color="white")
            ax.set_title("GC vs AT/AU", color="white")
            ax.tick_params(colors="white")
            ax.spines["bottom"].set_color("#37474F")
            ax.spines["left"].set_color("#37474F")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2.0, height + 0.5,
                        f"{height:.1f}%", ha="center", va="bottom", color="white", fontsize=10)
            st.pyplot(fig)
            plt.close(fig)

        # Amino acid frequency
        if a["protein"]:
            st.markdown("#### Amino Acid Frequency")
            aa_counts = {}
            for aa in a["protein"]:
                aa_counts[aa] = aa_counts.get(aa, 0) + 1

            sorted_aas = sorted(aa_counts.items(), key=lambda x: -x[1])
            labels = [x[0] for x in sorted_aas]
            vals = [x[1] for x in sorted_aas]

            fig, ax = plt.subplots(figsize=(max(8, len(labels) * 0.6), 4))
            fig.patch.set_facecolor("#0E1117")
            ax.set_facecolor("#1A1F2E")
            ax.bar(labels, vals, color="#4CAF50", width=0.7)
            ax.set_xlabel("Amino Acid (1-letter)", color="white")
            ax.set_ylabel("Frequency", color="white")
            ax.set_title("Amino Acid Distribution", color="white")
            ax.tick_params(colors="white")
            ax.spines["bottom"].set_color("#37474F")
            ax.spines["left"].set_color("#37474F")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            st.pyplot(fig)
            plt.close(fig)

        # Sequence statistics table
        st.markdown("#### Sequence Statistics")
        stat_rows = [
            {"Metric": "Sequence Type", "Value": a["seq_type"]},
            {"Metric": "Strand Type", "Value": a["strand_type"] if a["seq_type"] == "DNA" else "N/A"},
            {"Metric": "Total Length (bp)", "Value": stats["length"]},
            {"Metric": "GC Content", "Value": f"{stats['gc_content']}%"},
            {"Metric": "Number of Codons", "Value": len(a["codons"])},
            {"Metric": "Amino Acids (excl. stop)", "Value": len(a["protein"])},
        ]
        for base, cnt in sorted(stats["counts"].items()):
            pct = (cnt / stats["length"] * 100) if stats["length"] else 0
            stat_rows.append({"Metric": f"Base {base} count", "Value": f"{cnt} ({pct:.1f}%)"})

        st.dataframe(stat_rows, use_container_width=True, hide_index=True)

# ═══════════════════════════════════════════════════════
# PAGE 4 — HISTORY
# ═══════════════════════════════════════════════════════
elif page == "History":
    st.markdown('<div class="main-title">Analysis History</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">All saved sequence analyses and protein searches</div>', unsafe_allow_html=True)

    tab_seq, tab_prot = st.tabs(["Sequences", "Protein Searches"])

    with tab_seq:
        search = st.text_input("Search sequences", placeholder="Search by type, protein...")
        records = get_all_sequences(search)

        if not records:
            st.info("No sequence records yet.")
        else:
            csv_data = records_to_csv(records)
            st.download_button("Export CSV", csv_data, "sequences.csv", "text/csv")
            st.markdown(f"**{len(records)} record(s)**")

            for rec in records:
                with st.expander(f"#{rec['id']} | {rec['sequence_type']} | {rec['timestamp'][:19]}"):
                    c1, c2 = st.columns(2)
                    with c1:
                        st.markdown(f"**Type:** {rec['sequence_type']}")
                        st.markdown(f"**Strand:** {rec['strand_type']}")
                        st.markdown(f"**Length:** {rec['length']} bp")
                        st.markdown(f"**GC Content:** {rec['gc_content']}%")
                    with c2:
                        st.markdown(f"**Source:** {rec['source']}")
                        st.markdown(f"**Protein:** {truncate_seq(rec['protein'] or '', 40)}")
                        st.markdown(f"**mRNA (preview):** {truncate_seq(rec['mrna'] or '', 40)}")
                    if st.button("Delete", key=f"delseq_{rec['id']}"):
                        delete_sequence(rec["id"])
                        st.success("Deleted.")
                        st.rerun()

    with tab_prot:
        prot_records = get_protein_searches()
        if not prot_records:
            st.info("No protein search history yet.")
        else:
            csv_p = records_to_csv(prot_records)
            st.download_button("Export CSV", csv_p, "protein_searches.csv", "text/csv")
            st.markdown(f"**{len(prot_records)} search(es)**")
            for pr in prot_records:
                with st.expander(f"#{pr['id']} | {pr['query_protein']} → {pr['result_name']} | {pr['timestamp'][:19]}"):
                    st.markdown(f"**Query:** {pr['query_protein']}")
                    st.markdown(f"**Result:** {pr['result_name']}")
                    st.markdown(f"**Organism:** {pr['result_organism']}")
                    st.markdown(f"**Accession:** {pr['accession']}")
                    st.markdown(f"**Function:** {pr['result_function']}")
