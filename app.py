# app.py
# Streamlit app: Filter an MGF file by a SMILES list from a table (CSV/TSV).
# - Upload MGF
# - Upload table containing a SMILES column
# - Auto-detect delimiter and SMILES column
# - Stream-filter MGF blocks (memory-friendly)
# - Download filtered MGF
# - Plot SMILES as an RDKit structure grid (with compound_name if present)

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from typing import Iterable, Optional, Tuple, List, Set

import pandas as pd
import streamlit as st
from PIL import Image

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKit_AVAILABLE = True
except Exception:
    RDKit_AVAILABLE = False


# -----------------------------
# Utilities
# -----------------------------

COMMON_DELIMS = [",", "\t", ";", "|"]


def sniff_delimiter(file_bytes: bytes, default: str = ",") -> str:
    """Try to detect CSV delimiter. Falls back to default if sniffing fails."""
    sample = file_bytes[:8192].decode("utf-8", errors="ignore")
    if not sample.strip():
        return default

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=COMMON_DELIMS)
        return dialect.delimiter
    except Exception:
        # simple heuristic: choose delimiter with max occurrences in first lines
        lines = [ln for ln in sample.splitlines()[:20] if ln.strip()]
        if not lines:
            return default
        counts = {d: sum(ln.count(d) for ln in lines) for d in COMMON_DELIMS}
        best = max(counts, key=counts.get)
        return best if counts[best] > 0 else default


def detect_smiles_column(columns: Iterable[str]) -> Optional[str]:
    """Heuristic SMILES column detection."""
    cols = list(columns)
    lower_map = {c.lower().strip(): c for c in cols}

    # strong matches
    for key in ("smiles", "canonical_smiles", "compound_smiles", "structure_smiles"):
        if key in lower_map:
            return lower_map[key]

    # broader contains
    for c in cols:
        if "smiles" in c.lower():
            return c

    return None


def detect_compound_name_column(columns: Iterable[str]) -> Optional[str]:
    """Heuristic compound name column detection."""
    cols = list(columns)
    lower_map = {c.lower().strip(): c for c in cols}

    # strong matches
    for key in ("compound_name", "compoundname", "name", "compound"):
        if key in lower_map:
            return lower_map[key]

    # broader contains
    for c in cols:
        cl = c.lower()
        if "compound" in cl and "name" in cl:
            return c

    return None


def load_table_df(table_bytes: bytes) -> Tuple[pd.DataFrame, str]:
    """Load the full table as a dataframe, returning (df, detected_sep)."""
    sep = sniff_delimiter(table_bytes, default=",")
    buf = io.BytesIO(table_bytes)
    df = pd.read_csv(buf, sep=sep, engine="python")
    return df, sep


def load_smiles_set(table_bytes: bytes) -> Tuple[Set[str], str, str]:
    """
    Load a table (CSV/TSV) and return:
      - smiles_set
      - detected delimiter
      - detected column name used
    """
    df, sep = load_table_df(table_bytes)
    smiles_col = detect_smiles_column(df.columns)
    if not smiles_col:
        raise ValueError(
            "No SMILES column found. Expected a column named like: SMILES, compound_smiles, structure_smiles."
        )

    smiles = df[smiles_col].dropna().astype(str).map(lambda s: s.strip())
    smiles_set = {s for s in smiles if s}

    if not smiles_set:
        raise ValueError(f"SMILES column '{smiles_col}' was found, but it contains no valid values.")

    return smiles_set, sep, smiles_col


# -----------------------------
# Filtering core
# -----------------------------

@dataclass
class FilterReport:
    total_blocks: int
    matched_blocks: int
    blocks_with_no_smiles: int
    unique_smiles_in_table: int
    unique_smiles_matched: int


def iter_filtered_mgf_bytes(
    mgf_text_stream: io.TextIOBase,
    smiles_allowed: Set[str],
) -> Tuple[bytes, FilterReport, Set[str]]:
    """
    Stream-parse an MGF and return the filtered MGF bytes + report + matched_smiles_set.
    This avoids storing all blocks in memory.

    Rules:
    - Extract SMILES= within each BEGIN/END IONS block (first occurrence).
    - Keep block only if SMILES exactly matches a value in smiles_allowed.
    """
    out = io.StringIO()
    in_block = False
    current_lines: List[str] = []
    current_smiles: Optional[str] = None

    total = 0
    matched = 0
    no_smiles = 0
    matched_smiles_set: Set[str] = set()

    for raw in mgf_text_stream:
        line = raw.rstrip("\n")
        stripped = line.strip()

        if stripped == "BEGIN IONS":
            in_block = True
            current_lines = [line + "\n"]
            current_smiles = None
            continue

        if stripped == "END IONS" and in_block:
            current_lines.append(line + "\n")
            total += 1

            if not current_smiles:
                no_smiles += 1

            if current_smiles and current_smiles in smiles_allowed:
                out.writelines(current_lines)
                if not current_lines[-1].endswith("\n"):
                    out.write("\n")
                matched += 1
                matched_smiles_set.add(current_smiles)

            in_block = False
            current_lines = []
            current_smiles = None
            continue

        if in_block:
            current_lines.append(line + "\n")
            if current_smiles is None and stripped.startswith("SMILES="):
                current_smiles = stripped.split("=", 1)[1].strip()

    report = FilterReport(
        total_blocks=total,
        matched_blocks=matched,
        blocks_with_no_smiles=no_smiles,
        unique_smiles_in_table=len(smiles_allowed),
        unique_smiles_matched=len(matched_smiles_set),
    )

    return out.getvalue().encode("utf-8"), report, matched_smiles_set


# -----------------------------
# Streamlit UI
# -----------------------------

st.set_page_config(page_title="MGF Filter by SMILES", layout="wide")
st.title("MGF Filter by SMILES (Stream-Filter GNPS Library Blocks)")

st.markdown(
    """
Upload:
1) An **MGF** file containing blocks with `BEGIN IONS ... END IONS` and (ideally) a `SMILES=` line inside each block.
2) A **CSV/TSV** table containing a **SMILES** column (optionally `compound_name`).

The app will keep only spectra whose `SMILES=` exactly matches one of the SMILES in your table.
"""
)

with st.sidebar:
    st.header("Inputs")
    mgf_file = st.file_uploader("Upload MGF", type=["mgf", "txt"], accept_multiple_files=False)
    table_file = st.file_uploader(
        "Upload table (CSV/TSV) with SMILES column",
        type=["csv", "tsv", "txt"],
        accept_multiple_files=False,
    )

    st.divider()
    st.header("Options")
    st.checkbox("Strict SMILES match (exact string)", value=True, disabled=True)
    run_btn = st.button("Run filtering", type="primary", use_container_width=True)

if not mgf_file or not table_file:
    st.info("Upload both files in the sidebar to begin.")
    st.stop()

if not run_btn:
    st.stop()

# Load SMILES set
try:
    table_bytes = table_file.getvalue()
    smiles_allowed, detected_sep, smiles_col = load_smiles_set(table_bytes)
except Exception as e:
    st.error(f"Failed to read the table: {e}")
    st.stop()

st.success(
    f"Loaded {len(smiles_allowed):,} unique SMILES from column **{smiles_col}** "
    f"(delimiter detected: `{repr(detected_sep)}`)."
)

# Stream-filter MGF
mgf_bytes = mgf_file.getvalue()

with st.spinner("Filtering MGF blocks (streaming)..."):
    text_stream = io.StringIO(mgf_bytes.decode("utf-8", errors="ignore"))
    filtered_bytes, report, matched_smiles_set = iter_filtered_mgf_bytes(text_stream, smiles_allowed)

# Results
col1, col2, col3, col4, col5 = st.columns(5)
col1.metric("Total spectra (blocks)", f"{report.total_blocks:,}")
col2.metric("Matched spectra", f"{report.matched_blocks:,}")
col3.metric("Blocks without SMILES", f"{report.blocks_with_no_smiles:,}")
col4.metric("Unique SMILES (table)", f"{report.unique_smiles_in_table:,}")
col5.metric("Unique SMILES matched", f"{report.unique_smiles_matched:,}")

st.divider()

# Warnings / notes
if report.total_blocks == 0:
    st.warning("No MGF blocks were detected (no 'BEGIN IONS' / 'END IONS' pairs found). Check the input file format.")
elif report.matched_blocks == 0:
    st.warning(
        "No blocks matched the SMILES list. Common causes:\n"
        "- The MGF blocks do not contain `SMILES=` lines\n"
        "- SMILES strings differ in formatting (salts, stereochemistry, canonicalization)\n"
        "- The table SMILES column is not what you expected\n"
    )
elif report.blocks_with_no_smiles > 0:
    st.info(
        "Some blocks had no `SMILES=` line and were automatically excluded. "
        "If this is unexpected, confirm your MGF source includes SMILES annotations."
    )

# Download
output_name = f"filtered_{mgf_file.name.rsplit('.', 1)[0]}.mgf"
st.download_button(
    label="Download filtered MGF",
    data=filtered_bytes,
    file_name=output_name,
    mime="chemical/x-mgf",
    use_container_width=True,
)

# Preview (first N lines)
with st.expander("Preview filtered MGF (first 200 lines)"):
    preview = filtered_bytes.decode("utf-8", errors="ignore").splitlines()[:200]
    st.code("\n".join(preview) if preview else "(empty)", language="text")

with st.expander("Preview SMILES loaded from table (first 50)"):
    smi_preview = sorted(smiles_allowed)[:50]
    st.code("\n".join(smi_preview), language="text")


# -----------------------------
# RDKit grid
# -----------------------------

st.divider()
st.subheader("SMILES structure matrix (RDKit)")

if not RDKit_AVAILABLE:
    st.warning(
        "RDKit is not available in this environment. "
        "Install it (e.g., conda install -c conda-forge rdkit) to enable structure plotting."
    )
else:
    try:
        df_table, _sep2 = load_table_df(table_bytes)
        smiles_col2 = detect_smiles_column(df_table.columns)
        if not smiles_col2:
            st.warning("Could not detect a SMILES column for structure plotting.")
            st.stop()

        name_col = detect_compound_name_column(df_table.columns)  # may be None

        only_matched = st.checkbox("Plot only SMILES that matched spectra in the filtered MGF", value=True)
        max_mols = st.slider("Max molecules to display", min_value=12, max_value=240, value=60, step=12)
        mols_per_row = st.slider("Molecules per row", min_value=3, max_value=10, value=5, step=1)

        df_plot = df_table.copy()
        df_plot[smiles_col2] = df_plot[smiles_col2].astype(str).map(lambda s: s.strip())

        if only_matched:
            df_plot = df_plot[df_plot[smiles_col2].isin(matched_smiles_set)]

        df_plot = df_plot[df_plot[smiles_col2].astype(bool)]
        df_plot = df_plot.drop_duplicates(subset=[smiles_col2], keep="first")
        df_plot = df_plot.head(max_mols)

        if df_plot.empty:
            st.info("No SMILES available to plot under the current selection.")
        else:
            smiles_list = df_plot[smiles_col2].tolist()

            if name_col and name_col in df_plot.columns:
                legends = df_plot[name_col].fillna("").astype(str).tolist()
                legends = [leg.strip() if leg.strip() else smi for leg, smi in zip(legends, smiles_list)]
            else:
                legends = smiles_list

            # truncate legends a bit (keeps the grid readable)
            legends = [s[:30] + "…" if len(s) > 31 else s for s in legends]

            mols: List = []
            mol_legends: List[str] = []
            invalid = 0

            for smi, leg in zip(smiles_list, legends):
                m = Chem.MolFromSmiles(smi)
                if m is None:
                    invalid += 1
                    continue
                mols.append(m)
                mol_legends.append(leg)

            if not mols:
                st.warning("All SMILES in the selected set failed RDKit parsing.")
            else:
                if invalid > 0:
                    st.info(f"Skipped {invalid} invalid SMILES (RDKit could not parse).")

                img = Draw.MolsToGridImage(
                    mols,
                    molsPerRow=mols_per_row,
                    subImgSize=(260, 220),
                    legends=mol_legends,
                    useSVG=False,
                )

                if isinstance(img, Image.Image):
                    st.image(img, caption="Structure grid (SMILES → 2D RDKit depiction)", use_container_width=True)
                else:
                    st.write(img)

            with st.expander("Table used for plotting (preview)"):
                show_cols = [c for c in [name_col, smiles_col2] if c]
                st.dataframe(df_plot[show_cols].head(50))

    except Exception as e:
        st.error(f"Failed to plot structures: {e}")
