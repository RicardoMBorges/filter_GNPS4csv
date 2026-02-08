
```markdown
# GNPS MGF Filter & Structure Viewer (Streamlit App)

This Streamlit application filters **GNPS-style MGF libraries** using a list of **SMILES** provided in a CSV/TSV table and optionally visualizes the corresponding **chemical structures as a matrix/grid**.

It is designed for:
- GNPS library curation
- Class-level or priority-list extraction
- Teaching / exploratory metabolomics
- Lightweight deployment (no RDKit required)

---

## ✨ What this app does

### 1. Filters an MGF file by SMILES
- Reads an MGF file block-by-block (`BEGIN IONS … END IONS`)
- Extracts `SMILES=` annotations inside each block
- Keeps only spectra whose SMILES **exactly match** the SMILES in your table
- Outputs a **new filtered MGF** for download

### 2. Reports useful statistics
- Total spectra in the original MGF
- Number of matched spectra
- Number of spectra without SMILES
- Number of unique SMILES provided
- Number of unique SMILES actually found in the MGF

### 3. Visualizes chemical structures (no RDKit)
- Builds a **matrix of structures** from SMILES
- Uses an online structure resolver (CACTUS, NCI)
- Displays **compound names** if your table has a `compound_name` column
- Allows filtering to only SMILES that matched spectra

---

## 📁 Required input files

### 1️⃣ MGF file
Your MGF **must**:
- Use standard `BEGIN IONS` / `END IONS` blocks
- Contain a line like:
```

SMILES=CC1=CC=CC=C1

````
inside each block (GNPS libraries usually do)

### 2️⃣ CSV / TSV table with SMILES
Minimum requirement:
- A column named one of:
- `SMILES`
- `compound_smiles`
- `structure_smiles`

Optional (recommended):
- A column named:
- `compound_name`

Example:

```csv
compound_name,SMILES
Quercetin,OC1=CC=C2C(=O)C(O)=C(O)C=C2C1
Kaempferol,OC1=CC=C2C(=O)C(O)=C(O)C=C2C1O
````

---

## 🚀 How to run the app

### 1️⃣ Create a Python environment

```bash
python -m venv mgf_env
source mgf_env/bin/activate  # Linux/macOS
mgf_env\Scripts\activate     # Windows
```

### 2️⃣ Install dependencies

```bash
pip install streamlit pandas pillow
```

> ⚠️ RDKit is **not required**.

### 3️⃣ Run Streamlit

```bash
streamlit run app.py
```

Your browser will open automatically.

---

## 🧪 Step-by-step tutorial (inside the app)

1. **Upload the MGF file**

   * Example: `GNPS-LIBRARY.mgf`

2. **Upload the CSV/TSV priority list**

   * Example: `GNPS_Shikimates_Phenylpropanoids_Flavonoids.csv`

3. Click **“Run filtering”**

4. Review:

   * Metrics (how many spectra matched)
   * Warnings (missing SMILES, no matches, etc.)

5. **Download the filtered MGF**

   * File name: `filtered_<original>.mgf`

6. Scroll down to **SMILES structure matrix**

   * Adjust:

     * number of molecules
     * grid size
     * matched-only toggle
   * Inspect chemical diversity visually

---

## 🧬 About structure visualization (important)

This app **does not use RDKit**.

Instead:

* It sends each SMILES to the **CACTUS Structure Resolver**
* Receives a rendered PNG image
* Displays it in a grid

### Requirements:

* 🌐 Internet access from the Streamlit environment

### Limitations:

* If a SMILES is invalid → image not shown
* If the resolver is down → image not shown
* This is **visual inspection only**, not cheminformatics analysis

---

## 🔬 Scientific design choices

* **Exact SMILES matching**

  * Avoids false positives
  * Reproducible and transparent
* **Stream-based MGF parsing**

  * Handles very large GNPS libraries
  * Low memory footprint
* **No dependency on cheminformatics toolkits**

  * Easy deployment (Streamlit Cloud, Docker, shared servers)

---

## 🧩 Typical use cases

* Extract flavonoids, alkaloids, terpenes from GNPS libraries
* Build class-specific MGF subsets
* Teaching fragmentation + structure relationships
* Preparing curated libraries for:

  * Molecular networking
  * MassQL benchmarking
  * Class-level annotation studies

---

## 🛠 Troubleshooting

### ❌ “No blocks matched”

Possible causes:

* MGF does not contain `SMILES=`
* SMILES formatting differs (salts, stereochemistry)
* Wrong column detected in the CSV

### ❌ “RDKit not available”

This is expected.
The app **intentionally does not require RDKit**.

### ❌ Structure images missing

* Internet blocked
* Invalid SMILES
* Temporary resolver outage

---

## 📜 License & reuse

This code is intended for:

* Academic research
* Teaching
* Method development

Feel free to adapt, fork, and extend.

If you use it in a publication or course, attribution is appreciated.

---

## 🙌 Credits

Designed for GNPS-based metabolomics workflows
with a focus on **chemical reasoning, transparency, and pedagogy**.



