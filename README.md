# 🧬 Oligo-Barcode-Profiler

**Oligo-Barcode-Profiler** is a Python-based tool designed to analyze barcoded oligonucleotide libraries from paired-end Illumina sequencing data. In such libraries, each oligonucleotide is tagged with a short, random barcode - typically 8–12 bases long - that uniquely identifies the molecule or the plasmid it originated from. These barcodes are used to track individual molecules, quantify amplification biases, or group sequencing reads from the same original DNA source.

Following each barcode is a short sequence context, often designed to test the effect of local sequence environment on specific biochemical events - such as nucleotide modifications or repair. This tool focuses on quantifying nucleotide edits, such as **U→C correction events**, by extracting the barcode and a surrounding 3-base context from millions of sequencing reads and comparing forward and reverse read pairs to validate true edits.

This project was developed to support research involving **sequence context-dependent mutation profiling**, such as those used in error-correction modeling in oligo libraries.

---

## 💡 What does this project do?

- Parses compressed or uncompressed paired-end FASTQ files.
- Extracts 9 bp random barcodes followed by a 3 bp sequence (where the middle base is the expected base, e.g., U).
- Compares forward and reverse reads to validate base calls.
- Identifies base conversions (and optionally other mutations).
- Groups reads by barcode and aggregates per-context mutation frequencies.
- Outputs:
  - Mutation summary tables as CSV files
  - Plots (mutation rates, barcode abundance, context bias)
  - Optional interactive HTML report for exploring results visually

---

## 📥 Input

- Paired-end FASTQ files (`sample_R1.fastq(.gz)`, `sample_R2.fastq(.gz)`)
- Parameters:
  - Known anchor sequences (anchor1, anchor2, anchor3)
  - Allowed number of mismatches per anchor (--anch1-mm, --anch2-mm, --anch3-mm, default = 1 if not specified)
  - Barcode length (e.g., 9 bp)
  - Context length (e.g., 3 bp)
  - Expected base (e.g., `U`) at position 2 of the 3-bp context
- Note: This tool expects input FASTQ files to be pre-processed, including quality control (QC), adapter trimming, and read filtering if needed. Both files in each pair (R1 and R2) must contain matching, synchronized reads.

---

## 📤 Output

- `mutation_summary.csv`: Per-barcode mutation stats
- `context_profile.csv`: Aggregated mutation frequency by 3-bp context
- `mutation_plots/`: PDF or PNG plots of mutation patterns
- `report.html`: Optional interactive dashboard with Plotly or Dash

---

## ⚙️ Technical details

### How to install

You will need Python 3.8+ and a few bioinformatics/data science packages:

```bash
git clone https://github.com/ElizabethPopov/Oligo-Barcode-Profiler.git
cd Oligo-Barcode-Profiler
pip install -r requirements.txt
```

**Dependencies:**
- Biopython
- pandas
- matplotlib / seaborn
- plotly (optional)
- gzip (for fastq.gz handling)

### How to run

```bash
python analyze_barcodes.py \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --anchor1 AGCTTG \
  --anchor2 TAG \
  --anchor3 CUTGGTC \
  --anch1-mm 2 \
  --anch2-mm 1 \
  --anch3-mm 2 \
  --expected-base U \
  --barcode-length 9 \
  --context-length 3 \
  --output-dir results/
```

**Note:**
- The program will search for reads containing the pattern:  
  `[anchor1][barcode][anchor2][context][anchor3]`
- For example, with `--anchor1 AGCTTG`, `--anchor2 TAG` and `--anchor3 CUTGGTC`, the tool expects:
  - A 9-bp barcode immediately after `AGCTTG`
  - Followed by the anchor sequence `TAG`
  - Then a 3-bp context (e.g., AUG)
  - Finally, the anchor sequence `CUTGGTC` appears immediately after the context
  - In this example, anchor1 and anchor3 allow up to 2 mismatches, while anchor2 allows only 1 mismatch. These thresholds can be adjusted to match your data quality and tolerance.
- The middle base of the context (position 2) is compared to the expected base (`U`) to detect conversions (e.g., U→C).

---

## ✅ Example use case

You have 24 plasmid samples containing 9-bp barcoded oligos with a U at the center of a 3-bp context. This tool helps you determine how frequently that U was corrected to C, and whether that frequency depends on the surrounding sequence.

---

## 🧪 Tests

```bash
pytest tests/
```

Test suite includes:
- Synthetic read parsing
- Barcode and context extraction
- Mutation frequency calculations

---

## 📚 Attribution

This project was developed as part of the [Python Programming course](https://github.com/Code-Maven/wis-python-course-2025-03) at [Weizmann School of Science](https://www.weizmann.ac.il/pages/).
