# 🧬 Oligo-Barcode-Profiler

**Oligo-Barcode-Profiler** is a Python-based tool designed to analyze barcoded oligonucleotide libraries from paired-end Illumina sequencing data. In such libraries, each oligonucleotide is tagged with a short, random barcode - typically 8–12 bases long - that uniquely identifies the molecule or the plasmid it originated from. These barcodes are used to track individual molecules, quantify amplification biases, or group sequencing reads from the same original DNA source.

Following each barcode is a short sequence context, often designed to test the effect of local sequence environment on specific biochemical events - such as nucleotide modifications or repair. This tool focuses on quantifying nucleotide edits, such as **U→C correction events**, by extracting the barcode and a surrounding 3-base context from millions of sequencing reads and comparing forward and reverse read pairs to validate true edits.

This project was developed to support research involving **sequence context-dependent mutation profiling**, such as those used in error-correction modeling in oligo libraries.

---

## 💡 What does this project do?

- Parses uncompressed paired-end FASTQ files.
- Extracts 9 bp random barcodes followed by a 3 bp sequence (where the middle base is the original base of interest, e.g., T for U→C correction).
- Compares forward and reverse reads to validate base calls.
- Identifies base corrections (e.g., U→C) based on sequence context.
- Groups reads by barcode and aggregates per-context mutation frequencies.
- Outputs:
  - Summary tables as CSV files (barcode and context counts before and after filtering, context distribution, correction rates)
  - Plots (barcode count distribution, correction vs. no-correction rates)

---

## 📥 Input

- Paired-end FASTQ files (`sample1_R1.fastq(.fq)`, `sample1_R2.fastq(.fq)`)
- Parameters:
  - Known anchor sequences (`--anchor1`, `--anchor2`, `--anchor3`)
  - Original context before the correction event (`--context`)
  - Corrected context (`--corrected-context`)
  - Allowed number of mismatches per anchor  
    (`--anch1-mm`, `--anch2-mm`, `--anch3-mm`; defaults: 2, 1, 2)
  - Minimum percent for context validation  
    (`--min-pct`; default: 40%)
  - Barcode length (`--barcode-length`; default: 9 bp)
  - Context length (`--context-length`; default: 3 bp)
  - Output directory for saving results (`--output-dir`; default: output)
- **Note:** This tool expects input FASTQ files to be pre-processed (e.g., adapter trimming, QC, filtering).  
  R1 and R2 files must contain **synchronized read pairs**.

---

## 📤 Output

Output files are saved to the directory specified by `--output-dir`, and include:
- `<sample>_barcode_count_summary_unfiltered.csv`: Raw counts per barcode after paired-read validation
- `<sample>_barcode_count_distribution.png`: Plot of barcode count distribution
- `<sample>_barcode_count_summary_validated.csv`: Barcode counts after filtering for high-confidence context assignments (based on --min-pct)
- `<sample>_context_distribution_percent.csv`: Aggregated context distribution
- `<sample>_correction_summary.csv`: Context correction statistics
- `<sample>_correction_summary.png`: Barplot of corrected vs. uncorrected contexts

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
- regex
- matplotlib
- seaborn
- pytest (for tests only)

### How to run

```bash
python analyze_barcodes.py \
  --r1 sample1_R1.fastq \
  --r2 sample1_R2.fastq \
  --min-pct 40 \
  --context CTA \
  --corrected-context CCA \
  --anchor1 CGTAC \
  --anchor2 TTCGA \
  --anchor3 GGACATT \
  --anch1-mm 2 \            # default
  --anch2-mm 1 \            # default
  --anch3-mm 2 \            # default
  --barcode-length 9 \      # default
  --context-length 3 \      # default
  --output-dir output/      # default
```

**Note:**
- Only uncompressed `.fastq` or `.fq` files are supported currently.
- The program will search for reads containing the pattern:  
  `[anchor1][barcode][anchor2][context][anchor3]`
- For example, with `--anchor1 CGTAC`, `--anchor2 TTCGA` and `--anchor3 GGACATT`, the tool expects:
  - A 9-bp barcode immediately after `CGTAC`
  - Followed by the anchor sequence `TTCGA`
  - Then a 3-bp context (e.g., CTA)
  - Finally, the anchor sequence `GGACATT` appears immediately after the context
  - In this example, anchor1 and anchor3 allow up to 2 mismatches, while anchor2 allows only 1 mismatch. These thresholds can be adjusted to match your data quality and tolerance.
- The original context (uncorrected) and the corrected context (`CCA`) are calculated to detect correction rates (in this example, U→C correction).

### Average runtime

On a typical modern laptop (e.g., Apple M1 Pro), analyzing one sample with ~3 million read pairs takes approximately **77–100 seconds**.  
Actual runtime may vary depending on read file size, the complexity of the search pattern, disk speed, and CPU performance.

---

## ✅ Example use case

You have 24 plasmid samples containing 9-bp barcoded oligos with a U at the center of a 3-bp context. This tool helps you determine how frequently that U was corrected to C by searching for contexts that contain T, and whether that frequency depends on the surrounding sequence.

---

## 🚀 Quick Start: Try It Out with Example Data

This repository includes a small synthetic dataset (`example_R1.fastq`, `example_R2.fastq`) to help you test the tool immediately without needing large FASTQ files.

### 🔧 Step-by-step:

```bash
# Clone the repository
git clone https://github.com/ElizabethPopov/Oligo-Barcode-Profiler.git
cd Oligo-Barcode-Profiler

# Install required packages
pip install -r requirements.txt

# Run the pipeline on example data
python analyze_barcodes.py \
  --r1 data/example_R1.fastq \
  --r2 data/example_R2.fastq \
  --min-pct 40 \
  --context CTA \
  --corrected-context CCA \
  --anchor1 CGTAC \
  --anchor2 TTCGA \
  --anchor3 GGACATT \
  --barcode-length 9 \
  --context-length 3 \
  --output-dir test_outputs/ \
  --verbose
```

Output files will appear in test_outputs/, including:
- `example_barcode_count_summary_unfiltered.csv`
- `example_barcode_count_summary_validated.csv`
- `example_context_distribution_percent.csv`
- `example_correction_summary.csv`
- `example_barcode_count_distribution.png`
- `example_correction_summary.png`

---

## 🧪 Tests

Run all tests with:
```bash
pytest tests/
```

Test suite includes:
- Validation of barcode and context extraction from R1 and R2 reads
- Read parsing and FASTQ structural checks
- Anchor mismatch handling
- Summary statistics and correction events calculations
- Exception handling (e.g., malformed reads, missing files, invalid input parameters)
- Visual output generation (barcode and correction summary plots tested for successful creation)

---

## 📊 Output from Real Dataset
To demonstrate real-world results, the `real_data_outputs/` folder includes output files generated from a sample with ~3 million read pairs:
- `real_data_outputs/S1_amp_S1_context_distribution_percent.csv`
- `real_data_outputs/S1_amp_S1_correction_summary.csv`
- `real_data_outputs/S1_amp_S1_correction_summary.png`
- `real_data_outputs/S1_amp_S1_barcode_count_distribution.png`
These reflect actual correction and context distribution statistics from an experimental oligo library. You can use them to:
- Compare results against your own data
- Explore the output format and file structure

**Note:**
I did **not** include the following files due to their large size:
- `example_outputs/S1_amp_S1_barcode_count_summary_unfiltered.csv`
- `example_outputs/S1_amp_S1_barcode_count_summary_validated.csv`

The included files should still be sufficient to give you a solid understanding of the output. 

---

## 🧪 Example Output from Sample Data

The `example_outputs/` folder contains output files generated by running the tool on a small test dataset (`data/`) included in this repository.

Included example files:
- `example_barcode_count_summary_unfiltered.csv`
- `example_barcode_count_summary_validated.csv`
- `example_context_distribution_percent.csv`
- `example_correction_summary.csv`
- `example_barcode_count_distribution.png`
- `example_correction_summary.png`

These outputs were generated from a lightweight synthetic dataset designed for testing.

You can reproduce them by running:
```bash
python analyze_barcodes.py \
  --r1 data/example_R1.fastq \
  --r2 data/example_R2.fastq \
  --min-pct 40 \
  --context CTA \
  --corrected-context CCA \
  --anchor1 CGTAC \
  --anchor2 TTCGA \
  --anchor3 GGACATT \
  --output-dir example_outputs/
```

---

## 📚 Attribution

This project was developed as part of the [Python Programming course](https://github.com/Code-Maven/wis-python-course-2025-03) at [Weizmann School of Science](https://www.weizmann.ac.il/pages/).
