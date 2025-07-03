import os
import sys
import pytest
import regex as re
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fastq_reader import process_paired_reads, load_validated_reads
from barcode_parser import compile_pattern


def test_process_paired_reads(tmp_path):
    # Create test FASTQ files
    r1_path = tmp_path / "S1_amp_S1_R1.fastq"
    r2_path = tmp_path / "S1_amp_S1_R2.fastq"

    # Matching read with anchors/barcode/context
    with open(r1_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/1\nCGTACGCATCTAAATTCGACCAGGACATT\n+\n" + "I" * 29 + "\n")

    with open(r2_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/2\nAATGTCCTGGTCGAATTTAGATGCGTACG\n+\n" + "I" * 29 + "\n")

    # Compile regex pattern
    pattern = compile_pattern("CGTAC", "TTCGA", "GGACATT", 9, 3, 2, 1, 2)

    # Call function
    matched = process_paired_reads(str(r1_path), str(r2_path), pattern)

    # Assertions
    assert isinstance(matched, list)
    assert len(matched) == 1
    assert matched[0]['barcode'] == "GCATCTAAA"
    assert matched[0]['context'] == "CCA"


def test_load_validated_reads(tmp_path):
    r1_path = tmp_path / "S1_amp_S1_R1.fastq"
    r2_path = tmp_path / "S1_amp_S1_R2.fastq"

    with open(r1_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/1\nCGTACGCATCTAAATTCGACCAGGACATT\n+\n" + "I" * 29 + "\n")

    with open(r2_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/2\nAATGTCCTGGTCGAATTTAGATGCGTACG\n+\n" + "I" * 29 + "\n")

    pattern = compile_pattern("CGTAC", "TTCGA", "GGACATT", 9, 3, 2, 1, 2)

    sample_name, reads = load_validated_reads(str(r1_path), str(r2_path), pattern)

    assert sample_name == "S1_amp_S1"
    assert isinstance(reads, list)
    assert reads and 'barcode' in reads[0]


def test_mismatched_headers(tmp_path):
    r1_path = tmp_path / "S1_amp_S1_R1.fastq"
    r2_path = tmp_path / "S1_amp_S1_R2.fastq"

    with open(r1_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/1\nCGTACGCATCTAAATTCGACCAGGACATT\n+\n" + "I" * 29 + "\n")

    with open(r2_path, 'w') as f:
        f.write("@F350043104L1C001R00100001973/2\nAATGTCCTGGTCGAATTTAGATGCGTACG\n+\n" + "I" * 29 + "\n") # Different last digit in the header

    pattern = compile_pattern("CGTAC", "TTCGA", "GGACATT", 9, 3, 2, 1, 2)
    matched = process_paired_reads(str(r1_path), str(r2_path), pattern)

    assert matched == []


def test_incomplete_fastq(tmp_path):
    r1_path = tmp_path / "S1_amp_S1_R1.fastq"
    r2_path = tmp_path / "S1_amp_S1_R2.fastq"

    with open(r1_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/1\nCGTACGCATCTAAATTCGACCAGGACATT\n+\n")  # No quality line

    with open(r2_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/2\nAATGTCCTGGTCGAATTTAGATGCGTACG\n+\n" + "I" * 29 + "\n")

    pattern = compile_pattern("CGTAC", "TTCGA", "GGACATT", 9, 3, 2, 1, 2)

    with pytest.raises(ValueError, match=f"Malformed FastQ entry detected in file '{r1_path}' or '{r2_path}'."):
        process_paired_reads(str(r1_path), str(r2_path), pattern)


def test_no_matching_pattern(tmp_path):
    r1_path = tmp_path / "S1_amp_S1_R1.fastq"
    r2_path = tmp_path / "S1_amp_S1_R2.fastq"

    with open(r1_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/1\nAAACGACGGCCAGTGAGCGCGCGTAATACGACTCACTATAGGGCGAATTGGGTACCCGGTCCCCCCTCGAGGTCGACCGTATCGATAAGCTTGATATCGAATTCTCATTGCTACAGCCGTGGAGCTCCAGCTTTTGTTCCCTTTAGTGAG\n+\n" + "I" * 150 + "\n")
    with open(r2_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/2\nGCCAAGCGCGCAATTAACCCTCACTAAAGGGAACAAAAGCTGGAGCTCCACCGCTGTAGCAATGAGAATTCGATATCAAGCTTATCGATACCGTCGACCTCCAGGGAGGGCCCGGTACCCAATTCGCCCTATACTGAGTCGTATTACGCG\n+\n" + "I" * 150 + "\n")

    pattern = compile_pattern("CGTAC", "TTCGA", "GGACATT", 9, 3, 2, 1, 2)
    matched = process_paired_reads(str(r1_path), str(r2_path), pattern)

    assert matched == []


def test_invalid_sample_name(tmp_path):
    r1_path = tmp_path / "sample_r1.fastq"
    r2_path = tmp_path / "sample_r2.fastq"

    with open(r1_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/1\nCGTACGCATCTAAATTCGACCAGGACATT\n+\n" + "I" * 29 + "\n")
    with open(r2_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/2\nAATGTCCTGGTCGAATTTAGATGCGTACG\n+\n" + "I" * 29 + "\n")

    pattern = compile_pattern("CGTAC", "TTCGA", "GGACATT", 9, 3, 2, 1, 2)

    with pytest.raises(ValueError, match="Could not infer sample name"):
        load_validated_reads(str(r1_path), str(r2_path), pattern)


def test_missing_files():
    pattern = compile_pattern("AAA", "TTT", "GGG", 5, 3, 0, 0, 0)
    with pytest.raises(FileNotFoundError, match="Neither R1 nor R2 files found"):
        process_paired_reads("nonexistent_R1.fastq", "nonexistent_R2.fastq", pattern)


def test_missing_r1_file(tmp_path):
    pattern = compile_pattern("AAA", "TTT", "GGG", 5, 3, 0, 0, 0)
    r2 = tmp_path / "sample_R2.fastq"
    r2.write_text("@header\nACGT\n+\n!!!!\n")
    with pytest.raises(FileNotFoundError, match="R1 file not found"):
        process_paired_reads("nonexistent_R1.fastq", str(r2), pattern)


def test_missing_r2_file(tmp_path):
    pattern = compile_pattern("AAA", "TTT", "GGG", 5, 3, 0, 0, 0)
    r1 = tmp_path / "sample_R1.fastq"
    r1.write_text("@header\nACGT\n+\n!!!!\n")
    with pytest.raises(FileNotFoundError, match="R2 file not found"):
        process_paired_reads(str(r1), "nonexistent_R2.fastq", pattern)


def test_sample_name_mismatch(tmp_path):
    r1_path = tmp_path / "sample1_R1.fastq"
    r2_path = tmp_path / "sample2_R2.fastq"

    with open(r1_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/1\nCGTACGCATCTAAATTCGACCAGGACATT\n+\n" + "I" * 29 + "\n")
    with open(r2_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/2\nAATGTCCTGGTCGAATTTAGATGCGTACG\n+\n" + "I" * 29 + "\n")

    pattern = compile_pattern("CGTAC", "TTCGA", "GGACATT", 9, 3, 2, 1, 2)

    with pytest.raises(ValueError, match="R1 and R2 filenames appear to come from different samples"):
        load_validated_reads(str(r1_path), str(r2_path), pattern)


def test_invalid_fastq_extension(tmp_path):
    r1_path = tmp_path / "sample_R1.txt"
    r2_path = tmp_path / "sample_R2.fastq"

    with open(r1_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/1\nCGTACGCATCTAAATTCGACCAGGACATT\n+\n" + "I" * 29 + "\n")
    with open(r2_path, 'w') as f:
        f.write("@F350043104L1C001R00100001972/2\nAATGTCCTGGTCGAATTTAGATGCGTACG\n+\n" + "I" * 29 + "\n")

    pattern = compile_pattern("CGTAC", "TTCGA", "GGACATT", 9, 3, 2, 1, 2)

    with pytest.raises(ValueError, match="Expected FastQ files"):
        load_validated_reads(str(r1_path), str(r2_path), pattern)