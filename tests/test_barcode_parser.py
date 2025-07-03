import os
import sys
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from barcode_parser import compile_pattern, extract_barcode_and_context, reverse_complement, extract_barcode_and_context_R2


def test_compile_pattern_with_mismatches():
    pattern = compile_pattern("AAA", "TTT", "GGG", 5, 3, 1, 1, 1)
    seq = "AAGCCCGTATTACAGCG"
    match = pattern.search(seq)
    assert match is not None
    assert match.group('barcode') == "CCCGT"
    assert match.group('context') == "ACA"


def test_extract_barcode_and_context():
    pattern = compile_pattern("AAA", "TTT", "GGG", 5, 3, 0, 0, 0)
    read = "AAACCCGTTTTACAGGG"
    result = extract_barcode_and_context(read, pattern)
    assert result is not None
    assert result['barcode'] == "CCCGT"
    assert result['context'] == "ACA"
    assert result['full_sequence'] == "AAACCCGTTTTACAGGG"


def test_extract_barcode_and_context_R2():
    from Bio.Seq import Seq
    pattern = compile_pattern("AAA", "TTT", "GGG", 5, 3, 0, 0, 0)
    read = "CCCTGGAAAGTCGGTTTGTACG"     # Reverse complement of "CCCGTAAACAGTTTGGG"
    result = extract_barcode_and_context_R2(read, pattern)
    assert result is not None
    assert result['barcode'] == "CCGAC"
    assert result['context'] == "CCA"


def test_compile_pattern_invalid_anchor():
    with pytest.raises(ValueError, match="anchor1 contains invalid characters"):
        compile_pattern("AAX", "TTT", "GGG", 5, 3, 0, 0, 0)


def test_compile_pattern_negative_mismatch():
    with pytest.raises(ValueError, match="mm1 must be a non-negative integer"):
        compile_pattern("AAA", "TTT", "GGG", 5, 3, -1, 0, 0)