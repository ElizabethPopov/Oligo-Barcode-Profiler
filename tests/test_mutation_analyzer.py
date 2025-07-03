import os
import sys
import pandas as pd
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from mutation_analyzer import get_validated_contexts, compute_correction_stats, get_context_distribution

def test_get_validated_contexts():
    context_counts = {'CTA': 7, 'CCA': 3}
    result = get_validated_contexts(context_counts, min_pct=60)
    assert len(result) == 1
    assert result[0]['context'] == 'CTA'
    assert result[0]['count'] == 7
    assert result[0]['percent'] == 100.0


def test_get_validated_contexts_multiple_pass():
    context_counts = {'CTA': 6, 'CCA': 4}
    result = get_validated_contexts(context_counts, min_pct=30)
    assert len(result) == 2
    assert {r['context'] for r in result} == {'CTA', 'CCA'}


def test_compute_correction_stats():
    df = pd.DataFrame([
        {'sample': 'Sample1', 'barcode': 'AAAAAAAAG', 'contexts': ['CTA', 'CCA'], 'context_counts': {'CTA': 7, 'CCA': 3}, 'total_count': 10, 'R1_sequences': set(), 'R2_sequences': set(), 'validated_contexts': [{'context': 'CTA', 'count': 7, 'percent': 70.0}]},
        {'sample': 'Sample1', 'barcode': 'AAAAAAAAG', 'contexts': ['CCA', 'CTA'], 'context_counts': {'CTA': 7, 'CCA': 3}, 'total_count': 10, 'R1_sequences': set(), 'R2_sequences': set(), 'validated_contexts': [{'context': 'CCA', 'count': 3, 'percent': 30.0}]},
        {'sample': 'Sample1', 'barcode': 'TTTTTTTTT', 'contexts': ['CTA', 'CCA'], 'context_counts': {'CTA': 7, 'CCA': 3}, 'total_count': 10, 'R1_sequences': set(), 'R2_sequences': set(), 'validated_contexts': [{'context': 'CTA', 'count': 7, 'percent': 70.0}]},
    ])

    result = compute_correction_stats(df, 'CTA', 'CCA', threshold=40)
    assert result['correction_events'] == 1
    assert result['no_correction_events'] == 2
    assert result['correction_%'] == 33.33
    assert result['no_correction_%'] == 66.67
    assert result['total_barcodes'] == 3


def test_compute_correction_stats_no_validated():
    df = pd.DataFrame([
        {'sample': 'Sample1', 'barcode': 'AAAAAAAAG', 'contexts': ['CCA', 'CGA', 'GCA'], 'context_counts': {'CCA': 3, 'CGA': 3, 'GCA': 3}, 'total_count': 9, 'R1_sequences': set(), 'R2_sequences': set(), 'validated_contexts': []},
    ])
    result = compute_correction_stats(df, 'CTA', 'CCA', threshold=40)
    assert result['correction_events'] == 0
    assert result['no_correction_events'] == 0
    assert result['total_barcodes'] == 0
    assert result['correction_%'] == 0
    assert result['no_correction_%'] == 0


def test_get_context_distribution():
    df = pd.DataFrame([
        {'sample': 'Sample1', 'barcode': 'AAAAAAAAG', 'contexts': ['CTA', 'CCA'], 'context_counts': {'CTA': 7, 'CCA': 3}, 'total_count': 10, 'R1_sequences': set(), 'R2_sequences': set(), 'validated_contexts': [{'context': 'CTA', 'count': 7, 'percent': 70.0}]},
        {'sample': 'Sample1', 'barcode': 'AAAAAAAAG', 'contexts': ['CCA', 'CTA'], 'context_counts': {'CTA': 7, 'CCA': 3}, 'total_count': 10, 'R1_sequences': set(), 'R2_sequences': set(), 'validated_contexts': [{'context': 'CCA', 'count': 3, 'percent': 30.0}]},
        {'sample': 'Sample1', 'barcode': 'TTTTTTTTT', 'contexts': ['CTA', 'CCA'], 'context_counts': {'CTA': 7, 'CCA': 3}, 'total_count': 10, 'R1_sequences': set(), 'R2_sequences': set(), 'validated_contexts': [{'context': 'CTA', 'count': 7, 'percent': 70.0}]},
    ])

    dist = get_context_distribution(df, 'validated_contexts')
    assert dist['CTA'] == 14
    assert dist['CCA'] == 3


def test_get_context_distribution_with_empty():
    df = pd.DataFrame([
        {'sample': 'Sample1', 'barcode': 'AAAAAAAAG', 'contexts': ['CTA', 'CCA'], 'context_counts': {'CTA': 7, 'CCA': 3}, 'total_count': 10, 'R1_sequences': set(), 'R2_sequences': set(), 'validated_contexts': []},
        {'sample': 'Sample1', 'barcode': 'AAAAAAAAG', 'contexts': ['CCA', 'CTA'], 'context_counts': {'CTA': 7, 'CCA': 3}, 'total_count': 10, 'R1_sequences': set(), 'R2_sequences': set(), 'validated_contexts': [{'context': 'CCA', 'count': 3, 'percent': 30.0}]},
        {'sample': 'Sample1', 'barcode': 'TTTTTTTTT', 'contexts': ['CTA', 'CCA'], 'context_counts': {'CTA': 7, 'CCA': 3}, 'total_count': 10, 'R1_sequences': set(), 'R2_sequences': set(), 'validated_contexts': [{'context': 'CTA', 'count': 7, 'percent': 70.0}]},
    ])

    dist = get_context_distribution(df, 'validated_contexts')
    assert dist['CTA'] == 7
    assert dist['CCA'] == 3


def test_get_validated_contexts_invalid_min_pct():
    with pytest.raises(ValueError):
        get_validated_contexts({"A": 10}, -5)


def test_compute_correction_stats_invalid_threshold():
    df = pd.DataFrame({'validated_contexts': [[]]})
    with pytest.raises(ValueError):
        compute_correction_stats(df, 'CTA', 'CCA', -10)


def test_compute_correction_stats_invalid_context_types():
    df = pd.DataFrame({'validated_contexts': [[]]})
    
    with pytest.raises(ValueError, match="reference_context and corrected_context must be strings."):
        compute_correction_stats(df, reference_context=123, corrected_context="CCA", threshold=50)

    with pytest.raises(ValueError, match="reference_context and corrected_context must be strings."):
        compute_correction_stats(df, reference_context="CTA", corrected_context=None, threshold=50)