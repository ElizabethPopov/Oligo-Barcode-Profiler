import os
import sys
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from summary_builder import build_barcode_summary

def test_build_barcode_summary():
    data = {
        "sample1": [
            {"barcode": "AAAAAAAAG", "context": "CTA"},
            {"barcode": "AAAAAAAAG", "context": "CCA"},
            {"barcode": "TTTTTTTTT", "context": "CTA"},
            {"barcode": "TTTTTTTTT", "context": "CTA"}
        ]
    }
    df = build_barcode_summary(data)
    assert len(df) == 2  # two unique barcodes
    assert df['total_count'].sum() == 4 # total counts of all barcodes
    assert 'context_counts' in df.columns
    
    # Check context counts for each barcode
    row1 = df[df['barcode'] == 'AAAAAAAAG'].iloc[0]
    row2 = df[df['barcode'] == 'TTTTTTTTT'].iloc[0]
    assert row1['context_counts']['CTA'] == 1
    assert row1['context_counts']['CCA'] == 1
    assert row2['context_counts']['CTA'] == 2

    # Check the total count for each barcode
    assert row1['total_count'] == 2  # 1 CTA + 1 CCA for this barcode
    assert row2['total_count'] == 2  # 2 CTA for this barcode


def test_build_barcode_summary_invalid_input():
    with pytest.raises(TypeError):
        build_barcode_summary("Error in read processing")

def test_build_barcode_summary_missing_fields():
    data = {}  # Empty data
    with pytest.raises(ValueError):
        build_barcode_summary(data)
