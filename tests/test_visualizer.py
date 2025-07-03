import os
import sys
import pandas as pd
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from visualizer import plot_barcode_count_distribution, plot_context_correction_summary


def test_plot_barcode_count_distribution_single(tmp_path):
    df = pd.DataFrame({
        'barcode': ['A', 'B', 'C', 'D'],
        'total_count': [1, 2, 15, 200]
    })
    save_path = tmp_path / "barcode_dist.png"
    plot_barcode_count_distribution(df, save_path)
    assert save_path.exists()


def test_plot_context_correction_summary(tmp_path):
    correction_result = {'correction_%': 75.0, 'no_correction_%': 25.0}
    original_context = 'CTA'
    corrected_context = 'CCA'
    save_path = tmp_path / "correction_summary.png"

    plot_context_correction_summary(correction_result, original_context, corrected_context, save_path)
    assert save_path.exists()