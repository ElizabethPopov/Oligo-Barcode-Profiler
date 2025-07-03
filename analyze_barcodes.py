import argparse
import textwrap
import sys
import os
import pandas as pd
from fastq_reader import load_validated_reads
from barcode_parser import compile_pattern
from summary_builder import build_barcode_summary
from mutation_analyzer import get_validated_contexts, compute_correction_stats, get_context_distribution
from visualizer import plot_barcode_count_distribution, plot_context_correction_summary


def parse_args():
    parser = argparse.ArgumentParser(
        prog='Oligo-Barcode-Profiler',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            Analyze barcoded oligonucleotide libraries from paired-end Illumina sequencing data
            -----------------------------------------------------------------------------------
                        The program will search for reads containing the pattern:
                            [anchor1][barcode][anchor2][context][anchor3]
            ''')
    )

    parser.add_argument('--r1', type=str, required=True, help='FASTQ file of forward reads')
    parser.add_argument('--r2', type=str, required=True, help='FASTQ file of reverse reads')
    parser.add_argument('--min-pct', type=float, default=40, help='Minimum percent for context validation. Default = 40')
    parser.add_argument('--context', type=str, required=True, help='Original context sequence (e.g., CTA)')
    parser.add_argument('--corrected-context', type=str, required=True, help='Corrected version of the context (e.g., CCA)')
    parser.add_argument('--anchor1', type=str, required=True, help='Anchor before the barcode')
    parser.add_argument('--anchor2', type=str, required=True, help='Anchor between barcode and context')
    parser.add_argument('--anchor3', type=str, required=True, help='Anchor after the context')
    parser.add_argument('--anch1-mm', type=int, default=2, help='Mismatches allowed in anchor1. Default = 2')
    parser.add_argument('--anch2-mm', type=int, default=1, help='Mismatches allowed in anchor2. Default = 1')
    parser.add_argument('--anch3-mm', type=int, default=2, help='Mismatches allowed in anchor3. Default = 2')
    parser.add_argument('--barcode-length', type=int, default=9, help='Length of the barcode. Default = 9')
    parser.add_argument('--context-length', type=int, default=3, help='Length of the context. Default = 3')
    parser.add_argument('--output-dir', type=str, default='output', help='Directory to store output files (default: output)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Strip whitespace from all DNA sequence arguments
    for attr in ['anchor1', 'anchor2', 'anchor3', 'context', 'corrected_context']:
        value = getattr(args, attr)
        if isinstance(value, str):
            setattr(args, attr, value.strip())

    # Collect all validation errors
    errors = []

    # Validate arguments before proceeding
    if not (0 <= args.min_pct <= 100):
        errors.append(f"--min-pct must be between 0 and 100. Got: {args.min_pct}")
    if len(args.context) != args.context_length:
        errors.append(f"--context must match --context-length. Context length: {args.context_length}. Got: '{args.context}'")
    if len(args.corrected_context) != args.context_length:
        errors.append(f"--corrected-context must match --context-length. Context length: {args.context_length}. Got: '{args.corrected_context}'")
    if not (isinstance(args.barcode_length, int) and args.barcode_length > 0):
        errors.append(f"--barcode-length must be a positive integer. Got: {args.barcode_length}")
    if not (isinstance(args.context_length, int) and args.context_length > 0):
        errors.append(f"--context-length must be a positive integer. Got: {args.context_length}")
    
    # Validate anchors and contexts
    dna_params = [
        ('--anchor1', args.anchor1),
        ('--anchor2', args.anchor2),
        ('--anchor3', args.anchor3),
        ('--context', args.context),
        ('--corrected-context', args.corrected_context)
    ]
    
    for name, seq in dna_params:
        print(f"Validating {name}: '{seq}'")  # Add this line
        if not (isinstance(seq, str) and len(seq) > 0):
            errors.append(f"{name} must be a non-empty string.")
        if not all(c in 'ACGT' for c in seq):
            errors.append(f"{name} can only contain A, C, G, T characters. Got: '{seq}'")

    # Validate all mm arguments in a loop
    int_params = [
        ('--anch1-mm', args.anch1_mm),
        ('--anch2-mm', args.anch2_mm),
        ('--anch3-mm', args.anch3_mm)
    ]
    for name, value in int_params:
        if not (isinstance(value, int) and value >= 0):
            errors.append(f"{name} must be a non-negative integer. Got: {value}")
    
    # If any errors, print them all and exit
    if errors:
        raise ValueError("Argument validation failed.\n" + "\n".join(errors))

    # Print the header if verbose mode is enabled
    if args.verbose:
        print("\n       === Oligo-Barcode-Profiler ===")
        print("              Pattern searched:")
        print("  [anchor1][barcode][anchor2][context][anchor3]")
        print('-' * 60)

    os.makedirs(args.output_dir, exist_ok=True)

    # Compile pattern
    compiled_pattern = compile_pattern(
        args.anchor1, args.anchor2, args.anchor3,
        args.barcode_length, args.context_length,
        args.anch1_mm, args.anch2_mm, args.anch3_mm
    )

    # Check uf I want this
    if args.verbose:
        print(f"Compiled pattern: {compiled_pattern.pattern}")
        print(f"Minimum context validation percentage: {args.min_pct}%")
        print(f"Context: {args.context}, Corrected Context: {args.corrected_context}")
        print(f"Barcode length: {args.barcode_length}, Context length: {args.context_length}")
        print(f"Anchors: [{args.anchor1}], [{args.anchor2}], [{args.anchor3}]")
        print(f"Output directory: {args.output_dir}")
        print('-' * 60)

    # Load and validate reads
    sample_name, reads = load_validated_reads(args.r1, args.r2, compiled_pattern)
    dict_processed_reads = {sample_name: reads}
    sample_name = list(dict_processed_reads.keys())[0]

    # Create summary
    df_barcode_summary = build_barcode_summary(dict_processed_reads)
    
    # Check if the summary is empty and exist if no reads were validated
    if df_barcode_summary.empty:
        print("No validated reads were found. Exiting.")
        sys.exit(1)

    df_barcode_summary.to_csv(os.path.join(args.output_dir, f"{sample_name}_barcode_count_summary_unfiltered.csv"), index=False)

    if args.verbose:
        print(f"Unique barcodes in {sample_name}: {df_barcode_summary['barcode'].nunique()}")
        print(f"Total barcode counts: {df_barcode_summary['total_count'].sum()}")
        print('-' * 60)

    # Plot barcode distribution
    plot_barcode_count_distribution(
        df_barcode_summary,
        save_path=os.path.join(args.output_dir, f"{sample_name}_barcode_count_distribution.png")
    )

    # Validate contexts
    df_barcode_summary['validated_contexts'] = df_barcode_summary.apply(
        lambda row: get_validated_contexts(row['context_counts'], min_pct=args.min_pct), axis=1
    )
    
    df_barcode_summary = df_barcode_summary[df_barcode_summary['validated_contexts'].apply(len) > 0]
    
    # Check if the summary is empty after context validation
    if df_barcode_summary.empty:
        print("No validated contexts were found. Exiting.")
        sys.exit(1)
    
    df_barcode_summary.to_csv(os.path.join(args.output_dir, f"{sample_name}_barcode_count_summary_validated.csv"), index=False)

    # Context distribution
    context_dist = get_context_distribution(df_barcode_summary, 'validated_contexts')
    context_pct = pd.Series(context_dist) / sum(context_dist.values()) * 100
    context_pct.sort_values(ascending=False).to_csv(
        os.path.join(args.output_dir, f"{sample_name}_context_distribution_percent.csv")
    )

    if args.verbose:
        print("\nTop 10 context distributions across all reads (%):")
        print(context_pct.sort_values(ascending=False).head(10))
        print('-' * 60)

    # Correction summary
    correction_result = compute_correction_stats(
        df_barcode_summary,
        reference_context=args.context,
        corrected_context=args.corrected_context,
        threshold=args.min_pct
    )

    pd.DataFrame([correction_result]).to_csv(
        os.path.join(args.output_dir, f"{sample_name}_correction_summary.csv"), index=False
    )

    if args.verbose:
        print("=== Correction Summary ===")
        for k, v in correction_result.items():
            print(f"{k}: {v}")
        print('-' * 60)

    plot_context_correction_summary(
        correction_result,
        original_context=args.context,
        corrected_context=args.corrected_context,
        save_path=os.path.join(args.output_dir, f"{sample_name}_correction_summary.png")
    )

if __name__ == "__main__":
    import time
    import sys
    start_time = time.time()

    try:
        main()
    except ValueError as e1:
        print(f"\nError: {e1}")
        sys.exit(1)
    except FileNotFoundError as e2:
        print(f"\nFile Error: {e2}")
        sys.exit(1)
    except Exception as e3:
        print(f"\nUnexpected Error: {e3}")
        sys.exit(1)

    end_time = time.time()
    print(f"\nTotal runtime: {end_time - start_time:.2f} seconds")