import pandas as pd
from collections import defaultdict, Counter

### Store the validated reads from each sample in a combined df ###

def build_barcode_summary(dict_processed_reads):
    '''
    Builds a summary DataFrame from processed validated reads.
    Args:
        dict_processed_reads: Dictionary with sample names as keys and lists of processed reads as values.
    Returns:
        A DataFrame with columns: sample, barcode, contexts, context_counts, total_count, R1_sequences, R2_sequences.
    '''

    if not isinstance(dict_processed_reads, dict):
        raise TypeError("Error in read processing. Check the input data.")
    if not dict_processed_reads:
        raise ValueError("No reads were found. Check the anchor sequences, barcode length and context length parameters.")
    if not all(isinstance(v, list) for v in dict_processed_reads.values()):
        raise ValueError("All values in dict_processed_reads must be lists of reads.")

    # Count context occurrences per (sample, barcode)
    barcode_context_counts = defaultdict(lambda: defaultdict(Counter))
    total_counts = defaultdict(Counter)
    
    # Lookup for (sample, barcode) > sets of R1 and R2 sequences
    barcode_seq_sets = defaultdict(lambda: {'R1': set(), 'R2': set()})

    for sample, reads in dict_processed_reads.items():
        for read in reads:
            barcode = read['barcode']
            context = read['context']
            barcode_context_counts[sample][barcode][context] += 1
            total_counts[sample][barcode] += 1
            # Add sequences to the sets
            barcode_seq_sets[(sample, barcode)]['R1'].add(read.get('sequence_R1', ''))
            barcode_seq_sets[(sample, barcode)]['R2'].add(read.get('sequence_R2', ''))

    # Build a structured list of rows for the DataFrame, including sets of R1 and R2 sequences
    rows = []
    for sample, barcode_dict in barcode_context_counts.items():
        for barcode, context_counter in barcode_dict.items():
            total = sum(context_counter.values())
            context_list = list(context_counter.keys())
            seq_sets = barcode_seq_sets.get((sample, barcode), {'R1': set(), 'R2': set()})
            row = {
                'sample': sample,
                'barcode': barcode,
                'contexts': context_list,
                'context_counts': dict(context_counter),
                'total_count': total,
                'R1_sequences': seq_sets['R1'],
                'R2_sequences': seq_sets['R2']
            }
            rows.append(row)

    return pd.DataFrame(rows)
