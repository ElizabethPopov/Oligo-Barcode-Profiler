### Function to filter high-confidence contexts ###
def get_validated_contexts(context_counts, min_pct):
    '''
    Identifies which contexts appear frequently enough (>= min_pct of the reads),
    and calculates percent out of the total validated contexts.
    
    Input:
        - context_counts: a dictionary like {'CCA': 120, 'GCA': 5}
        - min_pct: the minimum percent cutoff (default is 2%) for filtering noise
        - This can be applied on each row in a selected df.

    Returns:
    A list of dicts (one per validated context) that contains:
        - context sequence
        - % of context out of all the validated context counts
        - context counts with a specific barcode
      {
        'context': 'CCA',
        'percent: 87.5,
        'count': 320
      }
    '''
    if not isinstance(min_pct, (int, float)) or not (0 <= min_pct <= 100):
        raise ValueError("min_pct must be a number between 0 and 100.")
    if not isinstance(context_counts, dict):
        raise ValueError("context_counts must be a dictionary.")
    if not context_counts:
        return []

    total_contexts = sum(context_counts.values())
    if total_contexts == 0:
        return []

    # Filter by min_pct using the original total
    initial_validated = []
    for context, count in context_counts.items():
        pct = (count / total_contexts) * 100
        if pct >= min_pct:
            initial_validated.append({'context': context, 'count': count})

    # Recalculate percent out of the total validated contexts
    validated_total = sum(d['count'] for d in initial_validated)
    if validated_total == 0:
        return []

    validated = []
    for d in initial_validated:
        context = d['context']
        count = d['count']
        pct_validated = (count / validated_total) * 100
        validated.append({
            'context': context,
            'percent': round(pct_validated, 2),
            'count': count
        })

    return validated

# Count the number of times each context appears across all reads in the data
def get_context_distribution(df,col):
    '''
    Returns a Counter object with the distribution of contexts in the specified column.
    The column should contain lists of contexts (e.g., ['CCA', 'GCA']) or strings that can be evaluated
    to lists (e.g., "['CCA', 'GCA']").
    '''
    from collections import Counter
    import pandas as pd

    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame.")
    
    if col not in df.columns:
        raise ValueError(f"Column '{col}' does not exist in the DataFrame.")

    context_counter = Counter()
    for row in df[col]:
        if isinstance(row, list):
            for d in row:
                context = d.get('context')
                count = d.get('count', 0)   # Get context count, default to 0 if not present
                if context is not None:
                    context_counter[context] += count
    return context_counter  


def compute_correction_stats(df, reference_context, corrected_context, threshold, validated_col='validated_contexts'):
    '''
    Computes the correction statistics for a single sample DataFrame.
    Args:
        df: DataFrame containing the sample data.
        reference_context: The context that is considered as no correction.
        corrected_context: The context that is considered as correction.
        threshold: The minimum percentage threshold to consider a context as dominant.
        validated_col: The column containing validated contexts (default is 'validated_contexts').
    Returns:
        A dictionary with the following keys:
            - 'sample': The sample name.
            - 'total_barcodes': Total number of barcodes processed (with correction/no correction events).
            - 'correction_events': Number of correction events.
            - 'no_correction_events': Number of no correction events.
            - 'correction_%': Percentage of correction events.
            - 'no_correction_%': Percentage of no correction events.
    '''
    import pandas as pd
    
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame.")
    
    if not isinstance(reference_context, str) or not isinstance(corrected_context, str):
        raise ValueError("reference_context and corrected_context must be strings.")
    
    if not isinstance(threshold, (int, float)) or not (0 <= threshold <= 100):
        raise ValueError("threshold must be a number between 0 and 100.")

    correction = 0
    no_correction = 0

    for row in df[validated_col]:
        if isinstance(row, list):
            if len(row) == 1:
                ctx = row[0]['context']
                if ctx == reference_context:
                    no_correction += 1
                elif ctx == corrected_context:
                    correction += 1
            elif len(row) > 1:
                dominant = None
                for d in row:
                    ctx = d['context']
                    pct = d.get('avg_percent', 0)
                    if ctx == reference_context and pct >= threshold:
                        no_correction += 1
                        dominant = 'no_correction'
                        break
                if dominant != 'no_correction':
                    correction += 1  # fallback to correction if reference context is not dominant

    total = correction + no_correction
    correction_pct = (correction / total * 100) if total > 0 else 0
    no_correction_pct = (no_correction / total * 100) if total > 0 else 0

    return {
        'sample': df['sample'].iloc[0] if 'sample' in df.columns else 'sample',
        'total_barcodes': total,
        'correction_events': correction,
        'no_correction_events': no_correction,
        'correction_%': round(correction_pct, 2),
        'no_correction_%': round(no_correction_pct, 2)
    }