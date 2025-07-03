# Plot:
# - Barcode count distribution in the sample
# - Frequency of barcode groups
# - Correction rate per context

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

### Plotting barcode count distribution ###
def plot_barcode_count_distribution(df_barcode_summary, save_path):
    '''
    Plot distribution of barcodes across count bins for a single sample.
    X-axis = count bins, Y-axis = number of barcodes.
    Args:
        df_barcode_summary: DataFrame with columns: sample, barcode, contexts, context_counts, total_count, R1_sequences, R2_sequences.
        save_path: Path to save the plot image.
    '''
    if df_barcode_summary.empty:
        print("No validated reads were found. Exiting.")
        return
    if not isinstance(df_barcode_summary, pd.DataFrame):
        raise TypeError("Expected a DataFrame for barcode summary.")
     
    df = df_barcode_summary.copy()

    bins = [0, 1, 2, 10, 21, 51, 101, 501, 1001, float('inf')]
    labels = ['1', '2', '3–10', '11–21', '22–51', '52–101', '102–501', '502–1001', '>1001']

    df['count_bin'] = pd.cut(df['total_count'], bins=bins, labels=labels, right=True)

    # Count how many barcodes fall into each bin
    bin_counts = df['count_bin'].value_counts(sort=False)

    print("Barcode count per bin:")
    print(bin_counts)

    # Plot
    ax = bin_counts.plot(kind='bar', figsize=(10, 5), rot=45)
    plt.title("Barcode Count Distribution")
    plt.xlabel("Read Count Bin")
    plt.ylabel("Number of Barcodes")
    plt.tight_layout()
    
    # Save
    plt.savefig(save_path)
    plt.close()


### Plotting context correction summary ###
def plot_context_correction_summary(correction_result, original_context, corrected_context, save_path):
    """
    Plot and save a bar chart showing corrected vs. uncorrected context percentages.
    Args:
        correction_result: Dictionary with correction statistics.
        original_context: The context that is considered as no correction.
        corrected_context: The context that is considered as correction.
        save_path: Path to save the plot image.
    """
    if not isinstance(correction_result, dict):
        raise TypeError("Expected a dictionary for correction_result.")
    if not isinstance(original_context, str) or not isinstance(corrected_context, str):
        raise ValueError("reference_context and corrected_context must be strings.")
    
    # Prepare data as a list of tuples with string labels
    data = [
        (f'Corrected\n({corrected_context})', float(correction_result.get('correction_%', 0))),
        (f'Not Corrected\n({original_context})', float(correction_result.get('no_correction_%', 0)))
    ]
    plot_df = pd.DataFrame(data, columns=['Event Type', 'Percentage'])

    # Plot
    plt.figure(figsize=(8, 5))
    sns.barplot(data=plot_df, x='Event Type', y='Percentage', hue='Event Type', palette='Set2', legend=False)
    plt.title("Context Correction Summary")
    plt.ylabel("Percentage of Barcodes")
    plt.xlabel("Correction Category")
    plt.ylim(0, 100)
    plt.tight_layout()

    # Save
    plt.savefig(save_path)
    plt.close()